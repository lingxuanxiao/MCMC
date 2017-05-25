#Program MCMC
#Markov Chain Monte-Carlo Method
#Author: Xiaohua_Ye
#======================================================================================
#Public
#1. The probability density function of cauchy distribution
#1.1 Cauchy Distribution:
Cauchy = function(theta, mu = 0, gam = 1){return(1/(pi*gam*(1+((theta-mu)/gam)^2)))}
#1.2 Exponential Distribution
Exponential = function(theta, lam = 1){return(ifelse(theta > 0, lam*exp(-lam*theta), 0))}
#1.3 2 Dim Normal Distribution
twoDimNormal = function(theta1, theta2, mu1 = 0, mu2 = 0, sigma1 = 1, sigma2 = 1, rou = 0){return(exp(-((theta1-mu1)^2/sigma1^2+(theta2-mu2)^2/sigma2^2-2*rou*(theta1-mu1)*(theta2-mu2)/sigma1/sigma2)/(2*sqrt(1-rou*rou)))/(2*pi*sigma1*sigma2*sqrt(1-rou*rou)))}
#1.4 Conditional Distribution of 2 Dim Normal Distribution
conTwoDimNormal = function(theta1, mu1 = 0, mu2 = 0, sigma1 = 1, sigma2 = 1, rou = 0){return(rnorm(1, rou*sigma2/sigma1*(theta1-mu1)+mu2, sigma2*sqrt(1-rou*rou)))}
#2. Initialization
#2.1 Define public variable and assign it
Time = 20000; Theta = numeric(20001); Seed = 8888; Theta[1] = 1
#2.2 Initialization function
Init = function(T = 20000, Begin = 1, seed = 8888){Time <<- T; Theta <<- numeric(T+1); Seed <<- seed; Theta[1] <<- Begin}
#2.3 drawing colors
jet.colors = colorRampPalette(c('#000080', '#0000FF', '#0080FF', '#00FFFF', '#80FF80', '#FFFF00', '#FF8000', '#FF0000', '#800000'))
#======================================================================================
#MCMC
#Method: Metropolis
#Proposal Distribution: Normal Distribution
#Target Distribution: Cauchy Distribution
Metropolis = function(aptoticRandomSeed = TRUE){
	if(aptoticRandomSeed){set.seed(Seed)}
	for(i in 1:Time+1){
		thetaStar = Theta[i-1] + rnorm(1)
		u = runif(1)
		Theta[i] <<- ifelse(u <= Cauchy(thetaStar)/Cauchy(Theta[i-1]), thetaStar, Theta[i-1])
	}
}
#--------------------------------------------------------------------------------------
#Plot and analysis
Init()
Metropolis()
pdf('C:\\Users\\Feixiao_L\\Desktop\\pdf\\Metropolis.pdf')
plot(1:10000, round(Theta[10001:20000+1], 2), 'l', ylim = c(-30, 30), xlab = 'Index', ylab = 'Random number from MC', main = 'Trace of MCMC')
plot(1:10000, abs(fft(round(Theta[10001:20000+1], 2))), 'l', xlab = 'Cycle', ylab = 'Amplitude', main = 'Spectrum of MCMC', ylim = c(0, 1000))
hist(Theta[10001:20000+1], c(-30:31-0.5), freq = FALSE, main = 'Histogram of MCMC', xlab = 'Value')
points(-3000:3000/100, Cauchy(-3000:3000/100), 'l', lty = 2, col = 'dark blue')
plot(density(Theta[10001:20000+1]), 'l', col = 'dark red', xlim = c(-10, 10), xlab = 'Value', ylab = 'Density', main = 'Power Spectrum of MCMC')
points(-1000:1000/100, Cauchy(-1000:1000/100), 'l', lty = 2, col = 'dark blue')
dev.off()
#======================================================================================
#MCMC
#Method: Metropolis-Hastings
#Proposal Distribution: Chi-Squared Distribution
#Target Distribution: Exponential Distribution
metropolisHastings = function(aptoticRandomSeed = TRUE){
	if(aptoticRandomSeed){set.seed(Seed)}
	for(i in 1:Time+1){
		thetaStar = rchisq(1, df = Theta[i-1])
		u = runif(1)
		Theta[i] <<- ifelse(u <= Exponential(thetaStar)/Exponential(Theta[i-1])*dchisq(Theta[i-1], df = thetaStar)/dchisq(thetaStar, df = Theta[i-1]), thetaStar, Theta[i-1])
	}
}
#--------------------------------------------------------------------------------------
#Plot and analysis
Init()
metropolisHastings()
pdf('C:\\Users\\Feixiao_L\\Desktop\\pdf\\MetropolisHastings.pdf')
plot(1:10000, round(Theta[10001:20000+1], 2), 'l', ylim = c(0, 10), xlab = 'Index', ylab = 'Random number from MC', main = 'Trace of MCMC')
plot(1:10000, abs(fft(round(Theta[10001:20000+1], 2))), 'l', xlab = 'Cycle', ylab = 'Amplitude', main = 'Spectrum of MCMC', ylim = c(0, 1000))
hist(Theta[10001:20000+1], c(0:40/4), freq = FALSE, main = 'Histogram of MCMC', xlab = 'Value', ylim = c(0, 1))
points(0:1000/100, Exponential(0:1000/100), 'l', lty = 2, col = 'dark blue')
plot(density(Theta[10001:20000+1]), 'l', col = 'dark red', xlim = c(-1, 8), ylim = c(0, 1), xlab = 'Value', ylab = 'Density', main = 'Power Spectrum of MCMC')
points(-100:800/100, Exponential(-100:800/100), 'l', lty = 2, col = 'dark blue')
dev.off()
#======================================================================================
#MCMC
#Method: Blockwise Metropolis-Hastings
#Proposal Distribution: 2 Dim Uniform Distribution
#Target Distribution: 2 Dim Normal Distribution, mu1 = 3, mu2 = 0, sigma1 = 2, sigma2 = 1, rou = 0.7
blockMetroHastings = function(aptoticRandomSeed = TRUE){
	if(aptoticRandomSeed){set.seed(Seed)}
	#initialization Theta
	Theta <<- matrix(0, 2, Time+1)
	for(i in 1:Time+1){
		thetaStar = t(Theta[, i-1]) + runif(2, -1, 1)
		u = runif(1)
		alpha = twoDimNormal(thetaStar[1], thetaStar[2], mu1 = 3, mu2 = 0, sigma1 = 2, sigma2 = 1, rou = 0.7)/twoDimNormal(Theta[1, i-1], Theta[2, i-1], mu1 = 3, mu2 = 0, sigma1 = 2, sigma2 = 1, rou = 0.7)
		Theta[1, i] <<- ifelse(u <= alpha, thetaStar[1], Theta[1, i-1])
		Theta[2, i] <<- ifelse(u <= alpha, thetaStar[2], Theta[2, i-1])
	}
}
#--------------------------------------------------------------------------------------
#Plot and analysis
Init()
blockMetroHastings()
pdf('C:\\Users\\Feixiao_L\\Desktop\\pdf\\BlockwiseMH.pdf')
plot(1:10000, round(Theta[1, 10001:20000+1], 2), 'l', ylim = c(-6, 12), xlab = 'Index', ylab = 'Random number from MC', main = 'Marginal Trace of MCMC (x axis)')
plot(1:10000, round(Theta[2, 10001:20000+1], 2), 'l', ylim = c(-5, 5), xlab = 'Index', ylab = 'Random number from MC', main = 'Marginal Trace of MCMC (y axis)')
plot(round(Theta[1, 10001:20000+1], 2), round(Theta[2, 10001:20000+1], 2), 'l', xlim = c(-6, 12), ylim = c(-5, 5), xlab = 'x axis', ylab = 'y axis', main = 'Trace of MCMC')
plot(1:10000, abs(fft(round(Theta[1, 10001:20000+1], 2))), 'l', xlab = 'Cycle', ylab = 'Amplitude', main = 'Marginal Spectrum of MCMC (x axis)', ylim = c(0, 1000))
plot(1:10000, abs(fft(round(Theta[2, 10001:20000+1], 2))), 'l', xlab = 'Cycle', ylab = 'Amplitude', main = 'Marginal Spectrum of MCMC (y axis)', ylim = c(0, 1000))
hist(Theta[1, 10001:20000+1], 30, freq = FALSE, main = 'Marginal Histogram of MCMC (x axis)', xlab = 'Value', ylim = c(0, 0.2))
points(-600:1200/100, dnorm(-600:1200/100, 3, 2), 'l', lty = 2, col = 'dark blue')
hist(Theta[2, 10001:20000+1], 30, freq = FALSE, main = 'Marginal Histogram of MCMC (y axis)', xlab = 'Value', ylim = c(0, 0.4))
points(-500:500/100, dnorm(-500:500/100, 0, 1), 'l', lty = 2, col = 'dark blue')
plot(density(Theta[1, 10001:20000+1]), 'l', col = 'dark red', xlim = c(-6, 12), ylim = c(0, 0.2), xlab = 'Value', ylab = 'Density', main = 'Marginal Power Spectrum of MCMC (x axis)')
points(-600:1200/100, dnorm(-600:1200/100, 3, 2), 'l', lty = 2, col = 'dark blue')
plot(density(Theta[2, 10001:20000+1]), 'l', col = 'dark red', xlim = c(-5, 5), ylim = c(0, 0.4), xlab = 'Value', ylab = 'Density', main = 'Marginal Power Spectrum of MCMC (y axis)')
points(-500:500/100, dnorm(-500:500/100, 0, 1), 'l', lty = 2, col = 'dark blue')
#stat = matrix(0, 181, 101)
#for(i in 1:181){for(j in 1:101){stat[i, j] = twoDimNormal((i-61)/10, (j-51)/10, mu1 = 3, mu2 = 0, sigma1 = 2, sigma2 = 1, rou = 0.7)}}
#filled.contour(-60:120/10, -50:50/10, stat, color = jet.colors, nlevels = 100, xlab = 'x axis', ylab = 'y axis', main = 'Distribution of MCMC random number', key.title = title('Density'))
#filled.contour(-60:120/10, -50:50/10, matrix(0, 181, 101), color = colorRampPalette(c('#FFFFFF', '#FFFFFF')), nlevels = 100, plot.axes = {axis(1); axis(2); points(Theta[1, 10001:20000+1], Theta[2, 10001:20000+1], pch = 20, cex = 0.01)}, xlab = 'x axis', ylab = 'y axis', main = 'Distribution of MCMC random number', key.title = title('Density'))
dev.off()
#======================================================================================
#Method: Componentwise Metropolis-Hastings
#Proposal Distribution: Uniform Distribution
#Target Distribution: 2 Dim Normal Distribution, mu1 = 3, mu2 = 0, sigma1 = 2, sigma2 = 1, rou = 0.7
componentMetroHastings = function(aptoticRandomSeed = TRUE){
	if(aptoticRandomSeed){set.seed(Seed)}
	#initialization Theta
	Theta <<- matrix(0, 2, Time+1)
	for(i in 1:Time+1){
		thetaStar1 = Theta[1, i-1] + runif(1, -1, 1)
		thetaStar2 = Theta[2, i-1] + runif(1, -1, 1)
		u1 = runif(1)
		u2 = runif(1)
		alpha1 = twoDimNormal(thetaStar1, Theta[2, i-1], mu1 = 3, mu2 = 0, sigma1 = 2, sigma2 = 1, rou = 0.7)/twoDimNormal(Theta[1, i-1], Theta[2, i-1], mu1 = 3, mu2 = 0, sigma1 = 2, sigma2 = 1, rou = 0.7)
		alpha2 = twoDimNormal(Theta[1, i-1], thetaStar2, mu1 = 3, mu2 = 0, sigma1 = 2, sigma2 = 1, rou = 0.7)/twoDimNormal(Theta[1, i-1], Theta[2, i-1], mu1 = 3, mu2 = 0, sigma1 = 2, sigma2 = 1, rou = 0.7)
		Theta[1, i] <<- ifelse(u1 <= alpha1, thetaStar1, Theta[1, i-1])
		Theta[2, i] <<- ifelse(u2 <= alpha2, thetaStar2, Theta[2, i-1])
	}
}
#--------------------------------------------------------------------------------------
#Plot and analysis
Init()
componentMetroHastings()
pdf('C:\\Users\\Feixiao_L\\Desktop\\pdf\\ComponentwiseMH.pdf')
plot(1:10000, round(Theta[1, 10001:20000+1], 2), 'l', ylim = c(-6, 12), xlab = 'Index', ylab = 'Random number from MC', main = 'Marginal Trace of MCMC (x axis)')
plot(1:10000, round(Theta[2, 10001:20000+1], 2), 'l', ylim = c(-5, 5), xlab = 'Index', ylab = 'Random number from MC', main = 'Marginal Trace of MCMC (y axis)')
plot(round(Theta[1, 10001:20000+1], 2), round(Theta[2, 10001:20000+1], 2), 'l', xlim = c(-6, 12), ylim = c(-5, 5), xlab = 'x axis', ylab = 'y axis', main = 'Trace of MCMC')
plot(1:10000, abs(fft(round(Theta[1, 10001:20000+1], 2))), 'l', xlab = 'Cycle', ylab = 'Amplitude', main = 'Marginal Spectrum of MCMC (x axis)', ylim = c(0, 1000))
plot(1:10000, abs(fft(round(Theta[2, 10001:20000+1], 2))), 'l', xlab = 'Cycle', ylab = 'Amplitude', main = 'Marginal Spectrum of MCMC (y axis)', ylim = c(0, 1000))
hist(Theta[1, 10001:20000+1], 30, freq = FALSE, main = 'Marginal Histogram of MCMC (x axis)', xlab = 'Value', ylim = c(0, 0.2))
points(-600:1200/100, dnorm(-600:1200/100, 3, 2), 'l', lty = 2, col = 'dark blue')
hist(Theta[2, 10001:20000+1], 30, freq = FALSE, main = 'Marginal Histogram of MCMC (y axis)', xlab = 'Value', ylim = c(0, 0.4))
points(-500:500/100, dnorm(-500:500/100, 0, 1), 'l', lty = 2, col = 'dark blue')
plot(density(Theta[1, 10001:20000+1]), 'l', col = 'dark red', xlim = c(-6, 12), ylim = c(0, 0.2), xlab = 'Value', ylab = 'Density', main = 'Marginal Power Spectrum of MCMC (x axis)')
points(-600:1200/100, dnorm(-600:1200/100, 3, 2), 'l', lty = 2, col = 'dark blue')
plot(density(Theta[2, 10001:20000+1]), 'l', col = 'dark red', xlim = c(-5, 5), ylim = c(0, 0.4), xlab = 'Value', ylab = 'Density', main = 'Marginal Power Spectrum of MCMC (y axis)')
points(-500:500/100, dnorm(-500:500/100, 0, 1), 'l', lty = 2, col = 'dark blue')
#stat = matrix(0, 181, 101)
#for(i in 1:181){for(j in 1:101){stat[i, j] = twoDimNormal((i-61)/10, (j-51)/10, mu1 = 3, mu2 = 0, sigma1 = 2, sigma2 = 1, rou = 0.7)}}
#filled.contour(-60:120/10, -50:50/10, stat, color = jet.colors, nlevels = 100, xlab = 'x axis', ylab = 'y axis', main = 'Distribution of MCMC random number', key.title = title('Density'))
#filled.contour(-60:120/10, -50:50/10, matrix(0, 181, 101), color = colorRampPalette(c('#FFFFFF', '#FFFFFF')), nlevels = 100, plot.axes = {axis(1); axis(2); points(Theta[1, 10001:20000+1], Theta[2, 10001:20000+1], pch = 20, cex = 0.01)}, xlab = 'x axis', ylab = 'y axis', main = 'Distribution of MCMC random number', key.title = title('Density'))
dev.off()
#======================================================================================
#Method: Gibbs Sampling
#Proposal Distribution: Conditional Distribution of Target Distribution
#Target Distribution: 2 Dim Normal Distribution, mu1 = 3, mu2 = 0, sigma1 = 2, sigma2 = 1, rou = 0.7
Gibbs = function(aptoticRandomSeed = TRUE){
	if(aptoticRandomSeed){set.seed(Seed)}
	#initialization Theta
	Theta <<- matrix(0, 2, Time+1)
	for(i in 1:Time+1){
		Theta[2, i] <<- conTwoDimNormal(Theta[1, i-1], mu1 = 3, mu2 = 0, sigma1 = 2, sigma2 = 1, rou = 0.7)
		Theta[1, i] <<- conTwoDimNormal(Theta[2, i], mu2 = 3, mu1 = 0, sigma2 = 2, sigma1 = 1, rou = 0.7)
	}
}
#--------------------------------------------------------------------------------------
#Plot and analysis
Init()
Gibbs()
pdf('C:\\Users\\Feixiao_L\\Desktop\\pdf\\Gibbs.pdf')
plot(1:10000, round(Theta[1, 10001:20000+1], 2), 'l', ylim = c(-6, 12), xlab = 'Index', ylab = 'Random number from MC', main = 'Marginal Trace of MCMC (x axis)')
plot(1:10000, round(Theta[2, 10001:20000+1], 2), 'l', ylim = c(-5, 5), xlab = 'Index', ylab = 'Random number from MC', main = 'Marginal Trace of MCMC (y axis)')
plot(round(Theta[1, 10001:20000+1], 2), round(Theta[2, 10001:20000+1], 2), 'l', xlim = c(-6, 12), ylim = c(-5, 5), xlab = 'x axis', ylab = 'y axis', main = 'Trace of MCMC')
plot(1:10000, abs(fft(round(Theta[1, 10001:20000+1], 2))), 'l', xlab = 'Cycle', ylab = 'Amplitude', main = 'Marginal Spectrum of MCMC (x axis)', ylim = c(0, 1000))
plot(1:10000, abs(fft(round(Theta[2, 10001:20000+1], 2))), 'l', xlab = 'Cycle', ylab = 'Amplitude', main = 'Marginal Spectrum of MCMC (y axis)', ylim = c(0, 1000))
hist(Theta[1, 10001:20000+1], 30, freq = FALSE, main = 'Marginal Histogram of MCMC (x axis)', xlab = 'Value', ylim = c(0, 0.2))
points(-600:1200/100, dnorm(-600:1200/100, 3, 2), 'l', lty = 2, col = 'dark blue')
hist(Theta[2, 10001:20000+1], 30, freq = FALSE, main = 'Marginal Histogram of MCMC (y axis)', xlab = 'Value', ylim = c(0, 0.4))
points(-500:500/100, dnorm(-500:500/100, 0, 1), 'l', lty = 2, col = 'dark blue')
plot(density(Theta[1, 10001:20000+1]), 'l', col = 'dark red', xlim = c(-6, 12), ylim = c(0, 0.2), xlab = 'Value', ylab = 'Density', main = 'Marginal Power Spectrum of MCMC (x axis)')
points(-600:1200/100, dnorm(-600:1200/100, 3, 2), 'l', lty = 2, col = 'dark blue')
plot(density(Theta[2, 10001:20000+1]), 'l', col = 'dark red', xlim = c(-5, 5), ylim = c(0, 0.4), xlab = 'Value', ylab = 'Density', main = 'Marginal Power Spectrum of MCMC (y axis)')
points(-500:500/100, dnorm(-500:500/100, 0, 1), 'l', lty = 2, col = 'dark blue')
#stat = matrix(0, 181, 101)
#for(i in 1:181){for(j in 1:101){stat[i, j] = twoDimNormal((i-61)/10, (j-51)/10, mu1 = 3, mu2 = 0, sigma1 = 2, sigma2 = 1, rou = 0.7)}}
#filled.contour(-60:120/10, -50:50/10, stat, color = jet.colors, nlevels = 100, xlab = 'x axis', ylab = 'y axis', main = 'Distribution of MCMC random number', key.title = title('Density'))
#filled.contour(-60:120/10, -50:50/10, matrix(0, 181, 101), color = colorRampPalette(c('#FFFFFF', '#FFFFFF')), nlevels = 100, plot.axes = {axis(1); axis(2); points(Theta[1, 10001:20000+1], Theta[2, 10001:20000+1], pch = 20, cex = 0.01)}, xlab = 'x axis', ylab = 'y axis', main = 'Distribution of MCMC random number', key.title = title('Density'))
dev.off()