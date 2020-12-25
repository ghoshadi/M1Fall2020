# Following function simulates the Poisson branching process with Poi(lambda) as the progeny distribution.
GW.pois<-function(lambda){
	W = Z = 1
	while(min(Z) > 0)
		Z = c(Z, W <- rpois(1, W*lambda))
	return(Z)
}

# Below we study the extinction time (T_e) of the Poisson branching process with mean < 1.

N = 1000 # We generated the same process N no. of times to understand the stochastic behaviour of T_e
lambda = seq(0.1, 1, by = 0.02)
L = length(lambda)
T_e.mat = matrix(nrow = N, ncol = L) # matrix that will record the extinction times in each simulation
for(l in 1:L)
	for(k in 1:N) 
		T_e.mat[k,l] = length(GW.pois(lambda[l]))
T_e.mean = apply(T_e.mat, 2, mean)

plot(lambda, T_e.mean, pch = 1, ylab = "Avergae time of extinction", col = "blue")
lines(lambda, T_e.mean, lty = 1, col = "blue")
title(main = "How the average time of extinction E[T_e] of GWBP varies with lambda\n (where the progeny distribution is Poisson with mean lambda)", line = 1)

par(mfrow = c(3,3))
for(l in 1:9) plot(table(T_e.mat[ ,5*l-4]), ylab = "Freq", main = paste("lambda =", lambda[5*l-4]))
