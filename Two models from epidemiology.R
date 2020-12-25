#---------------------------------------------------------------------------
#   Solving the Kermack-McKendrick Model numerically (by Euler's method)
#---------------------------------------------------------------------------

Kerm.Mck<-function(n, y0, R, gamma = 1, ndays = 365, h = 0.2, plot.option = 0){
	x0 = n - y0
	lambda = R * gamma / x0
	days = seq(1, ndays, by = h)
	d = length(days)
	x = rep(x0, d); y = rep(y0, d); z = rep(0, d)
	for(i in 1:(d-1)){
		x[i+1] = x[i] - h * lambda*x[i]*y[i]
		z[i+1] = z[i] + h * gamma*y[i]
		y[i+1] = n - x[i+1] - z[i+1]
	}
	x = x/n; y = y/n; z = z/n
	if(plot.option == 0){
		plot(days, x, ylim = c(0, x0/n), xlab = "time (t)", ylab = "x(t), y(t), z(t) divided by n", type = "l", col = "green3", lwd = 2)
		lines(days, y, type = "l", col = "red", lwd = 2)
		lines(days, z, type = "l", col = "black", lwd = 2)
		if(min(abs(x - z)) > 0.3)		
			legend("topright", legend = c("x(t)", "y(t)", "z(t)"), col = c("green3", "red", "black"), lwd = 2)
		else
			legend("right", legend = c("x(t)", "y(t)", "z(t)"), col = c("green3", "red", "black"), lwd = 2)
		title(main = paste("R_0 = ", R, ",  y(0)/n =", y0/n, ",  gamma =", gamma), line = 1)
	}
	return(list(x = x, y = y, z = z))
}

par(mfrow = c(2, 3), family = "System Font")
Kerm.Mck(n = 1000, y0 = 50, R = 0.8, ndays = 20)
Kerm.Mck(n = 1000, y0 = 50, R = 1, ndays = 20)
Kerm.Mck(n = 1000, y0 = 50, R = 1.2, ndays = 20)
Kerm.Mck(n = 1000, y0 = 200, R = 0.8, ndays = 20)
Kerm.Mck(n = 1000, y0 = 200, R = 1, ndays = 20)
Kerm.Mck(n = 1000, y0 = 200, R = 1.2, ndays = 20)

Phase.Trans.KM<-function(n = 10^5, y0 = 200, R = 1.2, ndays = 100, h = 0.2){
	x0 = n - y0
	theta = R/x0
	days = seq(1, ndays, by = h)
	pts = seq(0, n, by = 0.1)
	w = pts[which.min(abs(n - x0 * exp(- theta * pts) - pts))]
	K = Kerm.Mck(n, y0, R, ndays = ndays, h = h, plot.option = 1)
	tstar = days[which.min(abs(K$x*n - 1/theta))]
	return(c(w = w, ymax = max(K$y), tstar = tstar))
}

Phase.Trans.KM()

n = 10^5; y0 = 200
Rseq = seq(0.8, 1.2, by = 0.01)
List = sapply(Rseq, function(r) Phase.Trans.KM(n, y0, R = r))
par(mfrow = c(1, 3), family = "System Font")
plot(unlist(List[1,])/n ~ Rseq, type = "l", lwd = 2, xlab = "Value of R_0", ylab = "Value of w/n = z(infty)/n", main = paste("How z(infty)/n changes with R_0\n where n =", n,", y(0) =", y0))
plot(unlist(List[2,]) ~ Rseq, type = "l", lwd = 2, xlab = "Value of R_0", ylab = "Value of y.max/n", main = paste("How y.max/n changes with R_0\n where n =", n,", y(0) =", y0))
plot(unlist(List[3,]) ~ Rseq, type = "l", lwd = 2, xlab = "Value of R_0", ylab = "Value of t.star", main = paste("How t.star changes with R_0\n where n =", n,", y(0) =", y0))

#---------------------------------------------------------------------------
#                      Simulating the Reed-Frost Model                     
#---------------------------------------------------------------------------

Reed.Frost<-function(n, y0, R, ndays = 365, plot.option = 0){
	x0 = n - y0
	p = R / x0
	days = 1:ndays
	x = rep(x0, ndays); y = rep(y0, ndays); z = rep(0, ndays)
	for(i in 1:(ndays-1)){
		y[i+1] = rbinom(1, x[i], 1 - (1 - p)^y[i])
		x[i+1] = x[i] - y[i+1]
		z[i+1] = z[i] + y[i]
	}
	x = x/n; y = y/n; z = z/n
	if(plot.option == 0){
		plot(days, x, ylim = c(0, 1), xlab = "time (t)", ylab = "x(t), y(t), z(t) divided by n", type = "l", col = "green3", lwd = 2)
		lines(days, y, type = "l", col = "red", lwd = 2)
		lines(days, z, type = "l", col = "black", lwd = 2)
		if(min(abs(x - z)) > 0.3)
			legend("right", legend = c("x(t)", "y(t)", "z(t)"), col = c("green3", "red", "black"), lwd = 2)
		else #so that legend does not collide with the curves
			legend("topright", legend = c("x(t)", "y(t)", "z(t)"), col = c("green3", "red", "black"), lwd = 2)
		title(main = paste("R_0 = ", R, ",  y(0)/n =", y0/n, ",  p =", round(p, 5)), line = 1)
	}
	return(list(x = x, y = y, z = z))
}

par(mfrow = c(2, 3), family = "System Font")
Reed.Frost(n = 1000, y0 = 50, R = 0.8, ndays = 20)
Reed.Frost(n = 1000, y0 = 50, R = 1, ndays = 20)
Reed.Frost(n = 1000, y0 = 50, R = 1.2, ndays = 20)
Reed.Frost(n = 1000, y0 = 200, R = 0.8, ndays = 20)
Reed.Frost(n = 1000, y0 = 200, R = 1, ndays = 20)
Reed.Frost(n = 1000, y0 = 200, R = 1.2, ndays = 20)

Phase.Trans.RF<-function(n = 10^5, y0 = 200, R = 1.2, ndays = 100, opn = 1){
	x0 = n - y0
	theta = R/x0
	days = seq(1, ndays, by = 1)
	pts = seq(0, n, by = 0.1)
	w = pts[which.min(abs(n - x0 * exp(- theta * pts) - pts))]
	K = Reed.Frost(n, y0, R, ndays = ndays, plot.option = opn)
	if(opn == 0) abline(h = w/n, lty = 2, col = "blue")
	tstar = days[which.min(abs(K$x*n - 1/theta))]
	return(c(w = w, ymax = max(K$y), tstar = tstar))
}

Phase.Trans.RF(opn = 0)

n = 10^5; y0 = 200
Rseq = seq(0.8, 1.2, by = 0.01)
List = sapply(Rseq, function(r) Phase.Trans.RF(n, y0, R = r))
par(mfrow = c(1, 3), family = "System Font")
plot(unlist(List[1,])/n ~ Rseq, type = "l", lwd = 2, xlab = "Value of R_0", ylab = "Value of w/n = z(infty)/n", main = paste("How z(infty)/n changes with R_0\n where n =", n,", y(0) =", y0))
plot(unlist(List[2,]) ~ Rseq, type = "l", lwd = 2, xlab = "Value of R_0", ylab = "Value of y.max/n", main = paste("How y.max/n changes with R_0\n where n =", n,", y(0) =", y0))
plot(unlist(List[3,]) ~ Rseq, type = "l", lwd = 2, xlab = "Value of R_0", ylab = "Value of t.star", main = paste("How t.star changes with R_0\n where n =", n,", y(0) =", y0))
