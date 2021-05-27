require(tseriesChaos);
require(colorspace)
require(viridis)

#palette(qualitative_hcl("Set 2", n=9));
palette(c("black", "blue", "red"))

embedd = tseriesChaos::embedd;
sim.cont = tseriesChaos::sim.cont;
lorenz.s0 = c(9.227021424263037729929, 15.117356048661351408668, 17.911512308502612000893);

if(length(lorenz.ts) == 2001)
	lorenz.ts = NULL;

graphics.off();
dev.new(width=0.6*12, height=0.6*8);

logistic = function(iter=40000, r=4){
	res = rep(0, iter);
	res[1] = runif(min=0, max=1, n=1);
	for(i in 2:iter)
		res[i] = r*res[i-1]*(1 - res[i-1]);
	res;
}

lorenz.real = function(N=1000, t.step=0.03){
	data = matrix(0, nrow=N, ncol=3);
	
	t.end = N * t.step;
	gambiarra = sim.cont(lorenz.syst, 0, t.end, t.step,
		 start.x=lorenz.s0, parms=c(10, 28, -8/3), obs.fun = function(x) list(x));

	for(i in 1:N){
		data[i,] = gambiarra[[i]][[1]];
	}

	data;
}

lorenz.embedded = function(N=1000, m=3, d=3, t.step=0.03){
	t.end = N * t.step;
	series = sim.cont(lorenz.syst, 0, t.end, t.step,
		 start.x=lorenz.s0, parms=c(10, 28, -8/3), obs.fun = function(x)
		 x[1]);

	data = embedd(series, m=m, d=d);
	data;
}

norm1 = function(vec){ sum(abs(vec)); }
norm2 = function(vec){ sqrt(sum(vec**2)); }
normInf = function(vec){ max(abs(vec)); }

embedd.custom = function(series, indices){
	subseries = list();
	for(l in indices){
		subseries[[length(subseries) + 1]] = series[l:length(series)];
	}

	minLength = length(series) - max(indices);
	result = NULL;
	for(i in 1:length(subseries)){
		result = cbind(result, subseries[[i]][1:minLength] );
	}

	result;
}

divergence.plot = function(data, box.x, box.y, box.z, from=20, to=3996, by=5, shuffle=FALSE, smooth=TRUE, plot=T){
	result = NULL;

	inBox = (data[,1] > box.x[1]) & (data[,1] < box.x[2]);
	inBox = inBox & (data[,2] > box.y[1]) & (data[,2] < box.y[2]);
	inBox = inBox & (data[,3] > box.z[1]) & (data[,3] < box.z[2]);
	
	if(shuffle){
		data = data[sample(1:nrow(data)),];
	}

	for(i in seq(from, to, by=by)){
		nrows = nrow(data);
		half1 = inBox[1:i];
		half2 = inBox[(i+1):(2*i)];

		freq1 = sum(half1) / length(half1);
		freq2 = sum(half2) / length(half2);

		result = rbind(result, c(i, freq1, freq2, abs(freq1 - freq2)));
	}
	if(plot){
		if(smooth){
			sm = smooth.spline(2*result[,1], result[,4], df=30);
			lines(sm, lwd=2, col="#00007755");
		} else {
			lines(2*result[,1], result[,4], lwd=2, col="#000000AA");
		}
	}
	result;
}

lorenz.single = function(data){
	N = nrow(data);
	m = 3;
	d = 3;

	xrange = range(data[,1]);
	yrange = range(data[,2]);
	zrange = range(data[,3]);

	#xrange = c(-19.11338, 19.1377);
	#yrange = c(-19.11338, 19.1377);
	#zrange = c(-19.11338, 19.1377); # Yes they're all equal

	# Real lorenz system
	#xrange = c(5.3, 14.1);
	#yrange = c(11.25, 19.3);
	#zrange = c(16.5, 24.35);

	n = 2;
	box.x = min(xrange) + (c(n-3, n) + 1) * diff(xrange) / 10;
	box.y = min(yrange) + (c(n-3, n) + 3) * diff(yrange) / 10;
	box.z = min(zrange) + (c(n-3, n) + 7) * diff(zrange) / 10;

	plot(0, 1, xlab="Sample size", ylab=expression(paste(list(histogram,divergence), phantom(aaa), epsilon, phantom(aaa))), type="l", xlim=c(0, N), ylim=c(1e-5, 1), #log="y",
			 main=paste("Divergences for the Lorenz Map (non-shuffled, rate 0.03s, m = ", m, ", d = ", d, ")."), col="#FFFFFF");

	rgl::plot3d(data);
	rgl::lines3d(t(rbind(box.x, box.y, box.z)), col=2, lwd=3);

	inBox = (data[,1] > box.x[1]) & (data[,1] < box.x[2]);
	inBox = inBox & (data[,2] > box.y[1]) & (data[,2] < box.y[2]);
	inBox = inBox & (data[,3] > box.z[1]) & (data[,3] < box.z[2]);
	
	#print(inBox[seq(1, 10000, by=97)]);
	freq = sum(inBox) / nrow(data);

	print(freq);

	result = divergence.plot(data, box.x, box.y, box.z,
							 to=N/2, by=N%/%1000, shuffle=F, smooth=F);

	x = seq(0, N * 1.1, length=1000);
	lines(x, sqrt(  (1 / (-x)) * log(0.05/2)  ), col=2)
	lines(x, sqrt(  (1 / (-x)) * log(0.001/2)  ), col=3)
	lines(x, sqrt(  (1 / (-x)) * log(1e-10/2)  ), col=4)

	legend("topright", c(expression(eta == 0.05), expression(eta == 0.001), expression(eta == 1e-10)), col=c(2, 3, 4), lwd=1)
	legend("topleft", "the probability to exceed the\ncolored lines is at most eta");

	#savePlot("histogram-rule-lorenz.png");
	result;
}

simulate = function(data, divisions=6, name="Lorenz", plot.average=F, newPlot=T, plotIdx=1){
	N = nrow(data);
	m = 3;
	d = 3;

	xrange = range(data[,1]);
	yrange = range(data[,2]);
	zrange = range(data[,3]);

	if(newPlot)
		plot(0, 1, xlab="sample size", ylab="Divergence", type="l", xlim=c(0, N), ylim=c(1e-2, 1.3), #log="y",
			 main=paste("Full histogram divergences (", name,")", sep=""), col="#FFFFFF");

	avgResult = NULL;
	for(i in seq(0, divisions-1))
	for(j in seq(0, divisions-1))
	for(k in seq(0, divisions-1)){
		#if( any(c(i, j, k) != c(0, 2, 0)) ) next;

		box.x = min(xrange) + c(i, i+1) * diff(xrange) / divisions;
		box.y = min(yrange) + c(j, j+1) * diff(yrange) / divisions;
		box.z = min(zrange) + c(k, k+1) * diff(zrange) / divisions;

		cat(i, j, k, sep="/");
		cat("\n");
		print(rbind(box.x, box.y, box.z));

		# if resampling
		#data = data[seq(1, nrow(data), by=5),]

		#rgl::plot3d(data);
		#rgl::lines3d(t(rbind(box.x, box.y, box.z)), col=2, lwd=3);

		inBox = (data[,1] > box.x[1]) & (data[,1] < box.x[2]);
		inBox = inBox & (data[,2] > box.y[1]) & (data[,2] < box.y[2]);
		inBox = inBox & (data[,3] > box.z[1]) & (data[,3] < box.z[2]);
		
		#print(inBox[seq(1, 10000, by=97)]);
		freq = sum(inBox) / nrow(data);

		#if(freq == 0)
		print(N / (divisions**3));
		print(sum(inBox));
		if(sum(inBox) < N / (divisions**2))
			next;

		if(plot.average){
			result = divergence.plot(data, box.x, box.y, box.z,
								 to=N/2, by=N%/%1000, shuffle=F, plot=F);

			if(is.null(avgResult)){
				avgResult = result[,c(1, 4)];
			} else {
				avgResult[,2] = (avgResult[,2] + result[,4]);
			}
		} else {
			result = divergence.plot(data, box.x, box.y, box.z,
								 to=N/2, by=N%/%1000, shuffle=F);
		}
	}

	if(plot.average){
		# Times two because it holds the size of one half only
		lines(2*avgResult[,1], avgResult[,2], lwd=2, col=plotIdx, lty=plotIdx);
		grid(col="#00000044");
	}

	#epsilon = function(eta, x){
	#	sqrt(
	#		 2 * log(eta/2) / (-2*x)
	#	) / 2;
	#}

	x = seq(0, N * 1.1, length=1000);
	#lines(x, epsilon(0.05, x), col=2);
	#lines(x, epsilon(0.001, x), col=3);
	#lines(x, epsilon(1e-10, x), col=4);
	#lines(x, epsilon(0.2, x), col=5);

	#legend("topright", c(expression(eta == 0.05), expression(eta == 0.001), expression(eta == 1e-10), expression(eta == 0.2)), col=c(2, 3, 4, 5), lwd=1);
	#legend("bottomleft", "the probability to exceed the\ncolored lines is at most eta");

	if(plot.average){
		return(avgResult);
	} else {
		return(result);
	}
}

data = lorenz.real(5000); # 1,500,000 uses up about 60 MB
data2 = logistic(5000); data2 = embedd(data2, m=3, d=1);
result = simulate(data, name="Logistic", divisions=30, plot.average=T);
result = simulate(data2, name="NOAA Vancouver Precipitation", divisions=30, plot.average=T, newPlot=F, plotIdx=2);

#data = read.csv("../../sunspot.csv", sep=";"); data = data[,4]; data = embedd(data, m=3, d=17);

data = read.csv("NOAA-vancouver-clean.csv"); data = data[,"TAVG"]; data = embedd.custom(data, c(1, 2, 365));
result = simulate(data, name="NOAA Vancouver Precipitation", divisions=30, plot.average=T, newPlot=F, plotIdx=3);

#data2 = read.csv("../../NOAA-vancouver-clean.csv"); data2 = data2[,"PRCP7"]; data2 = data2[!is.na(data)]; data2 = data2[seq(1, length(data), by=7)]; data2 = embedd.custom(data2, c(1, 2, 365%/%7));
#data2 = read.csv("../../NOAA-vancouver-clean.csv"); data2 = data2[,"PRCP7"]; data2 = data2[!is.na(data)]; data2 = embedd.custom(data2, c(1, 2, 365));
#data = read.csv("../../NOAA-vancouver-clean.csv"); data = data[,"PRCP7"]; data = data[!is.na(data)]; data = embedd.custom(data, c(1, 2, 365));
#result = simulate(data, name="Whatever", divisions=10, plot.average=T);
#result = simulate(data2, name="NOAA Vancouver Precipitation", divisions=10, plot.average=T, newPlot=F, plotIdx=2);

#result = lorenz.single(data);

legend("topright", c("Lorenz", "Logistic", "NOAA Temperatures"), col=1:3, lty=1:3, lwd=2);

#savePlot("histogram-rule-lorenz.png");
