[<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/banner.png" width="888" alt="Visit QuantNet">](http://quantlet.de/)

## [<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/qloqo.png" alt="Visit QuantNet">](http://quantlet.de/) **pricing_kernel_uniform_confidence_bands** [<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/QN2.png" width="60" alt="Visit QuantNet 2.0">](http://quantlet.de/)

```yaml

Name of Quantlet: 'pricing_kernel_uniform_confidence_bands'
Published in: 'Journal of Financial Econometrics'
Description: 'Pricing kernels implicit in option prices play a key role in assessing the risk aversion over equity returns. We deal with nonparametric estimation of the pricing kernel (PK) given by the ratio of the risk-neutral density estimator and the historical density (HD). The former density can be represented as the second derivative w.r.t. the European call option price function, which we estimate by nonparametric regression. HD is estimated nonparametrically too. In this framework, we develop the asymptotic distribution theory of the Empirical Pricing Kernel (EPK) in the L∞ sense. Particularly, to evaluate the overall variation of the pricing kernel, we develop a uniform confidence band of the EPK. Furthermore, as an alternative to the asymptotic approach, we propose a bootstrap confidence band. The developed theory is helpful for testing parametric specifications of pricing kernels and has a direct extension to estimating risk aversion patterns. The established results are assessed and compared in a Monte-Carlo study. As a real application, we test risk aversion over time induced by the EPK.'
Keywords: Pricing Kernel, Confidence Band, Uniform Confidence Band, State Price Density, Rookley, DAX, Bootstrap
Author: Wolfgang Karl Härdle, Yarema Okhrin, Weining Wang
Submitted: 20150301

```

### R Code
```r

setwd("d:/papers/epk/program/")

source("EPK_library.r")

ngrid = 200;
bandwidth = 0.08;
tau = 0.04722;

read.dates = FALSE;

#all = read.table("tick2006.txt")
#day1 = all[all[,1]=="28-02-2006",];
#day2 = all[all[,1]=="04-04-2006",];
#write.table(day1, file="28-02-2006.txt", row.names=F, col.names=F, quote=F)
#write.table(day2, file="04-04-2006.txt", row.names=F, col.names=F, quote=F)

XX  = read.table("tick2006.txt");
#XX  = read.table("C20010117.dat");
#XX = XX[c(2,1,3,5,6,8,7,4)];
#XX = cbind(rep("17-01-2001", length=dim(XX)[1]), XX)

cat("reading finished ....\n");

if (read.dates==T)
{
    date = as.vector(XX[,1]);
    date.dif = date[1];
    date.part=date;
    flag=1;

    while(flag==1)
    {
        date.part = date.part[date.part!=date.dif[length(date.dif)]];
        cat("current date ... ",    date.dif[length(date.dif)], "\n");
        if (length(date.part>0)) {date.dif=c(date.dif, date.part[1]);} else {flag=0;}    
    }
    write.table(date.dif, file="tick2006_dates.txt", row.names=F, col.names=F, quote=F);
    } else { date.dif= read.table("tick2006_dates.txt"); date.dif=date.dif[,1];}

#############################################

#for( iday in c(1:length(date.dif)))
{
iday=42;

day1 = XX[XX[,1]==date.dif[iday],]; 
day1 = day1[day1[,2]>0,];
tau = day1[1,4];
#name = "28-02-2006";
#day1=read.table(paste(name,".txt", sep=""));

day1.mat = day1[day1[,4]== tau,];
day1.call = day1.mat[day1.mat[,3]==1,];
day1.put = day1.mat[day1.mat[,3]==0,];

# compute the moneyness
day1.call[,9]=day1.call[,7]/day1.call[,5];
day1.put[,9]=day1.put[,7]/day1.put[,5];

# only out and at the money options
day1.call = day1.call[day1.call[,9]>=1,];
day1.put = day1.put[day1.put[,9]<=1,];

# put to calls
put2call = day1.put;
put2call[,6] = put2call[,6]+ put2call[,7] -  put2call[,5]*exp(- mean(put2call[,8])* tau);

put2call = put2call[order(put2call[,5]),];
day1.call = day1.call[order(day1.call[,5]),];
data = rbind(put2call,day1.call);

data = data[(data[,2]>0.05),];

# no - arbitrage condition
data = data[ ((data[,6]<data[,7]) && (data[,6]>data[,7]-data[,5]*exp(-data[,8]*data[,4])) ),]

#write.table(data, file=paste(name,"_cleaned.txt", sep=""), row.names=F, col.names=F, quote=F)

n = dim(data)[1];

price.median = median(data[,7]);

## regression for volatility 

volas = data[,c(2,9)];
nsample = dim(volas)[1]; 

mon.min = min(volas[,2]);mon.max = max(volas[,2]);
volas[,2] = (volas[,2]-mon.min)/(mon.max-mon.min);
mon.grid = seq(1/ngrid, 1, length=ngrid);
mon.grid.scaled =  (mon.min+mon.grid*(mon.max-mon.min));

for(i in c(1:ngrid))
{
    X = cbind(rep(1, length=nsample), volas[,2]-mon.grid[i], (volas[,2]-mon.grid[i])^2, (volas[,2]-mon.grid[i])^3);
    W = diag(K.h(volas[,2]-mon.grid[i], bandwidth));
    beta = t(solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% volas[,1]);
    if (i==1) {sigmas=beta[c(1:3)];} else {sigmas=rbind(sigmas, beta[c(1:3)]);}
}

sigmas[,3] = 2*sigmas[,3];

#### bands for the derivatives

bands.init = bands(2, bandwidth, sigmas, volas[,2], volas[,1], mon.grid);

#### applying rookley

fder = rookley(mon.grid.scaled, sigmas[,1], sigmas[,2]/(mon.max-mon.min), sigmas[,3]/(mon.max-mon.min)^2, mean(data[,8]), tau);

#### final computations

strike.grid  =  price.median / mon.grid.scaled;
d2fdX2 = (mon.grid.scaled^2*fder[,3]+2*mon.grid.scaled*fder[,2])/strike.grid^2;

#### bands for SPD

d1=(log(mon.grid.scaled)+tau*(mean(data[,8])+0.5*(sigmas[,1])^2)) / (sqrt(sigmas[,1])*sqrt(tau));
d2=d1-sqrt(sigmas[,1])*sqrt(tau);
dgds= exp(-mean(data[,8])*tau) * ( mon.grid.scaled^2 * dnorm(d1) / strike.grid^2 - exp(-mean(data[,8])*tau)*dnorm(d2)/mon.grid.scaled);

band.limit = bands.init$L[3] * sqrt(bands.init$V[,3]) * dgds ; #* exp(-mean(data[,8])*tau); 

SPD = new.env();
SPD$SPD = price.median^2 * exp(mean(data[,8])*tau) * d2fdX2 / mon.grid.scaled;
SPD$lo = SPD$SPD - abs(price.median^2 * exp(mean(data[,8])*tau) * d2fdX2 * band.limit / mon.grid.scaled);
SPD$up = SPD$SPD + abs(price.median^2 * exp(mean(data[,8])*tau) * d2fdX2 * band.limit / mon.grid.scaled);
SPD=as.list(SPD);

dax=read.table("dax_index.dat"); dax=dax[,2];
dax=dax[c((length(dax)-500):length(dax))];

dax.scaled = 1+log(dax[c(ceiling(tau*31/0.0833):length(dax))]/dax[c(1:(length(dax)+1-ceiling(tau*31/0.0833)))]);
kern=density(dax.scaled, bw="nrd0", from=mon.grid.scaled[1], to=mon.grid.scaled[ngrid], n=ngrid);

########### Black-Scholes

sigma.BS = mean((volas[volas[,1]<quantile(volas[,1], probs=0.75),])[,1]);
q.BS = exp(-(log(mon.grid.scaled) - (mean(data[,8])-sigma.BS/2)*tau)^2/(2*tau*sigma.BS))/(mon.grid.scaled*sqrt(2*3.1415926*sigma.BS*tau));
p.BS = dnorm(log(mon.grid.scaled), mean=mean(log(dax.scaled)), sd = sd(log(dax.scaled)))/mon.grid.scaled;

########### plotting

pdf(paste("q_",date.dif[iday],".pdf", sep=""), width=10, height=6, onefile=F);
matplot(mon.grid.scaled, cbind(SPD$SPD, SPD$lo, SPD$up, q.BS), type="l",lty=1, lwd=c(2,1,1,2),col=c(2, "blue3", "blue3", "green3"), xlab="moneyness", ylab="q(K)", cex=1.5, main=paste(date.dif[iday], " tau = ", tau, "\n"), xlim=c(.95, 1.15), ylim=c(0,20));
dev.off()

pdf(paste("EPK_",date.dif[iday],".pdf", sep=""), width=10, height=6, onefile=F);
matplot(mon.grid.scaled, cbind(cbind(SPD$SPD, SPD$lo, SPD$up)/kern$y, q.BS/p.BS), type="l",lty=1, lwd=c(2,1,1,2),col=c(2, "blue3", "blue3", "green3"), xlab="moneyness", ylab="EPK", cex=1.5, main=paste(date.dif[iday], " tau = ", tau, "\n"), xlim=c(.95, 1.15), ylim=c(0,10));
dev.off()

pdf(paste("p_",date.dif[iday],".pdf", sep=""), width=10, height=6, onefile=F);
matplot(mon.grid.scaled, cbind(kern$y, p.BS), type="l",lty=1, lwd=c(2,2),col=c(2, "green3"), xlab="moneyness", ylab="q(K)", cex=1.5, main=paste(date.dif[iday], " tau = ", tau, "\n"), xlim=c(.95, 1.15), ylim=c(0,10));
dev.off()

mod = lm(log(SPD$SPD)~ log(mon.grid.scaled));
util.coef = mod$coef;

if (exists("util.param")==0) {util.param = c(exp(util.coef[1]), -util.coef[2], cor(log(SPD$SPD), log(mon.grid.scaled))^2);
    } else {util.coef = cbind(util.coef,c(exp(util.coef[1]), -util.coef[2], cor(log(SPD$SPD), log(mon.grid.scaled))^2));}
    
} # end of the days

```

automatically created on 2021-05-18