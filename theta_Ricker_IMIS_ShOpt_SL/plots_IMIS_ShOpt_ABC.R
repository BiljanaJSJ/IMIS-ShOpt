###############################################################################
#generate plots
###############################################################################

#load the data
de1=get(load('LoadDataEps_theta_RickerData_newData_poiss_neweps_0.3_trynew.RData'))
obs.data=de1$obs.data

theta.true=obs.data$theta
eps=de1$eps
N0.true=de1$N0.true
timestep=length(de1$obs.data$x)


#load the samples obtained from the IMIS-ShOpt-ABC 
samples=get(load('samples.RData'))
X_all=samples$X_all

X_all=X_all[-((dim(X_all)[1]-999):dim(X_all)[1]),]
Weights=samples$Weights



#marginal posteriors for the model
#png('histogram_1000_3000_M2_new_r_0_5_1_sig2_10_6_thet_1_5_1.png',width=650,height = 500)
setEPS()
postscript("FIG7.eps",horizontal=FALSE, paper="special",height=18,width=19, colormodel = "cmyk", 
           family = "Helvetica")

m <- matrix(c(1,2,3,4,5,5),nrow = 3,ncol = 2,byrow = TRUE)
layout(mat = m,heights = c(0.8,0.8,0.2))
#par(mar = c(6.5,5,3,2.5))
par(mar = c(10,10,10,10))
plot(density(samples$resample[,1],adjust=1.5),lwd=2,xlim=range(samples$X_all[,1]),main=bquote(paste("Posterior of log" ~ r)),cex.axis=4,cex=4,xlab='log r',cex.main=4,cex.lab=4,mgp=c(6,2.5,0))
xx <- seq(range(samples$X_all[,1])[1],range(samples$X_all[,1])[2],length=200)
lines(xx, dnorm(xx,theta.true[1],1),lwd=2.5,col="darkgrey")
abline(v=mean(samples$resample[,1]),col='blue',lwd=2,lty=2)
abline(v=theta.true[1],col='red',lwd=2,lty=3)




densphi=density(samples$resample[,2],adjust=1.5)
plot(densphi,lwd=2,xlim=c(0,15),main=bquote(paste("Posterior of " ~ phi)),cex.axis=4,cex=4,xlab=expression(phi),cex.main=4,cex.lab=4,mgp=c(6,2.5,0))
xx <- seq(0,15,length=200)
lines(xx, dchisq(xx,df=4),lwd=2.5,col="darkgrey")
abline(v=mean(samples$resample[,2]),col='blue',lwd=2,lty=2)
abline(v=theta.true[2],col='red',lwd=2,lty=3)


densSig2=density(samples$resample[,3],adjust=1)
plot(densSig2,mgp=c(7,2,0), xlim=c(0,0.4), lwd=2,main=bquote(paste("Posterior of " ~ sigma[p]^2)),cex.axis=4,cex=4,xlab=expression(sigma[p]^2),cex.main=4,cex.lab=4)
xx <- seq(0,0.4,length=200)
lines(xx,10^2.5*dgamma(1/xx,shape=2, scale=1/0.05),lwd=2.5,col="darkgrey")
abline(v=mean(samples$resample[,3]),col='blue',lwd=2,lty=2)
abline(v=theta.true[3],col='red',lwd=2,lty=3)


plot(density(samples$resample[,4],adjust = 2),lwd=2,xlim=c(-3,5),main=bquote(paste("Posterior of log " ~ tilde(theta))),mgp=c(6,2.5,0),cex.axis=4,cex=4,xlab=expression(paste('log ',tilde(theta))),cex.main=4,cex.lab=4)
xx <- seq(-3,5,length=200)
lines(xx, dnorm(xx,theta.true[4],1),lwd=2.5,col="darkgrey")
abline(v=mean(samples$resample[,4]),col='blue',lwd=2,lty=2)
abline(v=theta.true[4],col='red',lwd=2,lty=3)

par(mar = c(0.8,0.8,0.8,0.8))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend(x = "bottom",inset = 0,
       legend = c('Posterior','Prior','Posterior mean','True value'), 
       col=c('black','darkgrey','blue','red'), 
       lty=c(1,1,2,3),lwd=2, cex=4, horiz = TRUE)


dev.off()









#weigths for model M2 over the full importance sampling distribution before the re-sampling
#png('Weights_M2.png',width=600,height = 500)
setEPS()
postscript("FIG8.eps",horizontal=FALSE, paper="special",height=18,width=19, colormodel = "cmyk", 
           family = "Helvetica")

m <- matrix(c(1,2,3,4,5,5),nrow = 3,ncol = 2,byrow = TRUE)
layout(mat = m,heights = c(0.8,0.8,0.2))
par(mar = rep(10,4))


plot(X_all[,1],Weights,xlab='log r',mgp=c(6,2.5,0),ylab='Weights',main=bquote(paste("Weights log" ~ r)),cex.axis=4,cex=2,cex.main=4,cex.lab=4)
abline(v=theta.true[1],col='red',lwd=2)
plot(X_all[,2],Weights,ylab='Weights',mgp=c(6,2.5,0),main=bquote(paste("Weights " ~ phi)),cex.axis=4,cex=4,xlab=expression(phi),cex.main=4,cex.lab=4)
abline(v=theta.true[2],col='red',lwd=2)
plot(X_all[,3],Weights,ylab='Weights',mgp=c(6,2.5,0),main=bquote(paste("Weights " ~ sigma[p]^2)),cex.axis=4,cex=4,xlab=expression(sigma[p]^2),cex.main=4,cex.lab=4)
abline(v=theta.true[3],col='red',lwd=2)
plot(X_all[,4],Weights,ylab='Weights',mgp=c(6,2.5,0),main=bquote(paste("Weights log" ~ tilde(theta))),cex.axis=4,cex=4,xlab=expression(paste('log ',tilde(theta))),cex.main=4,cex.lab=4)
abline(v=theta.true[4],col='red',lwd=2)


par(mar = c(0.8,0.8,0.8,0.8))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend(x = "bottom",inset = 0,
       legend = c('Weights','True parameter value'), 
       col=c('black','red'), lwd=5, cex=4, horiz = TRUE)

dev.off()




