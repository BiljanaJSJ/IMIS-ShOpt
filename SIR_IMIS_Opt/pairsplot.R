

##############################################################
#generate plot for IMIS_Opt on SIR-ODE
##############################################################

out_ls=get(load('IMIS_opt_WidePrior1000.RData'))

#png('Plot_IMIS_WidePrior1000.png',width=500,height=400)
setEPS()
postscript("FIG5.eps",horizontal=FALSE, paper="special",height=18,width=24, colormodel = "cmyk", 
           family = "Helvetica")
par(mfrow=c(2,2),mar=c(10,10,10,10))

h=hist(out_ls$theta[,2],nclass=25,plot=F)
h$counts=h$counts/sum(h$counts)
plot(h,mgp=c(6,2.5,0),cex.main=4,cex.axis=4,cex.lab=4,xlim=c(0.08,0.12),xlab=expression(alpha),col='cyan',ylab='Density',main=bquote(paste("Posterior of " ~ alpha)))
box(which = "plot", lty = "solid")

h1=hist(out_ls$theta[,1],nclass=15,plot=F)
h1$counts=h1$counts/sum(h1$counts)
plot(h1,cex.main=4,cex.axis=4,cex.lab=4,xlim=c(0.0005,0.00075),mgp=c(6,2.5,0),xlab=expression(beta),col='cyan',ylab='Density',main=bquote(paste("Posterior of " ~ beta)))
box(which = "plot", lty = "solid")
plot(1, type="n", cex.main=4,cex.axis=4,cex.lab=4,xlim=c(3,7),mgp=c(6,2.5,0),xlab='I(0)',ylab='Density',main=bquote(paste("Posterior of " ~ I(0))))
box(which = "plot", lty = "solid")
abline(v=unique(out_ls$theta[,4]),col='cyan',lwd=3)
plot(out_ls$theta[,1],out_ls$theta[,2],cex.main=4,cex.axis=4,cex.lab=4,mgp=c(6,2.5,0),xlim=c(0.0005,0.00075),ylim=c(0.08,0.12),xlab=expression(beta),ylab=expression(alpha),main=bquote(paste("Joint posterior," ~ alpha ," and " ~ beta)))
dev.off()


