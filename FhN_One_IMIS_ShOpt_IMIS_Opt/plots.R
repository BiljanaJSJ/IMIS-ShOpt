#################################################################
#Generate plots
#################################################################
#setwd('D:/Publications/IMIS-ShOpt/incremental-mixture-importance-submitted/codeSubmit/codeSubmit/FhN_One_IMIS_ShOpt_IMIS_Opt')
rm(list = ls(all = TRUE))
setwd('E:/IMISCode/codeSubmit_keep_all_files/FhN_One_IMIS_ShOpt_IMIS_Opt')
source("Two-stage-FhN-just-c-with-prior.R")                         # 2-stage functions
source("IMIS.opt.colloc.proc-3optimizers.general-no-touchups.R")    # General IMIS 3 optimizers function 
source("fhn-model-set-up-x0proc-just-c.R")                          # likelihood etc...
source("FhN-model-set-up-as-ode-x0proc-thetalik-justc.R")		        # basic FhN functions
source("makeSSElik.R")
source("makeSSEprocFHN.R")

library(doParallel)
library(CollocInfer)


output_fullModel=get(load('F:/IMISCode/codeSubmit_keep_all_files/FhN_fullModel_IMIS_ShOpt/IMIS_shopt_full_fhn_D30.RData'))
output_1parModelIMIS_shOpt=get(load('FhN_1Param_IMIS_Shopt_D30.RData'))
output_IMIS_opt=get(load('FhN_1Param_IMIS_Opt_D3.RData'))
                          

#get solution and the data around c=mean(output_IMIS_opt$resample)
times = seq(0,20,0.2) 
print("Note that the parameter labelled 'sd' is actually a variance.  will fix this eventually")
x0	= c(-1,1) 
names(x0) = c("V","R")
pars=mean(output_IMIS_opt$resample)
parnames =names(pars)=c("c")
fhn=make.FHN()
y_c11	= lsoda(x0,times,fhn$fn.ode,pars) 
y_c11	= y_c11[,2:3] 
y=output_IMIS_opt$data
#data_c11 = y_c11 + matrix(rnorm(dim(y_c11)[1]*2,0,sqrt(.05^2)),length(times),2)

cl <- makeCluster(4)
registerDoParallel(cl)
clusterCall(cl,function(x) {library(deSolve);library(CollocInfer);library(numDeriv);library(lokern)})
clusterExport(cl,varlist=list('IMIS.opt.colloc.3optimizers.general.no.touch.ups','d2negnormdp2',"make.fhn","%dopar%","foreach",'make.SSEproc.FHN',"neglogprior","neglogpriorx0","prior","likelihood",'times',"dnegnormdp",'make.SSElik',"dneglogpriordpar","lokerns","ksLqudratic",'simex.fun.justc','neq','der.fhn.justc','jac.fhn.justc','d2neglogpriordpar2'))
clusterExport(cl,varlist=ls())

#output_IMIS_opt$data=data
#save(output_IMIS_opt,file='FhN_1Param_IMIS_Opt_D3.RData')
cgrid=seq(0.2,20,length=1000)
loglik=sapply(cgrid,function(x) {likelihood(x,logs=TRUE,data=output_IMIS_opt$data)})
logpost=sapply(cgrid,function(x) {prior(x,logs=TRUE)+likelihood(x,logs=TRUE,data=output_IMIS_opt$data)})
log_prior=sapply(cgrid,function(x) {prior(x,logs=TRUE)})

#stopCluster(cl)



setEPS()
postscript("FIG1.eps",horizontal=FALSE, paper="special",height=22,width=24, colormodel = "cmyk", 
    family = "Helvetica")
  
#png('loglikPriorIMIS.png',height = 450,width=600)
#par(mfrow=c(2,2),oma=c(3,2,rep(0,2))+2,mar=c(1,2,3,1))
#par(mfrow=c(2,2),oma=c(3,2,rep(1,2))+0.005,mar=c(8,8,8,8)+0.5)
par(mfrow=c(3,2),oma=c(3,2,rep(1,2))+0.005,mar=c(8,8,8,8)+0.5)
#unnormalized log posterior

plot(cgrid,(logpost),cex.lab=4,cex.axis=4,cex.main=4,mgp=c(6,2.5,0),xlab='c',ylab='density',main='A. Unnormalized log posterior')


#Likelihood over a coarse grid
plot(cgrid,loglik,cex.lab=4,cex.axis=4,cex.main=4,mgp=c(6,2.5,0),xlab='c',ylab='density',main='B.Log likelihood')


#log prior
plot(cgrid,log_prior,cex.lab=4,cex.axis=4,cex.main=4,mgp=c(6,2.5,0),xlab='c',ylab='density',main='C.Log prior c~N(14,2)')

#IMIS-Opt posterior estimate


plot(density(output_IMIS_opt$resample),xlab='c',ylab='density',mgp=c(6,2.5,0),xlim=range(cgrid),main='D.IMIS-Opt posterior density',cex.axis=4,cex.main=4,cex.lab=4)
##IMIS-Opt posterior estimate zoomed in
##par(new=TRUE, oma=c(5,6,0,0.005))
plot(density(output_1parModelIMIS_shOpt$resample),xlab='c',ylab='density',mgp=c(6,2.5,0),xlim=range(cgrid),main='E.IMIS-ShOpt posterior density',cex.axis=4,cex.main=4,cex.lab=4)

par(new=TRUE, oma=c(9,13,1,0))
# ##par(new=TRUE, oma=c(5,6,0,0.005))
#
matLayout=matrix(0,6, 4, byrow = TRUE)
#matLayout[2,1]=1; matLayout[2,2]=1
#matLayout[3,1]=1; matLayout[3,2]=1
#matLayout[3,3]=1
matLayout[4,3]=1; #matLayout[4,2]=1
layout(matLayout)
plot(density(output_IMIS_opt$resample,adj=6), col='red',main='IMIS-Opt:zoom in',cex.axis=3,cex.main=4,cex.lab=3,xlab='',lwd=2,mgp=c(1.5,1,0),ylab='')

#dev.off()

# setEPS()
# postscript("FIG11.eps",horizontal=FALSE, paper="special",height=9,width=12, colormodel = "cmyk", 
# 					 family = "Helvetica")
# par(mfrow=c(1,1))
# par(mar=c(8,8,8,8)+0.5)#c(7,9,5,5))
#plot(density(output_1parModelIMIS_shOpt$resample),xlab='c',ylab='density',mgp=c(6,2.5,0),xlim=range(cgrid),main='E.IMIS-ShOpt posterior density',cex.axis=4,cex.main=4,cex.lab=4)

 par(new=TRUE, oma=c(9,4,0,12))
# 
# # matLayout=matrix(0,3, 3, byrow = TRUE)
# # matLayout[2,2]=1; matLayout[2,3]=1
# # matLayout[3,2]=1; matLayout[3,3]=1
# 
matLayout=matrix(0,6, 4, byrow = TRUE)
matLayout[6,2]=1

layout(matLayout)
plot(density(output_1parModelIMIS_shOpt$resample,adj=6), col='red',main='IMIS-ShOpt:zoom in',cex.axis=3,cex.main=4,cex.lab=3,xlab='',lwd=2,mgp=c(1.5,1,0),ylab='')

dev.off()



png('StateVariablesData.png',width=700,height = 400)

par(mfrow=c(1,2),oma=c(3,2,rep(0,2))+0.05,mar=c(2,1,3,1)+2)


plot(times,y_c11[,1],col='blue',main=paste('A.State variables and obs., \n c=',round(mean(output_IMIS_opt$resample),2),sep=''),cex.axis=2,cex.main=2,cex.lab=2,xlab='times',ylab='',lwd=3,lty=1)
lines(times,y_c11[,2],col='green',lwd=2)
points(times,output_IMIS_opt$data[,1],col='red',lwd=3)
points(times,output_IMIS_opt$data[,2],col='orange',lwd=2)
par(xpd=TRUE)
legend(-0.00015,-0.95,c('V','R',expression(Y[V]),expression(Y[R])),cex=0.85,lty = c(1, 1, NA,NA), pch = c(NA, NA,1,1),col=c('blue',"green","red",'orange'),lwd=c(3,2,3,2))


plot(times,y[,1],col='blue',main='B. State variables and obs., \n c=3',cex.axis=2,cex.main=2,cex.lab=2,xlab='times',ylab='',type='l',lwd=3)
lines(times,y[,2],col='green',lwd=2)
points(times,output_IMIS_opt$data[,1],col='red',lwd=3)
points(times,output_IMIS_opt$data[,2],col='orange',lwd=2)
par(xpd=TRUE)
legend(-0.00015,-1,c('V','R',expression(Y[V]),expression(Y[R])),cex=0.85,lty = c(1, 1, NA,NA), pch = c(NA, NA,1,1),col=c('blue',"green","red",'orange'),lwd=c(3,2,3,2))

dev.off()


#caclulate KL divergence

library(flexmix)
#evaluate the density of IMIS-ShOpt samples over the interval [2.5,3.5]
#use marginal likelihood obtained from the IMIS-ShOpt as normalizing constant
dsamples=density(output_1parModelIMIS_shOpt$resample, from=2.5, to=3.5)
normalizedsamples=dsamples$y/exp(output_1parModelIMIS_shOpt$stat[2,1])

#evaluate the theoretical density over the same interval
cgrid=seq(2.5,3.5,length=length(dsamples$y))
logpost=sapply(cgrid,function(x) {prior(x,logs=TRUE)+likelihood(x,logs=TRUE,data=output_IMIS_opt$data)})

#numerically integrate the target posterior to obtain the normalizing constant
normconstintegrand <- function(x) {exp(prior(x,logs=TRUE)+likelihood(x,logs=TRUE,data=output_IMIS_opt$data))}
normconst=integrate(normconstintegrand,lower=2.5,upper=3.5)
norm_post=exp(logpost)/normconst$value


stopCluster(cl)


plot(cgrid,norm_post)
lines(dsamples$x,normalizedsamples,col='red')




plot(dsamples$x,normalizedsamples)
plot(cgrid,norm_post)

KLdiv(cbind(norm_post,normalizedsamples))





PlotsResampledTraj=function(times=seq(0,20,0.2),output,title,filename){
  
  if (is.vector(output$resample)){
    getTraj=lapply(1:length(output$resample), function(x) lsoda(x0,times,fhn$fn.ode, output$resample[x] ))
    meanSol=lsoda(x0,times,fhn$fn.ode, mean(output$resample) )
  }else{
    getTraj=lapply(1:nrow(output$resample), function(x) lsoda(x0,times,fhn$fn.ode, output$resample[x,] ))
    meanSol=lsoda(x0,times,fhn$fn.ode, colMeans(output$resample) )
  }
  #png(filename)
  
  setEPS()
  postscript(filename,horizontal=FALSE, paper="special",height=18,width=24, colormodel = "cmyk", 
             family = "Helvetica")
  
  
  par(mfrow=c(1,1),mar=c(9,9,9,9))
  plot(times,getTraj[[1]][,'V'], col='grey',main=title,type='l',mgp=c(7,2.5,0),ylab='',xlab='time',cex.axis=4,cex.main=4,cex.lab=4,lwd=2)
  for (i in (2:length(getTraj))){
    lines(times,getTraj[[i]][,'V'],col='grey',lwd=5)
  }  
  lines(times,meanSol[,'V'],col='blue',lwd=3)
  points(seq(0,20,0.2),output$data[,1],pch=21,col='red',lwd=1.5)
  
 
 
  lines(times,getTraj[[1]][,'R'], col='grey',main='',type='l',lwd=3)
  for (i in (2:length(getTraj))){
    lines(times,getTraj[[i]][,'R'],col='grey',lwd=5)
  }  
  lines(times,meanSol[,'R'],col='darkgreen',lwd=2)
  points(seq(0,20,0.2),output$data[,2],pch=21,col='red',lwd=1.5)
  legend(-0.2,-1.05,c("V","R",'data'),cex=4,lty=c(1,1,NA),pch=c(NA,NA,21),col=c("blue","darkgreen",'red'),lwd=c(3,2,1.5))
  dev.off()
}




par(mfrow=c(1,1))

PlotsResampledTraj(times=seq(0,20,0.2),output=output_IMIS_opt,title='A.IMIS-Opt, one parameter FhN',filename="FIG2.eps")

PlotsResampledTraj(times=seq(0,20,0.2),output=output_1parModelIMIS_shOpt,title='B.IMIS-ShOpt, one parameter FhN',filename="FIG3.eps")

PlotsResampledTraj(times=seq(0,20,0.2),output=output_fullModel,title='C.IMIS-ShOpt, full FhN',filename="FIG4.eps")


#PlotsResampledTraj(times=seq(0,20,0.2),output=output_1parModelIMIS_shOpt,title='IMIS-ShOpt samples, FhN-ODE model',filename='Oneparam_IMIS_shOpt_Splunk.png')

setEPS()
postscript("FIG9.eps",horizontal=FALSE, paper="special",height=18,width=24, colormodel = "cmyk", 
           family = "Helvetica")

par(mfrow=c(1,1),mar=rep(5,4))
h1=hist(output_1parModelIMIS_shOpt$X_all[1:10000],breaks=50,plot=F)
h1$counts=h1$counts/sum(h1$counts)
h2=hist(output_1parModelIMIS_shOpt$X_all[10001:40000],breaks=30,plot=F)
h2$counts=h2$counts/sum(h2$counts)
rangey=max(range(h1$counts)[2],range(h2$counts)[2])


plot(h1,main='Importance sampling distribution at the end of the Shotgun optimization stage',xlab='c',ylab='Density', xlim=c(0,20),ylim=c(0,rangey),col='coral3',cex.main=3,cex.lab=2.5,cex.axis=2.5,mgp=c(3.5,1.5,0))
points(output_1parModelIMIS_shOpt$center[1],0,col='red',pch=21,bg=2,cex=3)
par(new=T)
plot(h2,main='',xlab='',ylab='',  xlim=c(0,20),ylim=c(0,rangey),cex.main=2,cex.lab=2.5,cex.axis=2.5,col='darkcyan',mgp=c(3.5,1.5,0))
points(output_1parModelIMIS_shOpt$center[2],0,col='darkcyan',pch=16,bg=2,cex=3)
points(output_1parModelIMIS_shOpt$center[3],0,col='darkcyan',pch=16,bg=2,cex=3)
lines(density(output_1parModelIMIS_shOpt$resample)$x,density(output_1parModelIMIS_shOpt$resample)$y/245,col='red',lwd=10);
#hist(output_1parModelIMIS_shOpt$resample,col='red', add=T, nbins=10);
text(13,0.47,'NLS, GP',cex=3,col='darkcyan')
text(3,0.25,'Two-Stage',cex=3,col='darkcyan')
text(17,0.13,'Initial \nimportance sampling \ndistribution',cex=3,col='coral3')

dev.off()

# plot(density(output_1parModelIMIS_shOpt$center[1:30]),cex.axis=4,cex.main=4,cex.lab=4,main='D. Shotgun opimization: discovered modes', xlim=c(0,20), xlab='c',mgp=c(6,1.5,0))
# text(12,0.12,'NLS, GP',cex=2.5)
# text(3,0.07,'Two-Stage',cex=2.5)
