#################################################################
#Generate plots
#################################################################

source("Two-stage-FhN-just-c-with-prior.R")                         # 2-stage functions
source("IMIS.opt.colloc.proc-3optimizers.general-no-touchups.R")    # General IMIS 3 optimizers function 
source("fhn-model-set-up-x0proc-just-c.R")                          # likelihood etc...
source("FhN-model-set-up-as-ode-x0proc-thetalik-justc.R")		        # basic FhN functions
source("makeSSElik.R")
source("makeSSEprocFHN.R")

library(doParallel)
library(CollocInfer)


output_fullModel=get(load('D:/Publications/IMIS-ShOpt/code/FhN_fullModel_IMIS_ShOpt/IMIS_shopt_full_fhn_D30.RData'))
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

stopCluster(cl)


setEPS()
postscript("FIG1.eps",horizontal=FALSE, paper="special",height=18,width=24, colormodel = "cmyk", 
    family = "Helvetica")
  
#png('loglikPriorIMIS.png',height = 450,width=600)
#par(mfrow=c(2,2),oma=c(3,2,rep(0,2))+2,mar=c(1,2,3,1))
par(mfrow=c(2,2),oma=c(3,2,rep(1,2))+0.005,mar=c(8,8,8,8)+0.5)
#unnormalized log posterior
plot(cgrid,logpost,cex.lab=4,cex.axis=4,cex.main=4,mgp=c(6,2.5,0),xlab='c',ylab='density',main='A.Un-normalized log posterior')

#Likelihood over a coarse grid
plot(cgrid,loglik,cex.lab=4,cex.axis=4,cex.main=4,mgp=c(6,2.5,0),xlab='c',ylab='density',main='B.Log likelihood')


#log prior
plot(cgrid,log_prior,cex.lab=4,cex.axis=4,cex.main=4,mgp=c(6,2.5,0),xlab='c',ylab='density',main='C.Log prior c~N(14,2)')
#IMIS-Opt posterior estimate
plot(density(output_IMIS_opt$resample),xlab='c',ylab='density',mgp=c(6,2.5,0),xlim=range(cgrid),main='D.IMIS-Opt posterior estimate',cex.axis=4,cex.main=4,cex.lab=4)
##IMIS-Opt posterior estimate zoomed in
##par(new=TRUE, oma=c(5,6,0,0.005))
par(new=TRUE, oma=c(8,11,1,0))
##par(new=TRUE, oma=c(5,6,0,0.005))

matLayout=matrix(0,4, 4, byrow = TRUE)
matLayout[4,3]=1
layout(matLayout)
plot(density(output_IMIS_opt$resample,adj=6), col='red',main='IMIS-Opt post.\n est.--zoomed in',cex.axis=3,cex.main=4,cex.lab=3,xlab='',lwd=2,mgp=c(1.5,1,0),ylab='')

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



# par(mfrow=c(1,1))
# plot(cgrid,log_prior,cex.axis=1.5,main='c). log prior c~N(14,2)')
# plot(density(output$resample-range(log_prior)[1]),ylim=range(log_prior),col='red')


#make plots of trajectories in all the points from the posterior distribution
#for full model, 1 param model IMIS-Shopt and 1 param model IMIS-opt
# times=seq(0,20,length.out=100)
# getTraj=lapply(1:length(output$resample), function(x) lsoda(x0,times,fhn$fn.ode, output$resample[x] ))
# meanSol=lsoda(x0,times,fhn$fn.ode, mean(output$resample) )
# 
# 
# par(mar=rep(2,4),mfrow=c(1,1))
# plot(times,getTraj[[1]][,'V'], col='grey',main='Resampled trajectories, model 1',type='l')
# for (i in (2:length(getTraj))){
#   lines(times,getTraj[[i]][,'V'],col='grey')
# }  
# lines(times,meanSol[,'V'],col='blue',lwd=1.5)
# points(seq(0,20,0.2),output$data[,1],pch=21,col='red')
# 
# 
# lines(times,getTraj[[1]][,'R'], col='grey',main='',type='l')
# for (i in (2:length(getTraj))){
#   lines(times,getTraj[[i]][,'R'],col='grey')
# }  
# lines(times,meanSol[,'R'],col='darkgreen',lwd=1.5)
# points(seq(0,20,0.2),output$data[,2],col='red',type = "p")
# legend(-0.77,-1.35,c("V","R",'data'),cex=0.75,lty=c(1,1),col=c("blue","darkgreen",'red'))
# 
# 








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

PlotsResampledTraj(times=seq(0,20,0.2),output=output_IMIS_opt,title='A.IMIS-Opt, Model 1',filename="FIG2.eps")

PlotsResampledTraj(times=seq(0,20,0.2),output=output_1parModelIMIS_shOpt,title='B.IMIS-ShOpt, Model 1',filename="FIG3.eps")

PlotsResampledTraj(times=seq(0,20,0.2),output=output_fullModel,title='C.IMIS-ShOpt, Model 2',filename="FIG4.eps")


#PlotsResampledTraj(times=seq(0,20,0.2),output=output_1parModelIMIS_shOpt,title='IMIS-ShOpt samples, FhN-ODE model',filename='Oneparam_IMIS_shOpt_Splunk.png')


