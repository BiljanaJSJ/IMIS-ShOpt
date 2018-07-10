####################################################################
#run IMIS-Opt for the SIR model
####################################################################


#source the file with functions
source('SIR_IMIS_Opt_functions.R')
data              = read.table("Eyam_time_SIR.csv",sep = ",")
data[data=="NaN"] = NA
data              = as.matrix(data)
colnames(data)    = c("time","S","I","R")
times             = data[,"time"]; 
dat               = data[,-1]
y                 = data[,3]
N=261

set.seed(345666)
library(parallel)
cl=makeCluster(4)
clusterExport(cl,varlist=ls(),envir = environment())
clusterCall(cl,function(x) {library(deSolve)})


t1=proc.time()[1]
out_ls=runIMIS_Opt(cl,N0=3*1000,D=3,B=1000,J=10000,niter=1000,
                       priorpars=c(1,1,1,1,N,5/N),N,dat,times,
                       par_true=c(0.0062,0.098,256,5))

(t=proc.time()[1]-t1)/60
save(out_ls,file='IMIS_opt_WidePrior1000.RData')  
stopCluster(cl)


