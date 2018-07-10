#########################################################################
#run the IMIS-ShOpt on the SIR model
#########################################################################

#source the file with functions to run the IMIS-ShOpt on the SIR data
#and read the data
rm(list=ls())
source('SIR_IMIS_ShOpt_functions.R')
data              = read.table("Eyam_time_SIR.csv",sep = ",")
data[data=="NaN"] = NA
data              = as.matrix(data)
colnames(data)    = c("time","S","I","R")
times             = data[,"time"]; 
dat               = data[,-1]
y                 = data[,3]
N=261


#call the library parallel, to run the particles in parallel
#and setup the cluster
library(parallel)
set.seed(345666)
cl=makeCluster(4)
clusterExport(cl,varlist=ls(),envir = environment())
clusterCall(cl,function(x) {library(deSolve)})

#run the algorithm
t1=proc.time()[1]
out_ls=runIMIS_ShOpt(cl,N0=3*1000,D=3,R=10,B=1000,J=10000,niter=1000,
                       priorpars=c(1,1,1,1,N,5/N),N,dat,times,
                       par_true=c(0.0062,0.098,256,5))
(t=proc.time()[1]-t1)/60
stopCluster(cl)



#save the results
save(out_ls,file='IMIS_Shopt_new_WidePrior100000iter.RData') 


