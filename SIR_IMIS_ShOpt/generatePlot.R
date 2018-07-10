####################################################################
#generate plots for the SIR model using samples from the IMIS-ShOpt
####################################################################

source('pairs2_function.R')
out_ls=get(load('samples.RData'))


panel.hist <- function(x, ...)
{ 
  usr <- par("usr"); on.exit(par(usr))
  
  par(usr = c(usr[1:2], 0,1.5) )
  
  if (unique(x) %in% c(1,2,3,4,5,6,7)){
  
    h <- hist(x, plot = FALSE,breaks=seq(range(x)[1],range(x)[2],by=0.2))  
  }else{
    
    h <- hist(x, plot = FALSE,breaks=150)
  }
  
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...)
}


#pdf('IMIS_shopOut_SIR_new_widePrior100000Iter.pdf',width=15,height=12,pointsize=15)
# png('IMIS_shopOut_SIR_new_widePrior100000Iter.png',width=1200,height = 700)
# pairs(out_ls$theta[,c(2,1,4)],labels=c(expression(alpha),expression(beta),expression(I(0)),expression(tau)), bg="light blue",
#       diag.panel=panel.hist,cex.labels = 2.5, font.labels=2,cex.axis=2,upper.panel=NULL,gap = 3)
# dev.off()


setEPS()
postscript("FIG6.eps",horizontal=FALSE, paper="special",height=15,width=19, colormodel = "cmyk", 
           family = "Helvetica")
pairs2(out_ls$theta[,c(2,1,4)],labels=c(expression(alpha),expression(beta),expression(I(0)),expression(tau)), bg="cyan",
       diag.panel=panel.hist,cex.labels = 4, font.labels=4,cex.axis=4,upper.panel=NULL, oma=c(2,3,2,2),gap=8,label.pos =0.9,mgp = c(5, 3, 1))
dev.off()

