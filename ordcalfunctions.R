##########################################################
# ORDINAL MODEL CALIBRATION AND DISCRIMINATION FUNCTIONS #
##########################################################

# CALIBRATION PER OUTCOME SEPARATELY
calout <- function(out,preds,k){
  cores = matrix(0,k,2)
  for (i in (1:k)){
    cores[i,1] = glm(out==i ~ 1, offset=logit(preds[,i]), family=binomial(link='logit'))$coefficients
    cores[i,2] = glm(out==i ~ logit(preds[,i]), family=binomial(link='logit'))$coefficients[2]
    
  }
  cores
}

# output = calout(out=cad$o5,preds=mlrpred,5)


# STANDARD COX RECALIBRATION WITH DICHOTOMIZED OUTCOMES
caldout <- function(out,preds,k){
  cores = matrix(0,k-1,2)
  for (i in (2:(k-1))){
    cores[i-1,1] = glm(out>=i ~ 1, offset=logit(rowSums(preds[,i:k])), family=binomial(link='logit'))$coefficients
    cores[i-1,2] = glm(out>=i ~ logit(rowSums(preds[,i:k])), family=binomial(link='logit'))$coefficients[2]
  }
  cores[k-1,1] = glm(out>=k ~ 1, offset=logit(preds[,k]), family=binomial(link='logit'))$coefficients
  cores[k-1,2] = glm(out>=k ~ logit(preds[,k]), family=binomial(link='logit'))$coefficients[2]
  cores
}

# output = caldout(out=cad$o5,preds=mlrpred,k=5)


# GET E
getE <- function(out,preds,k){
  Ec=preds
  for (i in (1:k)){Ec[,i]=i*Ec[,i]}
  rowSums(Ec)
}

# OVERALL RECALIBRATION BASED ON E
calE <- function(out,preds,k){
  Ec=preds
  for (i in (1:k)){Ec[,i]=i*Ec[,i]}
  E=rowSums(Ec)
  Ei = glm(as.numeric(out) ~ 1, offset=log(E), family=poisson(link="log"), na.action = "na.omit")$coefficients
  Es = glm(as.numeric(out) ~ log(E), family=poisson(link="log"))$coefficients[2]
  c(Ei,Es)
}

# output = calE(out=cad$o5,preds=mlrpred,k=5)


# OVERALL RECALIBRATION BASED ON E USING GLM
calglmE <- function(out,preds,k){
  Ec=preds
  for (i in (1:k)){Ec[,i]=i*Ec[,i]}
  E=rowSums(Ec)
  Ei = glm(as.numeric(out) ~ 1, offset=E, na.action = "na.omit")$coefficients
  Es = glm(as.numeric(out) ~ E)$coefficients[2]
  c(Ei,Es)
}


# CALIBRATION PLOT WITH POINTS
plotscatter3 <- function(preds,obs,cexv=0.8,cexaxisv=1.4,cexlabv=1.5,cexlegendv=1.2){
  plot(preds[,1],fitted(obs)[,1],type="p",pch=1,col="green",lwd=1,ylab="Observed proportion",xlab="Estimated probability",xlim=0:1,ylim=0:1,cex=cexv,cex.axis=cexaxisv,cex.lab=cexlabv)
  points(preds[,2],fitted(obs)[,2],pch=1,col="orange",lwd=1,cex=cexv)
  points(preds[,3],fitted(obs)[,3],pch=1,col="red",lwd=1,cex=cexv)
  lines(c(0,1),c(0,1),type="l",col="gray",lty=3,lwd=2) # plot the ideal diagonal line
  legend(0, 1, col=c("gray","green","orange","red"), c("Ideal", "Cat. 1","Cat. 2", "Cat. 3"), lty=c(3,NA,NA,NA), lwd=c(2,NA,NA,NA), pch=c(NA,1,1,1), box.col="black", cex=cexlegendv, seg.len=2)
}
plotscatter4 <- function(preds,obs,cexv=0.8,cexaxisv=1.4,cexlabv=1.5,cexlegendv=1.2){
  plot(preds[,1],fitted(obs)[,1],type="p",pch=1,col="green",lwd=1,ylab="Observed proportion",xlab="Estimated probability",xlim=0:1,ylim=0:1,cex=cexv,cex.axis=cexaxisv,cex.lab=cexlabv)
  points(preds[,2],fitted(obs)[,2],pch=1,col="orange",lwd=1,cex=cexv)
  points(preds[,3],fitted(obs)[,3],pch=1,col="red",lwd=1,cex=cexv)
  points(preds[,4],fitted(obs)[,4],pch=1,col="brown",lwd=1,cex=cexv)
  lines(c(0,1),c(0,1),type="l",col="gray",lty=3,lwd=2) # plot the ideal diagonal line
  legend(0, 1, col=c("gray","green","orange","red","brown"), c("Ideal", "Cat. 1","Cat. 2", "Cat. 3", "Cat. 4"), lty=c(3,NA,NA,NA,NA), lwd=c(2,NA,NA,NA,NA), pch=c(NA,1,1,1,1), box.col="black", cex=cexlegendv, seg.len=2)
}
plotscatter5 <- function(preds,obs,cexv=0.8,cexaxisv=1.4,cexlabv=1.5,cexlegendv=1.2){
  plot(preds[,1],fitted(obs)[,1],type="p",pch=1,col="green",lwd=1,ylab="Observed proportion",xlab="Estimated probability",xlim=0:1,ylim=0:1,cex=cexv,cex.axis=cexaxisv,cex.lab=cexlabv)
  points(preds[,2],fitted(obs)[,2],pch=1,col="orange",lwd=1,cex=cexv)
  points(preds[,3],fitted(obs)[,3],pch=1,col="red",lwd=1,cex=cexv)
  points(preds[,4],fitted(obs)[,4],pch=1,col="brown",lwd=1,cex=cexv)
  points(preds[,5],fitted(obs)[,5],pch=1,col="black",lwd=1,cex=cexv)
  lines(c(0,1),c(0,1),type="l",col="gray",lty=3,lwd=2) # plot the ideal diagonal line
  legend(0, 1, col=c("gray","green","orange","red","brown","black"), c("Ideal", "Cat. 1","Cat. 2", "Cat. 3", "Cat. 4", "Cat. 5"), lty=c(3,NA,NA,NA,NA,NA), lwd=c(2,NA,NA,NA,NA,NA), pch=c(NA,1,1,1,1,1), box.col="black", cex=cexlegendv, seg.len=2)
}

#plotscatter5(preds=mlrpred,obs=mlrvgamsmps4,mtxt="MLR")


# CALIBRATION PLOT WITH SMOOTHED LINES
plotlines3 <- function(preds,obs,cexaxisv=1.4,cexlabv=1.5,cexlegendv=1.2){
  wa1=smooth.spline(preds[,1],fitted(obs)[,1])
  wa2=smooth.spline(preds[,2],fitted(obs)[,2])
  wa3=smooth.spline(preds[,3],fitted(obs)[,3])
  plot(wa1$x, wa1$y,type="l",col="green",lwd=3,ylab="Observed proportion",xlab="Estimated probability",xlim=0:1,ylim=0:1,cex.axis=cexaxisv,cex.lab=cexlabv)
  lines(wa2$x, wa2$y,col="orange",lwd=3)
  lines(wa3$x, wa3$y,col="red",lwd=3)
  lines(c(0,1),c(0,1),type="l",col="gray",lty=3,lwd=2) # plot the ideal diagonal line
  legend(0, 1, col=c("gray","green","orange","red"), c("Ideal", "Cat. 1","Cat. 2", "Cat. 3"), lty=c(3,1,1,1,1), lwd=c(2,3,3,3,3), box.col="black", cex=cexlegendv, seg.len=2)
}
plotlines4 <- function(preds,obs,cexaxisv=1.4,cexlabv=1.5,cexlegendv=1.2){
  wa1=smooth.spline(preds[,1],fitted(obs)[,1])
  wa2=smooth.spline(preds[,2],fitted(obs)[,2])
  wa3=smooth.spline(preds[,3],fitted(obs)[,3])
  wa4=smooth.spline(preds[,4],fitted(obs)[,4])
  plot(wa1$x, wa1$y,type="l",col="green",lwd=3,ylab="Observed proportion",xlab="Estimated probability",xlim=0:1,ylim=0:1,cex.axis=cexaxisv,cex.lab=cexlabv)
  lines(wa2$x, wa2$y,col="orange",lwd=3)
  lines(wa3$x, wa3$y,col="red",lwd=3)
  lines(wa4$x, wa4$y,col="brown",lwd=3)
  lines(c(0,1),c(0,1),type="l",col="gray",lty=3,lwd=2) # plot the ideal diagonal line
  legend(0, 1, col=c("gray","green","orange","red","brown"), c("Ideal", "Cat. 1","Cat. 2", "Cat. 3", "Cat. 4"), lty=c(3,1,1,1,1), lwd=c(2,3,3,3,3), box.col="black", cex=cexlegendv, seg.len=2)
}
plotlines5 <- function(preds,obs,cexaxisv=1.4,cexlabv=1.5,cexlegendv=1.2){
  wa1=smooth.spline(preds[,1],fitted(obs)[,1])
  wa2=smooth.spline(preds[,2],fitted(obs)[,2])
  wa3=smooth.spline(preds[,3],fitted(obs)[,3])
  wa4=smooth.spline(preds[,4],fitted(obs)[,4])
  wa5=smooth.spline(preds[,5],fitted(obs)[,5])
  plot(wa1$x, wa1$y,type="l",col="green",lwd=3,ylab="Observed proportion",xlab="Estimated probability",xlim=0:1,ylim=0:1,cex.axis=cexaxisv,cex.lab=cexlabv)
  lines(wa2$x, wa2$y,col="orange",lwd=3)
  lines(wa3$x, wa3$y,col="red",lwd=3)
  lines(wa4$x, wa4$y,col="brown",lwd=3)
  lines(wa5$x, wa5$y,col="black",lwd=3)
  lines(c(0,1),c(0,1),type="l",col="gray",lty=3,lwd=2) # plot the ideal diagonal line
  legend(0, 1, col=c("gray","green","orange","red","brown","black"), c("Ideal", "Cat. 1","Cat. 2", "Cat. 3", "Cat. 4", "Cat. 5"), lty=c(3,1,1,1,1,1), lwd=c(2,3,3,3,3,3), box.col="black", cex=cexlegendv, seg.len=2)
}

# plotlines5(preds=mlrpred,obs=mlrvgamsmps4,mtxt="MLR")


# DICHOTOMIZED CALIBRATION PLOT WITH POINTS
oplotscatter3 <- function(preds,obs,cexv=0.8,cexaxisv=1.4,cexlabv=1.5,cexlegendv=1.2){
  plot(preds[,2]+preds[,3],fitted(obs)[,2]+fitted(obs)[,3],type="p",pch=1,col="orange",lwd=1,ylab="Observed proportion",xlab="Estimated probability",xlim=0:1,ylim=0:1,cex=cexv,cex.axis=cexaxisv,cex.lab=cexlabv)
  points(preds[,3],fitted(obs)[,3],pch=1,col="red",lwd=1,cex=cexv)
  lines(c(0,1),c(0,1),type="l",col="gray",lty=3,lwd=2) # plot the ideal diagonal line
  legend(0, 1, col=c("gray","orange","red"), c("Ideal", "Cat. 2-3", "Cat. 3"), lty=c(3,NA,NA), lwd=c(2,NA,NA), pch=c(NA,1,1), box.col="black", cex=cexlegendv, seg.len=2)
}
oplotscatter4 <- function(preds,obs,cexv=0.8,cexaxisv=1.4,cexlabv=1.5,cexlegendv=1.2){
  plot(preds[,2]+preds[,3]+preds[,4],fitted(obs)[,2]+fitted(obs)[,3]+fitted(obs)[,4],type="p",pch=1,col="orange",lwd=1,ylab="Observed proportion",xlab="Estimated probability",xlim=0:1,ylim=0:1,cex=cexv,cex.axis=cexaxisv,cex.lab=cexlabv)
  points(preds[,3]+preds[,4],fitted(obs)[,3]+fitted(obs)[,4],pch=1,col="red",lwd=1,cex=cexv)
  points(preds[,4],fitted(obs)[,4],pch=1,col="brown",lwd=1,cex=cexv)
  lines(c(0,1),c(0,1),type="l",col="gray",lty=3,lwd=2) # plot the ideal diagonal line
  legend(0, 1, col=c("gray","orange","red","brown"), c("Ideal", "Cat. 2-4", "Cat. 3-4", "Cat. 4"), lty=c(3,NA,NA,NA), lwd=c(2,NA,NA,NA), pch=c(NA,1,1,1), box.col="black", cex=cexlegendv, seg.len=2)
}
oplotscatter5 <- function(preds,obs,cexv=0.8,cexaxisv=1.4,cexlabv=1.5,cexlegendv=1.2){
  plot(preds[,2]+preds[,3]+preds[,4]+preds[,5],fitted(obs)[,2]+fitted(obs)[,3]+fitted(obs)[,4]+fitted(obs)[,5],type="p",pch=1,col="orange",lwd=1,ylab="Observed proportion",xlab="Estimated probability",xlim=0:1,ylim=0:1,cex=cexv,cex.axis=cexaxisv,cex.lab=cexlabv)
  points(preds[,3]+preds[,4]+preds[,5],fitted(obs)[,3]+fitted(obs)[,4]+fitted(obs)[,5],pch=1,col="red",lwd=1,cex=cexv)
  points(preds[,4]+preds[,5],fitted(obs)[,4]+fitted(obs)[,5],pch=1,col="brown",lwd=1,cex=cexv)
  points(preds[,5],fitted(obs)[,5],pch=1,col="black",lwd=1,cex=cexv)
  lines(c(0,1),c(0,1),type="l",col="gray",lty=3,lwd=2) # plot the ideal diagonal line
  legend(0, 1, col=c("gray","orange","red","brown","black"), c("Ideal", "Cat. 2-5", "Cat. 3-5", "Cat. 4-5", "Cat. 5"), lty=c(3,NA,NA,NA,NA), lwd=c(2,NA,NA,NA,NA), pch=c(NA,1,1,1,1), box.col="black", cex=cexlegendv, seg.len=2)
}


# DICHOTOMOUS CALIBRATION PLOT WITH SMOOTHED LINES
oplotlines3 <- function(preds,obs,cexaxisv=1.4,cexlabv=1.5,cexlegendv=1.2){
  wa2=smooth.spline(preds[,2]+preds[,3],fitted(obs)[,2]+fitted(obs)[,3])
  wa3=smooth.spline(preds[,3],fitted(obs)[,3])
  plot(wa2$x, wa1$y,type="l",col="orange",lwd=3,ylab="Observed proportion",xlab="Estimated probability",xlim=0:1,ylim=0:1,cex.axis=cexaxisv,cex.lab=cexlabv)
  lines(wa3$x, wa3$y,col="red",lwd=3)
  lines(c(0,1),c(0,1),type="l",col="gray",lty=3,lwd=2) # plot the ideal diagonal line
  legend(0, 1, col=c("gray","orange","red"), c("Ideal", "Cat. 2-3","Cat. 3"), lty=c(3,1,1), lwd=c(2,3,3), box.col="black", cex=cexlegendv, seg.len=2)
}
oplotlines4 <- function(preds,obs,cexaxisv=1.4,cexlabv=1.5,cexlegendv=1.2){
  wa2=smooth.spline(preds[,2]+preds[,3]+preds[,4],fitted(obs)[,2]+fitted(obs)[,3]+fitted(obs)[,4])
  wa3=smooth.spline(preds[,3]+preds[,4],fitted(obs)[,3]+fitted(obs)[,4])
  wa4=smooth.spline(preds[,4],fitted(obs)[,4])
  plot(wa2$x, wa2$y,type="l",col="orange",lwd=3,ylab="Observed proportion",xlab="Estimated probability",xlim=0:1,ylim=0:1,cex.axis=cexaxisv,cex.lab=cexlabv)
  lines(wa3$x, wa3$y,col="red",lwd=3)
  lines(wa4$x, wa4$y,col="brown",lwd=3)
  lines(c(0,1),c(0,1),type="l",col="gray",lty=3,lwd=2) # plot the ideal diagonal line
  legend(0, 1, col=c("gray","orange","red","brown"), c("Ideal", "Cat. 2-4","Cat. 3-4", "Cat. 4"), lty=c(3,1,1,1), lwd=c(2,3,3,3), box.col="black", cex=cexlegendv, seg.len=2)
}
oplotlines5 <- function(preds,obs,cexaxisv=1.4,cexlabv=1.5,cexlegendv=1.2){
  wa2=smooth.spline(preds[,2]+preds[,3]+preds[,4]+preds[,5],fitted(obs)[,2]+fitted(obs)[,3]+fitted(obs)[,4]+fitted(obs)[,5])
  wa3=smooth.spline(preds[,3]+preds[,4]+preds[,5],fitted(obs)[,3]+fitted(obs)[,4]+fitted(obs)[,5])
  wa4=smooth.spline(preds[,4]+preds[,5],fitted(obs)[,4]+fitted(obs)[,5])
  wa5=smooth.spline(preds[,5],fitted(obs)[,5])
  plot(wa2$x, wa2$y,type="l",col="orange",lwd=3,ylab="Observed proportion",xlab="Estimated probability",xlim=0:1,ylim=0:1,cex.axis=cexaxisv,cex.lab=cexlabv)
  lines(wa3$x, wa3$y,col="red",lwd=3)
  lines(wa4$x, wa4$y,col="brown",lwd=3)
  lines(wa5$x, wa5$y,col="black",lwd=3)
  lines(c(0,1),c(0,1),type="l",col="gray",lty=3,lwd=2) # plot the ideal diagonal line
  legend(0, 1, col=c("gray","orange","red","brown","black"), c("Ideal", "Cat. 2-5","Cat. 3-5", "Cat. 4-5", "Cat. 5"), lty=c(3,1,1,1,1), lwd=c(2,3,3,3,3), box.col="black", cex=cexlegendv, seg.len=2)
}


# ORDINAL C STATISTIC (ORC)
orc <- function(out,preds,k){
  library(DescTools) # Cstat
  Ec=preds
  for (i in (1:k)){
    Ec[,i]=i*Ec[,i]
  }
  E=rowSums(Ec)
  pwc = rep(NA,k*(k-1)*0.5)
  for (i in (2:k)){
    for (j in (1:(i-1))){
      pwc[((i-1)*(i-2)*0.5)+j]=Cstat(x=E[out==(j) | out==i],resp=out[out==(j) | out==i]==i) # c statistic 1 vs 2
    }
  }
  mean(pwc)
}


# ECI own function (original formula from Van Hoorde et al, J Biomed Inform 2015) - this has a range between 0 (perfect) and 100 (perfect prediction of the wrong categories)
eci_bvc <- function(calout,preds,k){
  (mean((preds-fitted(calout))*(preds-fitted(calout))))*(100*k/2)
}

# ECI own function, compared with random model which is more realistic - this has a range between 0 (perfect) and 1 (similar to random model)
# 
eci_rel <- function(calout,preds,k,outc){
  prevm=matrix((table(outc)/length(outc))[1:k],nrow=dim(preds)[1],ncol=k,byrow=T)
  ecir=mean((preds-prevm)*(preds-prevm))
  ecim=mean((preds-fitted(calout))*(preds-fitted(calout)))
  return(ecim/ecir)
}

# Mean squared prediction error
rmspe <- function(truep,preds){
  return(sqrt((mean((preds-as.matrix(truep))*(preds-as.matrix(truep))))))
}

# Multicategory Brier
mbrier <- function(outc,preds,k){
  outm=1*(outc==1)
  for (j in 2:k){
    outm = cbind(outm,1*(outc==j))
    }
  return((mean((preds-outm)*(preds-outm))))
}