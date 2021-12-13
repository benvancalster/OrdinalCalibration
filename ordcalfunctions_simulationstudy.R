############################################
# SIMULATION STUDY ORDINAL RISK PREDICTION #
############################################

# SELF WRITTEN FUNCTIONS
# AUTHOR: BEN VAN CALSTER


# functions to simulate data with normal prediction distributions conditional on Y (i.e. MLR holds)
getdata_k4p3 <- function(mu_x1,mu_x2,mu_x3,yprev,mult){
  y=sample(1:4, 100*mult, replace=TRUE,prob=yprev/100)
  x1=rnorm(100*mult,0,1)+((y==1)*mu_x1[1]+(y==2)*mu_x1[2]+(y==3)*mu_x1[3]+(y==4)*mu_x1[4])
  x2=rnorm(100*mult,0,1)+((y==1)*mu_x2[1]+(y==2)*mu_x2[2]+(y==3)*mu_x2[3]+(y==4)*mu_x2[4])
  x3=rnorm(100*mult,0,1)+((y==1)*mu_x3[1]+(y==2)*mu_x3[2]+(y==3)*mu_x3[3]+(y==4)*mu_x3[4])
  datax = as.data.frame(cbind(x1,x2,x3,y))
  names(datax)=c("x1","x2","x3","y")
  datax$y = ordered(datax$y)
  pdf1 = dmvnorm(x=datax[,1:3],mean=c(mu_x1[1],mu_x2[1],mu_x3[1]),sigma=diag(3))
  pdf2 = dmvnorm(x=datax[,1:3],mean=c(mu_x1[2],mu_x2[2],mu_x3[2]),sigma=diag(3))
  pdf3 = dmvnorm(x=datax[,1:3],mean=c(mu_x1[3],mu_x2[3],mu_x3[3]),sigma=diag(3))
  pdf4 = dmvnorm(x=datax[,1:3],mean=c(mu_x1[4],mu_x2[4],mu_x3[4]),sigma=diag(3))
  datax$ptrue1 = (yprev[1]*pdf1) / (yprev[1]*pdf1+yprev[2]*pdf2+yprev[3]*pdf3+yprev[4]*pdf4)
  datax$ptrue2 = (yprev[2]*pdf2) / (yprev[1]*pdf1+yprev[2]*pdf2+yprev[3]*pdf3+yprev[4]*pdf4)
  datax$ptrue3 = (yprev[3]*pdf3) / (yprev[1]*pdf1+yprev[2]*pdf2+yprev[3]*pdf3+yprev[4]*pdf4)
  datax$ptrue4 = (yprev[4]*pdf4) / (yprev[1]*pdf1+yprev[2]*pdf2+yprev[3]*pdf3+yprev[4]*pdf4)
  return(list(datax,c(mu_x1,mu_x2,mu_x3,yprev,mult)))
}
getdata_k4p3b <- function(px1,px2,px3,yprev,mult){
  y=sample(1:4, 100*mult, replace=TRUE,prob=yprev/100)
  x1=runif(100*mult,0,1)
  x2=runif(100*mult,0,1)
  x3=runif(100*mult,0,1)
  x1=(y==1)*(x1<=px1[1])+(y==2)*(x1<=px1[2])+(y==3)*(x1<=px1[3])+(y==4)*(x1<=px1[4])
  x2=(y==1)*(x2<=px2[1])+(y==2)*(x2<=px2[2])+(y==3)*(x2<=px2[3])+(y==4)*(x2<=px2[4])
  x3=(y==1)*(x3<=px3[1])+(y==2)*(x3<=px3[2])+(y==3)*(x3<=px3[3])+(y==4)*(x3<=px3[4])
  datax = as.data.frame(cbind(x1,x2,x3,y))
  names(datax)=c("x1","x2","x3","y")
  datax$y = ordered(datax$y)
  pdf1 = ((datax[,1]==1)*px1[1]+(datax[,1]==0)*(1-px1[1]))*((datax[,2]==1)*px2[1]+(datax[,2]==0)*(1-px2[1]))*((datax[,3]==1)*px3[1]+(datax[,3]==0)*(1-px3[1]))
  pdf2 = ((datax[,1]==1)*px1[2]+(datax[,1]==0)*(1-px1[2]))*((datax[,2]==1)*px2[2]+(datax[,2]==0)*(1-px2[2]))*((datax[,3]==1)*px3[2]+(datax[,3]==0)*(1-px3[2]))
  pdf3 = ((datax[,1]==1)*px1[3]+(datax[,1]==0)*(1-px1[3]))*((datax[,2]==1)*px2[3]+(datax[,2]==0)*(1-px2[3]))*((datax[,3]==1)*px3[3]+(datax[,3]==0)*(1-px3[3]))
  pdf4 = ((datax[,1]==1)*px1[4]+(datax[,1]==0)*(1-px1[4]))*((datax[,2]==1)*px2[4]+(datax[,2]==0)*(1-px2[4]))*((datax[,3]==1)*px3[4]+(datax[,3]==0)*(1-px3[4]))
  datax$ptrue1 = (yprev[1]*pdf1) / (yprev[1]*pdf1+yprev[2]*pdf2+yprev[3]*pdf3+yprev[4]*pdf4)
  datax$ptrue2 = (yprev[2]*pdf2) / (yprev[1]*pdf1+yprev[2]*pdf2+yprev[3]*pdf3+yprev[4]*pdf4)
  datax$ptrue3 = (yprev[3]*pdf3) / (yprev[1]*pdf1+yprev[2]*pdf2+yprev[3]*pdf3+yprev[4]*pdf4)
  datax$ptrue4 = (yprev[4]*pdf4) / (yprev[1]*pdf1+yprev[2]*pdf2+yprev[3]*pdf3+yprev[4]*pdf4)
  return(list(datax,c(px1,px2,px3,yprev,mult)))
}
getdata_k3p4 <- function(mu_x1,mu_x2,mu_x3,mu_x4,yprev,mult){
  y=sample(1:3, 100*mult, replace=TRUE,prob=yprev/100)
  x1=rnorm(100*mult,0,1)+((y==1)*mu_x1[1]+(y==2)*mu_x1[2]+(y==3)*mu_x1[3])
  x2=rnorm(100*mult,0,1)+((y==1)*mu_x2[1]+(y==2)*mu_x2[2]+(y==3)*mu_x2[3])
  x3=rnorm(100*mult,0,1)+((y==1)*mu_x3[1]+(y==2)*mu_x3[2]+(y==3)*mu_x3[3])
  x4=rnorm(100*mult,0,1)+((y==1)*mu_x4[1]+(y==2)*mu_x4[2]+(y==3)*mu_x4[3])
  datax = as.data.frame(cbind(x1,x2,x3,x4,y))
  names(datax)=c("x1","x2","x3","x4","y")
  datax$y = ordered(datax$y)
  pdf1 = dmvnorm(x=datax[,1:4],mean=c(mu_x1[1],mu_x2[1],mu_x3[1],mu_x4[1]),sigma=diag(4))
  pdf2 = dmvnorm(x=datax[,1:4],mean=c(mu_x1[2],mu_x2[2],mu_x3[2],mu_x4[2]),sigma=diag(4))
  pdf3 = dmvnorm(x=datax[,1:4],mean=c(mu_x1[3],mu_x2[3],mu_x3[3],mu_x4[3]),sigma=diag(4))
  datax$ptrue1 = (yprev[1]*pdf1) / (yprev[1]*pdf1+yprev[2]*pdf2+yprev[3]*pdf3)
  datax$ptrue2 = (yprev[2]*pdf2) / (yprev[1]*pdf1+yprev[2]*pdf2+yprev[3]*pdf3)
  datax$ptrue3 = (yprev[3]*pdf3) / (yprev[1]*pdf1+yprev[2]*pdf2+yprev[3]*pdf3)
  return(list(datax,c(mu_x1,mu_x2,mu_x3,mu_x4,yprev,mult)))
}
getdata_k3p4b <- function(px1,px2,px3,px4,yprev,mult){
  y=sample(1:3, 100*mult, replace=TRUE,prob=yprev/100)
  x1=runif(100*mult,0,1)
  x2=runif(100*mult,0,1)
  x3=runif(100*mult,0,1)
  x4=runif(100*mult,0,1)
  x1=(y==1)*(x1<=px1[1])+(y==2)*(x1<=px1[2])+(y==3)*(x1<=px1[3])
  x2=(y==1)*(x2<=px2[1])+(y==2)*(x2<=px2[2])+(y==3)*(x2<=px2[3])
  x3=(y==1)*(x3<=px3[1])+(y==2)*(x3<=px3[2])+(y==3)*(x3<=px3[3])
  x4=(y==1)*(x4<=px4[1])+(y==2)*(x4<=px4[2])+(y==3)*(x4<=px4[3])
  datax = as.data.frame(cbind(x1,x2,x3,x4,y))
  names(datax)=c("x1","x2","x3","x4","y")
  datax$y = ordered(datax$y)
  # MAKE PDFS FOR THE BINARY VARIABLES!!!
  pdf1 = ((datax[,1]==1)*px1[1]+(datax[,1]==0)*(1-px1[1]))*((datax[,2]==1)*px2[1]+(datax[,2]==0)*(1-px2[1]))*((datax[,3]==1)*px3[1]+(datax[,3]==0)*(1-px3[1]))*((datax[,4]==1)*px4[1]+(datax[,4]==0)*(1-px4[1]))
  pdf2 = ((datax[,1]==1)*px1[2]+(datax[,1]==0)*(1-px1[2]))*((datax[,2]==1)*px2[2]+(datax[,2]==0)*(1-px2[2]))*((datax[,3]==1)*px3[2]+(datax[,3]==0)*(1-px3[2]))*((datax[,4]==1)*px4[2]+(datax[,4]==0)*(1-px4[2]))
  pdf3 = ((datax[,1]==1)*px1[3]+(datax[,1]==0)*(1-px1[3]))*((datax[,2]==1)*px2[3]+(datax[,2]==0)*(1-px2[3]))*((datax[,3]==1)*px3[3]+(datax[,3]==0)*(1-px3[3]))*((datax[,4]==1)*px4[3]+(datax[,4]==0)*(1-px4[3]))
  datax$ptrue1 = (yprev[1]*pdf1) / (yprev[1]*pdf1+yprev[2]*pdf2+yprev[3]*pdf3)
  datax$ptrue2 = (yprev[2]*pdf2) / (yprev[1]*pdf1+yprev[2]*pdf2+yprev[3]*pdf3)
  datax$ptrue3 = (yprev[3]*pdf3) / (yprev[1]*pdf1+yprev[2]*pdf2+yprev[3]*pdf3)
  return(list(datax,c(px1,px2,px3,px4,yprev,mult)))
}

# function to apply MLR, OLR (cumulative logit with proportional odds), AC-PO, stereotype model
mlrfunc<-function(datafit,dataval,ECI="none",levels=4,modelformula=y ~ x1 + x2 + x3){
  # Fit the model
  acnp <- vglm(modelformula, family=acat(parallel=F), data=datafit)
  # Get predicted probabilities, the linear predictors, and the expected outcome E
  acnppred <- predictvglm(acnp,newdata=dataval,type="response")
  acnplpred <- predictvglm(acnp,newdata=dataval,type="link")
  # Standard Cox recalibration for each outcome
  acnpcalout = calout(out=dataval$y,preds=acnppred,k=levels)
  # Standard Cox recalibration with dichotomized outcomes (2/3 vs 1; and 3 vs 1/2)
  acnpcaldout = caldout(out=dataval$y,preds=acnppred,k=levels)
  # 'Architecture-specific' recalibration: predict outcome using linear predictor and the same architecture used the fit the model itself
  if (levels==4){
    acnprecali <- coefficients(vglm(dataval$y ~ 1, offset = acnplpred[,1:3], family=acat(parallel=F)))[c(1:3)]
    acnprecals <- coefficients(vglm(dataval$y ~ acnplpred[,1] + acnplpred[,2] + acnplpred[,3], constraints=list("(Intercept)"=diag(3),"acnplpred[, 1]"=rbind(1,0,0),"acnplpred[, 2]"=rbind(0,1,0),"acnplpred[, 3]"=rbind(0,0,1)),family=acat(parallel=F)))[c(4:6)]
  }
  if (levels==3){
    acnprecali <- coefficients(vglm(dataval$y ~ 1, offset = acnplpred[,1:2], family=acat(parallel=F)))[c(1:2)]
    acnprecals <- coefficients(vglm(dataval$y ~ acnplpred[,1] + acnplpred[,2], constraints=list("(Intercept)"=diag(2),"acnplpred[, 1]"=rbind(1,0),"acnplpred[, 2]"=rbind(0,1)),family=acat(parallel=F)))[c(3,4)]
  }
  # Overall recalibration based on E
  acnpcale = calE(out=dataval$y,preds=acnppred,k=levels)
  # ordinal c statistic
  acnpc = orc(out=dataval$y,preds=acnppred,k=levels)
  # rMSPE and Brier
  acnpmspe = rmspe(truep=dataval[,c((ncol(dataval)-levels+1):ncol(dataval))],preds=acnppred)
  acnpbr = mbrier(outc=dataval$y,preds=acnppred,k=levels)
  # ECI
  acnplp1=log(acnppred[,2]/acnppred[,1])
  acnplp2=log(acnppred[,3]/acnppred[,1])
  if (levels==4){
    acnplp3=log(acnppred[,4]/acnppred[,1])
  }
  if (ECI=="nonpar"){
    if (levels==4){
      acnpvgamsmps4 = vgam(dataval$y ~ sm.ps(acnplp1,df=4) + sm.ps(acnplp2,df=4) + sm.ps(acnplp3,df=4),family=multinomial(refLevel = 1))
      acnpECI = eci_bvc(calout=acnpvgamsmps4,preds=acnppred,k=levels)
      acnpECIr = eci_rel(calout=acnpvgamsmps4,preds=acnppred,k=levels,outc=dataval$y)
      return(list(c(acnpcalout[,1],acnpcalout[,2],acnpcaldout[,1],acnpcaldout[,2],acnprecali,acnprecals,acnpcale[c(1,2)],acnpmspe,acnpbr,acnpECI,acnpECIr,acnpc),
                  Coefficients(acnp),acnppred,acnpvgamsmps4)) 
    }
    if (levels==3){
      acnpvgamsmps4 = vgam(dataval$y ~ sm.ps(acnplp1,df=4) + sm.ps(acnplp2,df=4),family=multinomial(refLevel = 1))
      acnpECI = eci_bvc(calout=acnpvgamsmps4,preds=acnppred,k=levels)
      acnpECIr = eci_rel(calout=acnpvgamsmps4,preds=acnppred,k=levels,outc=dataval$y)
      return(list(c(acnpcalout[,1],acnpcalout[,2],acnpcaldout[,1],acnpcaldout[,2],acnprecali,acnprecals,acnpcale[c(1,2)],acnpmspe,acnpbr,acnpECI,acnpECIr,acnpc),
                  Coefficients(acnp),acnppred,acnpvgamsmps4)) 
    }
  }
  if (ECI=="rcs3"){
    if (levels==4){
      acnpvgamsmps4 = vglm(dataval$y ~ rcs(log(acnppred[,1]/(1-acnppred[,1])),3) + rcs(log(acnppred[,2]/(1-acnppred[,2])),3) + rcs(log(acnppred[,3]/(1-acnppred[,3])),3),family=multinomial(refLevel = 1))
      acnpECI = eci_bvc(calout=acnpvgamsmps4,preds=acnppred,k=levels)
      acnpECIr = eci_rel(calout=acnpvgamsmps4,preds=acnppred,k=levels,outc=dataval$y)
      return(list(c(acnpcalout[,1],acnpcalout[,2],acnpcaldout[,1],acnpcaldout[,2],acnprecali,acnprecals,acnpcale[c(1,2)],acnpmspe,acnpbr,acnpECI,acnpECIr,acnpc),
                  Coefficients(acnp),acnppred,acnpvgamsmps4)) 
    }
    if (levels==3){
      acnpvgamsmps4 = vglm(dataval$y ~ rcs(log(acnppred[,1]/(1-acnppred[,1])),3) + rcs(log(acnppred[,2]/(1-acnppred[,2])),3),family=multinomial(refLevel = 1))
      acnpECI = eci_bvc(calout=acnpvgamsmps4,preds=acnppred,k=levels)
      acnpECIr = eci_rel(calout=acnpvgamsmps4,preds=acnppred,k=levels,outc=dataval$y)
      return(list(c(acnpcalout[,1],acnpcalout[,2],acnpcaldout[,1],acnpcaldout[,2],acnprecali,acnprecals,acnpcale[c(1,2)],acnpmspe,acnpbr,acnpECI,acnpECIr,acnpc),
                  Coefficients(acnp),acnppred,acnpvgamsmps4)) 
    }
  }
  if(ECI=="none"){
    if (levels==4){
      return(list(c(acnpcalout[,1],acnpcalout[,2],acnpcaldout[,1],acnpcaldout[,2],acnprecali,acnprecals,acnpcale[c(1,2)],acnpmspe,acnpbr,acnpc),
                  Coefficients(acnp),acnppred))
    }
    if (levels==3){
      return(list(c(acnpcalout[,1],acnpcalout[,2],acnpcaldout[,1],acnpcaldout[,2],acnprecali,acnprecals,acnpcale[c(1,2)],acnpmspe,acnpbr,acnpc),
                  Coefficients(acnp),acnppred))
    }
  }
  #acnpcoeff=rbind(Coefficients(acnp)[c(1:3)],Coefficients(acnp)[c(4:6)],Coefficients(acnp)[c(7:9)],Coefficients(acnp)[c(10:12)])
}
olrfunc<-function(datafit,dataval,ECI="none",levels=4,modelformula=y ~ x1 + x2 + x3){
  # Fit the model
  olr <- vglm(modelformula, family=cumulative(parallel=T, reverse=T), data=datafit)
  # Get predicted probabilities, the linear predictors, and the expected outcome E
  olrpred <- predictvglm(olr,newdata=dataval,type="response")
  olrlpred <- predictvglm(olr,newdata=dataval,type="link")
  # Standard Cox recalibration for each outcome
  olrcalout = calout(out=dataval$y,preds=olrpred,k=levels)
  # Standard Cox recalibration with dichotomized outcomes (2/3 vs 1; and 3 vs 1/2)
  olrcaldout = caldout(out=dataval$y,preds=olrpred,k=levels)
  # 'Architecture-specific' recalibration: predict outcome using linear predictor and the same architecture used the fit the model itself
  if (levels==4){
    olrrecali <- coefficients(vglm(dataval$y ~ 1, offset = olrlpred[,1], family=cumulative(parallel=T, reverse=T)))[1]
    olrrecals <- coefficients(vglm(dataval$y ~ olrlpred[,1], family=cumulative(parallel=T, reverse=T)))[c(4)]
    olrrecali=c(olrrecali, coefficients(vglm(dataval$y ~ 1, offset = olrlpred[,2], family=cumulative(parallel=T, reverse=T)))[2])
    olrrecals=c(olrrecals, coefficients(vglm(dataval$y ~ olrlpred[,2], family=cumulative(parallel=T, reverse=T)))[c(4)])
    olrrecali=c(olrrecali, coefficients(vglm(dataval$y ~ 1, offset = olrlpred[,3], family=cumulative(parallel=T, reverse=T)))[3])
    olrrecals=c(olrrecals, coefficients(vglm(dataval$y ~ olrlpred[,3], family=cumulative(parallel=T, reverse=T)))[c(4)])
  }
  if (levels==3){
    olrrecali <- coefficients(vglm(dataval$y ~ 1, offset = olrlpred[,1], family=cumulative(parallel=T, reverse=T)))[1]
    olrrecals <- coefficients(vglm(dataval$y ~ olrlpred[,1], family=cumulative(parallel=T, reverse=T)))[c(3)]
    olrrecali=c(olrrecali, coefficients(vglm(dataval$y ~ 1, offset = olrlpred[,2], family=cumulative(parallel=T, reverse=T)))[2])
    olrrecals=c(olrrecals, coefficients(vglm(dataval$y ~ olrlpred[,2], family=cumulative(parallel=T, reverse=T)))[c(3)])
  }
  # Overall recalibration based on E
  olrcale = calE(out=dataval$y,preds=olrpred,k=levels)
  # ordinal c statistic
  olrc = orc(out=dataval$y,preds=olrpred,k=levels)
  # rMSPE and Brier
  olrmspe = rmspe(truep=dataval[,c((ncol(dataval)-levels+1):ncol(dataval))],preds=olrpred)
  olrbr = mbrier(outc=dataval$y,preds=olrpred,k=levels)
  # ECI  
  olrlp1=log(olrpred[,2]/olrpred[,1])
  olrlp2=log(olrpred[,3]/olrpred[,1])
  if (levels==4){
    olrlp3=log(olrpred[,4]/olrpred[,1])
  }
  if (ECI=="nonpar"){
    if (levels==4){
      olrvgamsmps4 = vgam(dataval$y ~ sm.ps(olrlp1,df=4) + sm.ps(olrlp2,df=4) + sm.ps(olrlp3,df=4),family=multinomial(refLevel = 1))
      olrECI = eci_bvc(calout=olrvgamsmps4,preds=olrpred,k=levels)
      olrECIr = eci_rel(calout=olrvgamsmps4,preds=olrpred,k=levels,outc=dataval$y)
      return(list(c(olrcalout[,1],olrcalout[,2],olrcaldout[,1],olrcaldout[,2],olrrecali,olrrecals,olrcale[c(1,2)],olrmspe,olrbr,olrECI,olrECIr,olrc),
                  Coefficients(olr),olrpred,olrvgamsmps4)) 
    }
    if (levels==3){
      olrvgamsmps4 = vgam(dataval$y ~ sm.ps(olrlp1,df=4) + sm.ps(olrlp2,df=4),family=multinomial(refLevel = 1))
      olrECI = eci_bvc(calout=olrvgamsmps4,preds=olrpred,k=levels)
      olrECIr = eci_rel(calout=olrvgamsmps4,preds=olrpred,k=levels,outc=dataval$y)
      return(list(c(olrcalout[,1],olrcalout[,2],olrcaldout[,1],olrcaldout[,2],olrrecali,olrrecals,olrcale[c(1,2)],olrmspe,olrbr,olrECI,olrECIr,olrc),
                  Coefficients(olr),olrpred,olrvgamsmps4)) 
    }
  }
  if (ECI=="rcs3"){
    if (levels==4){
      olrvgamsmps4 = vglm(dataval$y ~ rcs(log(olrpred[,1]/(1-olrpred[,1])),3) + rcs(log(olrpred[,2]/(1-olrpred[,2])),3) + rcs(log(olrpred[,3]/(1-olrpred[,3])),3),family=multinomial(refLevel = 1))
      olrECI = eci_bvc(calout=olrvgamsmps4,preds=olrpred,k=levels)
      olrECIr = eci_rel(calout=olrvgamsmps4,preds=olrpred,k=levels,outc=dataval$y)
      return(list(c(olrcalout[,1],olrcalout[,2],olrcaldout[,1],olrcaldout[,2],olrrecali,olrrecals,olrcale[c(1,2)],olrmspe,olrbr,olrECI,olrECIr,olrc),
                  Coefficients(olr),olrpred,olrvgamsmps4)) 
    }
    if (levels==3){
      olrvgamsmps4 = vglm(dataval$y ~ rcs(log(olrpred[,1]/(1-olrpred[,1])),3) + rcs(log(olrpred[,2]/(1-olrpred[,2])),3),family=multinomial(refLevel = 1))
      olrECI = eci_bvc(calout=olrvgamsmps4,preds=olrpred,k=levels)
      olrECIr = eci_rel(calout=olrvgamsmps4,preds=olrpred,k=levels,outc=dataval$y)
      return(list(c(olrcalout[,1],olrcalout[,2],olrcaldout[,1],olrcaldout[,2],olrrecali,olrrecals,olrcale[c(1,2)],olrmspe,olrbr,olrECI,olrECIr,olrc),
                  Coefficients(olr),olrpred,olrvgamsmps4)) 
    }
  }
  if (ECI=="none"){
    if (levels==4){
      return(list(c(olrcalout[,1],olrcalout[,2],olrcaldout[,1],olrcaldout[,2],olrrecali,olrrecals,olrcale[c(1,2)],olrmspe,olrbr,olrc),
                  Coefficients(olr),olrpred)) 
    }
    if (levels==3){
      return(list(c(olrcalout[,1],olrcalout[,2],olrcaldout[,1],olrcaldout[,2],olrrecali,olrrecals,olrcale[c(1,2)],olrmspe,olrbr,olrc),
                  Coefficients(olr),olrpred)) 
    }
  }
}
acpfunc<-function(datafit,dataval,ECI="none",levels=4,modelformula=y ~ x1 + x2 + x3){
  # Fit the model
  acp <- vglm(modelformula, family=acat(parallel=T), data=datafit)
  # Get predicted probabilities, the linear predictors, and the expected outcome E
  acppred <- predictvglm(acp,newdata=dataval,type="response")
  acplpred <- predictvglm(acp,newdata=dataval,type="link")
  # Standard Cox recalibration for each outcome
  acpcalout = calout(out=dataval$y,preds=acppred,k=levels)
  # Standard Cox recalibration with dichotomized outcomes (2/3 vs 1; and 3 vs 1/2)
  acpcaldout = caldout(out=dataval$y,preds=acppred,k=levels)
  # 'Architecture-specific' recalibration: predict outcome using linear predictor and the same architecture used the fit the model itself
  if (levels==4){
    acprecali <- coefficients(vglm(dataval$y ~ 1, offset = acplpred[,1], family=acat(parallel=T)))[1]
    acprecals <- coefficients(vglm(dataval$y ~ acplpred[,1], family=acat(parallel=T)))[c(4)]
    acprecali=c(acprecali, coefficients(vglm(dataval$y ~ 1, offset = acplpred[,2], family=acat(parallel=T)))[2])
    acprecals=c(acprecals, coefficients(vglm(dataval$y ~ acplpred[,2], family=acat(parallel=T)))[c(4)])
    acprecali=c(acprecali, coefficients(vglm(dataval$y ~ 1, offset = acplpred[,3], family=acat(parallel=T)))[3])
    acprecals=c(acprecals, coefficients(vglm(dataval$y ~ acplpred[,3], family=acat(parallel=T)))[c(4)])
  }
  if (levels==3){
    acprecali <- coefficients(vglm(dataval$y ~ 1, offset = acplpred[,1], family=acat(parallel=T)))[1]
    acprecals <- coefficients(vglm(dataval$y ~ acplpred[,1], family=acat(parallel=T)))[c(3)]
    acprecali=c(acprecali, coefficients(vglm(dataval$y ~ 1, offset = acplpred[,2], family=acat(parallel=T)))[2])
    acprecals=c(acprecals, coefficients(vglm(dataval$y ~ acplpred[,2], family=acat(parallel=T)))[c(3)])
  }
  # Overall recalibration based on E
  acpcale = calE(out=dataval$y,preds=acppred,k=levels)
  # ordinal c statistic
  acpc = orc(out=dataval$y,preds=acppred,k=levels)
  # rMSPE and Brier
  acpmspe = rmspe(truep=dataval[,c((ncol(dataval)-levels+1):ncol(dataval))],preds=acppred)
  acpbr = mbrier(outc=dataval$y,preds=acppred,k=levels)
  # ECI
  acplp1=log(acppred[,2]/acppred[,1])
  acplp2=log(acppred[,3]/acppred[,1])
  if (levels==4){
    acplp3=log(acppred[,4]/acppred[,1])
  }
  if (ECI=="nonpar"){
    if (levels==4){
      acpvgamsmps4 = vgam(dataval$y ~ sm.ps(acplp1,df=4) + sm.ps(acplp2,df=4) + sm.ps(acplp3,df=4),family=multinomial(refLevel = 1))
      acpECI = eci_bvc(calout=acpvgamsmps4,preds=acppred,k=levels)
      acpECIr = eci_rel(calout=acpvgamsmps4,preds=acppred,k=levels,outc=dataval$y)
      return(list(c(acpcalout[,1],acpcalout[,2],acpcaldout[,1],acpcaldout[,2],acprecali,acprecals,acpcale[c(1,2)],acpmspe,acpbr,acpECI,acpECIr,acpc),
                  Coefficients(acp),acppred,acpvgamsmps4)) 
    }
    if (levels==3){
      acpvgamsmps4 = vgam(dataval$y ~ sm.ps(acplp1,df=4) + sm.ps(acplp2,df=4),family=multinomial(refLevel = 1))
      acpECI = eci_bvc(calout=acpvgamsmps4,preds=acppred,k=levels)
      acpECIr = eci_rel(calout=acpvgamsmps4,preds=acppred,k=levels,outc=dataval$y)
      return(list(c(acpcalout[,1],acpcalout[,2],acpcaldout[,1],acpcaldout[,2],acprecali,acprecals,acpcale[c(1,2)],acpmspe,acpbr,acpECI,acpECIr,acpc),
                  Coefficients(acp),acppred,acpvgamsmps4)) 
    }
  }
  if (ECI=="rcs3"){
    if (levels==4){
      acpvgamsmps4 = vglm(dataval$y ~ rcs(log(acppred[,1]/(1-acppred[,1])),3) + rcs(log(acppred[,2]/(1-acppred[,2])),3) + rcs(log(acppred[,3]/(1-acppred[,3])),3),family=multinomial(refLevel = 1))
      acpECI = eci_bvc(calout=acpvgamsmps4,preds=acppred,k=levels)
      acpECIr = eci_rel(calout=acpvgamsmps4,preds=acppred,k=levels,outc=dataval$y)
      return(list(c(acpcalout[,1],acpcalout[,2],acpcaldout[,1],acpcaldout[,2],acprecali,acprecals,acpcale[c(1,2)],acpmspe,acpbr,acpECI,acpECIr,acpc),
                  Coefficients(acp),acppred,acpvgamsmps4)) 
    }
    if (levels==3){
      acpvgamsmps4 = vglm(dataval$y ~ rcs(log(acppred[,1]/(1-acppred[,1])),3) + rcs(log(acppred[,2]/(1-acppred[,2])),3),family=multinomial(refLevel = 1))
      acpECI = eci_bvc(calout=acpvgamsmps4,preds=acppred,k=levels)
      acpECIr = eci_rel(calout=acpvgamsmps4,preds=acppred,k=levels,outc=dataval$y)
      return(list(c(acpcalout[,1],acpcalout[,2],acpcaldout[,1],acpcaldout[,2],acprecali,acprecals,acpcale[c(1,2)],acpmspe,acpbr,acpECI,acpECIr,acpc),
                  Coefficients(acp),acppred,acpvgamsmps4)) 
    }
  }
  if (ECI=="none"){
    if (levels==4){
      return(list(c(acpcalout[,1],acpcalout[,2],acpcaldout[,1],acpcaldout[,2],acprecali,acprecals,acpcale[c(1,2)],acpmspe,acpbr,acpc),
                  Coefficients(acp),acppred)) 
    }
    if (levels==3){
      return(list(c(acpcalout[,1],acpcalout[,2],acpcaldout[,1],acpcaldout[,2],acprecali,acprecals,acpcale[c(1,2)],acpmspe,acpbr,acpc),
                  Coefficients(acp),acppred)) 
    }
  }
}
smfunc<-function(datafit,dataval,ECI="none",levels=4,modelformula=y ~ x1 + x2 + x3){
  # Fit the model
  sm=rrvglm(modelformula, multinomial(refLevel = 1), data = datafit)
  # Get predicted probabilities, the linear predictors, and the expected outcome E
  smpred <- predictvglm(sm,newdata=dataval,type="response")
  smlpred <- predictvglm(sm,newdata=dataval,type="link")
  # Standard Cox recalibration for each outcome
  smcalout = calout(out=dataval$y,preds=smpred,k=levels)
  # Standard Cox recalibration with dichotomized outcomes (2/3 vs 1; and 3 vs 1/2)
  smcaldout = caldout(out=dataval$y,preds=smpred,k=levels)
  # 'Architecture-specific' recalibration: predict outcome using linear predictor and the same architecture used the fit the model itself
  if (levels==4){
    smrecali <- coefficients(vglm(dataval$y ~ 1, offset = smlpred[,1:3], family=multinomial(refLevel = "1")))[c(1:3)]
    smrecals <- coefficients(vglm(dataval$y ~ smlpred[,1] + smlpred[,2] + smlpred[,3], constraints=list("(Intercept)"=diag(3),"smlpred[, 1]"=rbind(1,0,0),"smlpred[, 2]"=rbind(0,1,0),"smlpred[, 3]"=rbind(0,0,1)),family=multinomial(refLevel = "1")))[c(4:6)]
  }
  if (levels==3){
    smrecali <- coefficients(vglm(dataval$y ~ 1, offset = smlpred[,1:2], family=multinomial(refLevel = "1")))[c(1:2)]
    smrecals <- coefficients(vglm(dataval$y ~ smlpred[,1] + smlpred[,2], constraints=list("(Intercept)"=diag(2),"smlpred[, 1]"=rbind(1,0),"smlpred[, 2]"=rbind(0,1)),family=multinomial(refLevel = "1")))[c(3:4)]
  }
  # Overall recalibration based on E
  smcale = calE(out=dataval$y,preds=smpred,k=levels)
  # ordinal c statistic
  smc = orc(out=dataval$y,preds=smpred,k=levels)
  # rMSPE and Brier
  smmspe = rmspe(truep=dataval[,c((ncol(dataval)-levels+1):ncol(dataval))],preds=smpred)
  smbr = mbrier(outc=dataval$y,preds=smpred,k=levels)
  # ECI
  smlp1=log(smpred[,2]/smpred[,1])
  smlp2=log(smpred[,3]/smpred[,1])
  if (levels==4){
    smlp3=log(smpred[,4]/smpred[,1])
  }
  if (ECI=="nonpar"){
    if (levels==4){
      smvgamsmps4 = vgam(dataval$y ~ sm.ps(smlp1,df=4) + sm.ps(smlp2,df=4) + sm.ps(smlp3,df=4),family=multinomial(refLevel = 1))
      smECI = eci_bvc(calout=smvgamsmps4,preds=smpred,k=levels)
      smECIr = eci_rel(calout=smvgamsmps4,preds=smpred,k=levels,outc=dataval$y)
      return(list(c(smcalout[,1],smcalout[,2],smcaldout[,1],smcaldout[,2],smrecali,smrecals,smcale[c(1,2)],smmspe,smbr,smECI,smECIr,smc),
                  c(constraints(sm)$x1[c(1:(levels-1))],Coefficients(sm)),smpred,smvgamsmps4)) 
    }
    if (levels==3){
      smvgamsmps4 = vgam(dataval$y ~ sm.ps(smlp1,df=4) + sm.ps(smlp2,df=4),family=multinomial(refLevel = 1))
      smECI = eci_bvc(calout=smvgamsmps4,preds=smpred,k=levels)
      smECIr = eci_rel(calout=smvgamsmps4,preds=smpred,k=levels,outc=dataval$y)
      return(list(c(smcalout[,1],smcalout[,2],smcaldout[,1],smcaldout[,2],smrecali,smrecals,smcale[c(1,2)],smmspe,smbr,smECI,smECIr,smc),
                  c(constraints(sm)$x1[c(1:(levels-1))],Coefficients(sm)),smpred,smvgamsmps4)) 
    }
  }
  if (ECI=="rcs3"){
    if (levels==4){
      smvgamsmps4 = vglm(dataval$y ~ rcs(log(smpred[,1]/(1-smpred[,1])),3) + rcs(log(smpred[,2]/(1-smpred[,2])),3) + rcs(log(smpred[,3]/(1-smpred[,3])),3),family=multinomial(refLevel = 1))
      smECI = eci_bvc(calout=smvgamsmps4,preds=smpred,k=levels)
      smECIr = eci_rel(calout=smvgamsmps4,preds=smpred,k=levels,outc=dataval$y)
      return(list(c(smcalout[,1],smcalout[,2],smcaldout[,1],smcaldout[,2],smrecali,smrecals,smcale[c(1,2)],smmspe,smbr,smECI,smECIr,smc),
                  c(constraints(sm)$x1[c(1:(levels-1))],Coefficients(sm)),smpred,smvgamsmps4)) 
    }
    if (levels==3){
      smvgamsmps4 = vglm(dataval$y ~ rcs(log(smpred[,1]/(1-smpred[,1])),3) + rcs(log(smpred[,2]/(1-smpred[,2])),3),family=multinomial(refLevel = 1))
      smECI = eci_bvc(calout=smvgamsmps4,preds=smpred,k=levels)
      smECIr = eci_rel(calout=smvgamsmps4,preds=smpred,k=levels,outc=dataval$y)
      return(list(c(smcalout[,1],smcalout[,2],smcaldout[,1],smcaldout[,2],smrecali,smrecals,smcale[c(1,2)],smmspe,smbr,smECI,smECIr,smc),
                  c(constraints(sm)$x1[c(1:(levels-1))],Coefficients(sm)),smpred,smvgamsmps4)) 
    }
  }
  if (ECI=="none"){
    if (levels==4){
      return(list(c(smcalout[,1],smcalout[,2],smcaldout[,1],smcaldout[,2],smrecali,smrecals,smcale[c(1,2)],smmspe,smbr,smc),
                  c(constraints(sm)$x1[c(1:(levels-1))],Coefficients(sm)),smpred)) 
    }
    if (levels==3){
      return(list(c(smcalout[,1],smcalout[,2],smcaldout[,1],smcaldout[,2],smrecali,smrecals,smcale[c(1,2)],smmspe,smbr,smc),
                  c(constraints(sm)$x1[c(1:(levels-1))],Coefficients(sm)),smpred)) 
    }
  }
}

# quick functions to get calibration plots with subsample
plotlines3sub <- function(preds,obs,mtxt,subsample0=1000,spar1=NULL,spar2=NULL,spar3=NULL,spar4=NULL,cexleg=0.6,cexaxis=1,cexlab=1,showtext=0,showleg=1,cextext=1.5,texttoshow="Model??"){
  wa1=smooth.spline(preds[1:subsample0,1],fitted(obs)[1:subsample0,1],cv=F,spar=spar1)
  wa2=smooth.spline(preds[1:subsample0,2],fitted(obs)[1:subsample0,2],cv=F,spar=spar2)
  wa3=smooth.spline(preds[1:subsample0,3],fitted(obs)[1:subsample0,3],cv=F,spar=spar3)
  ref <- rbind(c(0,0),c(1,1)) # needed to plot the ideal diagonal line
  plot(ref,ref,type="l",col="gray",lwd=2,lty=3,main=mtxt,ylab="Observed proportion",xlab="Estimated probability",xlim=0:1,ylim=0:1,cex.axis=cexaxis,cex.lab=cexlab)
  lines(wa1$x, wa1$y,col="green",lwd=3)
  lines(wa2$x, wa2$y,col="orange",lwd=3)
  lines(wa3$x, wa3$y,col="red",lwd=3)
  if (showleg==1){
    legend(0, 1, col=c("gray","green","orange","red"), c("Ideal", "Cat. 1", "Cat. 2", "Cat. 3"), lty=c(3,1,1,1), lwd=c(2,3,3,3), box.col="black", cex=cexleg, seg.len=3)
  }
  if (showtext==1){
    text(1, 0, texttoshow, adj=c(1,0),cex=cextext)
  }
  #  lines(ref,ref,type="l",col="gray",lwd=2,lty=2) # plot the ideal diagonal line
}
plotlines4sub <- function(preds,obs,mtxt,subsample0=1000,spar1=NULL,spar2=NULL,spar3=NULL,spar4=NULL,cexleg=0.6,cexaxis=1,cexlab=1,showtext=0,showleg=1,cextext=1.5,texttoshow="Model??"){
  wa1=smooth.spline(preds[1:subsample0,1],fitted(obs)[1:subsample0,1],cv=F,spar=spar1)
  wa2=smooth.spline(preds[1:subsample0,2],fitted(obs)[1:subsample0,2],cv=F,spar=spar2)
  wa3=smooth.spline(preds[1:subsample0,3],fitted(obs)[1:subsample0,3],cv=F,spar=spar3)
  wa4=smooth.spline(preds[1:subsample0,4],fitted(obs)[1:subsample0,4],cv=F,spar=spar4)
  ref <- rbind(c(0,0),c(1,1)) # needed to plot the ideal diagonal line
  plot(ref,ref,type="l",col="gray",lwd=2,lty=3,main=mtxt,ylab="Observed proportion",xlab="Estimated probability",xlim=0:1,ylim=0:1,cex.axis=cexaxis,cex.lab=cexlab)
  lines(wa1$x, wa1$y,col="green",lwd=3)
  lines(wa2$x, wa2$y,col="orange",lwd=3)
  lines(wa3$x, wa3$y,col="red",lwd=3)
  lines(wa4$x, wa4$y,col="brown",lwd=3)
  if (showleg==1){
    legend(0, 1, col=c("gray","green","orange","red","brown"), c("Ideal", "Cat. 1", "Cat. 2", "Cat. 3", "Cat. 4"), lty=c(3,1,1,1,1), lwd=c(2,3,3,3,3), box.col="black", cex=cexleg, seg.len=3)
  }
  if (showtext==1){
    text(1, 0, texttoshow, adj=c(1,0),cex=cextext)
  }
  #  lines(ref,ref,type="l",col="gray",lwd=2,lty=2) # plot the ideal diagonal line
}
plotscatter3sub <- function(preds,obs,mtxt,subsample0=1000,cexleg=0.6,cexaxis=1,cexlab=1,showtext=0,showleg=1,cextext=1.5,texttoshow="Model??"){
  ref <- rbind(c(0,0),c(1,1)) # needed to plot the ideal diagonal line
  plot(ref,ref,type="l",col="gray",lwd=2,lty=1,main=mtxt,ylab="Observed proportion",xlab="Estimated probability",xlim=0:1,ylim=0:1,cex.axis=cexaxis,cex.lab=cexlab)
  points(preds[1:subsample0,1],fitted(obs)[1:subsample0,1],col="green",pch=20,lwd=0.5)
  points(preds[1:subsample0,2],fitted(obs)[1:subsample0,2],col="orange",pch=20,lwd=0.5)
  points(preds[1:subsample0,3],fitted(obs)[1:subsample0,3],col="red",pch=20,lwd=0.5)
  if (showleg==1){
    legend(0, 1, col=c("gray","green","orange","red"), c("Ideal", "Cat. 1", "Cat. 2", "Cat. 3"), lty=c(1,NA,NA,NA), lwd=c(2,1,1,1), pch=c(NA,20,20,20), box.col="black", cex=cexleg, seg.len=3)
  }
  if (showtext==1){
    text(1, 0, texttoshow, adj=c(1,0),cex=cextext)
  }
  #  lines(ref,ref,type="l",col="gray",lwd=2,lty=2) # plot the ideal diagonal line
}
plotscatter4sub <- function(preds,obs,mtxt,subsample0=1000,cexleg=0.6,cexaxis=1,cexlab=1,showtext=0,showleg=1,cextext=1.5,texttoshow="Model??"){
  ref <- rbind(c(0,0),c(1,1)) # needed to plot the ideal diagonal line
  plot(ref,ref,type="l",col="gray",lwd=2,lty=1,main=mtxt,ylab="Observed proportion",xlab="Estimated probability",xlim=0:1,ylim=0:1,cex.axis=cexaxis,cex.lab=cexlab)
  points(preds[1:subsample0,1],fitted(obs)[1:subsample0,1],col="green",pch=20,lwd=0.5)
  points(preds[1:subsample0,2],fitted(obs)[1:subsample0,2],col="orange",pch=20,lwd=0.5)
  points(preds[1:subsample0,3],fitted(obs)[1:subsample0,3],col="red",pch=20,lwd=0.5)
  points(preds[1:subsample0,4],fitted(obs)[1:subsample0,4],col="brown",pch=20,lwd=0.5)
  if (showleg==1){
    legend(0, 1, col=c("gray","green","orange","red","brown"), c("Ideal", "Cat. 1", "Cat. 2", "Cat. 3", "Cat. 4"), lty=c(1,NA,NA,NA,NA), lwd=c(2,1,1,1,1), pch=c(NA,20,20,20,20), box.col="black", cex=cexleg, seg.len=3)
  }
  if (showtext==1){
    text(1, 0, texttoshow, adj=c(1,0),cex=cextext)
  }
  #  lines(ref,ref,type="l",col="gray",lwd=2,lty=2) # plot the ideal diagonal line
}
plotlines3sub_dich <- function(preds,obs1,obs2,obs3,mtxt,subsample0=1000,spar1=NULL,spar2=NULL,spar3=NULL,spar4=NULL,cexleg=0.6,cexaxis=1,cexlab=1,showtext=0,showleg=1,cextext=1.5,texttoshow="Model??"){
  wa1=smooth.spline(preds[1:subsample0,1],predict(obs1,type="fitted")[1:subsample0],cv=F,spar=spar1)
  wa2=smooth.spline(preds[1:subsample0,2],predict(obs2,type="fitted")[1:subsample0],cv=F,spar=spar2)
  wa3=smooth.spline(preds[1:subsample0,3],predict(obs3,type="fitted")[1:subsample0],cv=F,spar=spar3)
  ref <- rbind(c(0,0),c(1,1)) # needed to plot the ideal diagonal line
  plot(ref,ref,type="l",col="gray",lwd=2,lty=3,main=mtxt,ylab="Observed proportion",xlab="Estimated probability",xlim=0:1,ylim=0:1,cex.axis=cexaxis,cex.lab=cexlab)
  lines(wa1$x, wa1$y,col="green",lwd=3)
  lines(wa2$x, wa2$y,col="orange",lwd=3)
  lines(wa3$x, wa3$y,col="red",lwd=3)
  if (showleg==1){
    legend(0, 1, col=c("gray","green","orange","red"), c("Ideal", "Cat. 1", "Cat. 2", "Cat. 3"), lty=c(3,1,1,1), lwd=c(2,3,3,3), box.col="black", cex=cexleg, seg.len=3)
  }
  if (showtext==1){
    text(1, 0, texttoshow, adj=c(1,0),cex=cextext)
  }
  #  lines(ref,ref,type="l",col="gray",lwd=2,lty=2) # plot the ideal diagonal line
}
plotlines4sub_dich <- function(preds,obs1,obs2,obs3,obs4,mtxt,subsample0=1000,spar1=NULL,spar2=NULL,spar3=NULL,spar4=NULL,cexleg=0.6,cexaxis=1,cexlab=1,showtext=0,showleg=1,cextext=1.5,texttoshow="Model??"){
  wa1=smooth.spline(preds[1:subsample0,1],predict(obs1,type="fitted")[1:subsample0],cv=F,spar=spar1)
  wa2=smooth.spline(preds[1:subsample0,2],predict(obs2,type="fitted")[1:subsample0],cv=F,spar=spar2)
  wa3=smooth.spline(preds[1:subsample0,3],predict(obs3,type="fitted")[1:subsample0],cv=F,spar=spar3)
  wa4=smooth.spline(preds[1:subsample0,4],predict(obs4,type="fitted")[1:subsample0],cv=F,spar=spar4)
  ref <- rbind(c(0,0),c(1,1)) # needed to plot the ideal diagonal line
  plot(ref,ref,type="l",col="gray",lwd=2,lty=3,main=mtxt,ylab="Observed proportion",xlab="Estimated probability",xlim=0:1,ylim=0:1,cex.axis=cexaxis,cex.lab=cexlab)
  lines(wa1$x, wa1$y,col="green",lwd=3)
  lines(wa2$x, wa2$y,col="orange",lwd=3)
  lines(wa3$x, wa3$y,col="red",lwd=3)
  lines(wa4$x, wa4$y,col="brown",lwd=3)
  if (showleg==1){
    legend(0, 1, col=c("gray","green","orange","red","brown"), c("Ideal", "Cat. 1", "Cat. 2", "Cat. 3", "Cat. 4"), lty=c(3,1,1,1,1), lwd=c(2,3,3,3,3), box.col="black", cex=cexleg, seg.len=3)
  }
  if (showtext==1){
    text(1, 0, texttoshow, adj=c(1,0),cex=cextext)
  }
  #  lines(ref,ref,type="l",col="gray",lwd=2,lty=2) # plot the ideal diagonal line
}
plotscatter3sub_dich <- function(preds,obs1,obs2,obs3,mtxt,subsample0=1000,cexleg=0.6,cexaxis=1,cexlab=1){
  ref <- rbind(c(0,0),c(1,1)) # needed to plot the ideal diagonal line
  plot(ref,ref,type="l",col="gray",lwd=2,lty=1,main=mtxt,ylab="Observed proportion",xlab="Estimated probability",xlim=0:1,ylim=0:1,cex.axis=cexaxis,cex.lab=cexlab)
  points(preds[1:subsample0,1],predict(obs1,type="fitted")[1:subsample0],col="green",pch=20,lwd=0.5)
  points(preds[1:subsample0,2],predict(obs2,type="fitted")[1:subsample0],col="orange",pch=20,lwd=0.5)
  points(preds[1:subsample0,3],predict(obs3,type="fitted")[1:subsample0],col="red",pch=20,lwd=0.5)
  legend(0, 1, col=c("gray","green","orange","red"), c("Ideal", "Cat. 1", "Cat. 2", "Cat. 3"), lty=c(1,NA,NA,NA), lwd=c(2,1,1,1), pch=c(NA,20,20,20), box.col="black", cex=cexleg, seg.len=3)
  #  lines(ref,ref,type="l",col="gray",lwd=2,lty=2) # plot the ideal diagonal line
}
plotscatter4sub_dich <- function(preds,obs1,obs2,obs3,obs4,mtxt,subsample0=1000,cexleg=0.6,cexaxis=1,cexlab=1){
  ref <- rbind(c(0,0),c(1,1)) # needed to plot the ideal diagonal line
  plot(ref,ref,type="l",col="gray",lwd=2,lty=1,main=mtxt,ylab="Observed proportion",xlab="Estimated probability",xlim=0:1,ylim=0:1,cex.axis=cexaxis,cex.lab=cexlab)
  points(preds[1:subsample0,1],predict(obs1,type="fitted")[1:subsample0],col="green",pch=20,lwd=0.5)
  points(preds[1:subsample0,2],predict(obs2,type="fitted")[1:subsample0],col="orange",pch=20,lwd=0.5)
  points(preds[1:subsample0,3],predict(obs3,type="fitted")[1:subsample0],col="red",pch=20,lwd=0.5)
  points(preds[1:subsample0,4],predict(obs4,type="fitted")[1:subsample0],col="brown",pch=20,lwd=0.5)
  legend(0, 1, col=c("gray","green","orange","red","brown"), c("Ideal", "Cat. 1", "Cat. 2", "Cat. 3", "Cat. 4"), lty=c(1,NA,NA,NA,NA), lwd=c(2,1,1,1,1), pch=c(NA,20,20,20,20), box.col="black", cex=cexleg, seg.len=3)
  #  lines(ref,ref,type="l",col="gray",lwd=2,lty=2) # plot the ideal diagonal line
}
plotvstruerisk3 <- function(subsample=1000,cexpts=0.1){
  par(mfrow=c(4,3),
      oma=c(5,5,2,0) + 0.0,
      mar=c(0,0,1,1) + 0.0)
  plot(c(0,1),c(0,1),type="l",col="red",xlim=c(0,1),ylim=c(0,1),cex=cexpts,axes=F)
  points(acnpres[[3]][1:subsample,1],datapo[[1]][1:subsample,6],cex=cexpts)
  axis(side=2,at=c(0,0.2,0.4,0.6,0.8,1))
  box(which="plot")
  text(1, 0, "MLR", adj=c(1,0),cex=1.4)
  text(1, 0.15, "Cat. 1", adj=c(1,0),cex=1.2)
  plot(c(0,1),c(0,1),type="l",col="red",xlim=c(0,1),ylim=c(0,1),cex=cexpts,axes=F)
  points(acnpres[[3]][1:subsample,2],datapo[[1]][1:subsample,7],cex=cexpts)
  box(which="plot")
  text(1, 0, "MLR", adj=c(1,0),cex=1.4)
  text(1, 0.15, "Cat. 2", adj=c(1,0),cex=1.2)
  plot(c(0,1),c(0,1),type="l",col="red",xlim=c(0,1),ylim=c(0,1),cex=cexpts,axes=F)
  points(acnpres[[3]][1:subsample,3],datapo[[1]][1:subsample,8],cex=cexpts)
  box(which="plot")
  text(1, 0, "MLR", adj=c(1,0),cex=1.4)
  text(1, 0.15, "Cat. 3", adj=c(1,0),cex=1.2)
  plot(c(0,1),c(0,1),type="l",col="red",xlim=c(0,1),ylim=c(0,1),cex=cexpts,axes=F)
  points(olrres[[3]][1:subsample,1],datapo[[1]][1:subsample,6],cex=cexpts)
  axis(side=2,at=c(0,0.2,0.4,0.6,0.8,1))
  box(which="plot")
  text(1, 0, "CL-PO", adj=c(1,0),cex=1.4)
  text(1, 0.15, "Cat. 1", adj=c(1,0),cex=1.2)
  plot(c(0,1),c(0,1),type="l",col="red",xlim=c(0,1),ylim=c(0,1),cex=cexpts,axes=F)
  points(olrres[[3]][1:subsample,2],datapo[[1]][1:subsample,7],cex=cexpts)
  box(which="plot")
  text(1, 0, "CL-PO", adj=c(1,0),cex=1.4)
  text(1, 0.15, "Cat. 2", adj=c(1,0),cex=1.2)
  plot(c(0,1),c(0,1),type="l",col="red",xlim=c(0,1),ylim=c(0,1),cex=cexpts,axes=F)
  points(olrres[[3]][1:subsample,3],datapo[[1]][1:subsample,8],cex=cexpts)
  box(which="plot")
  text(1, 0, "CL-PO", adj=c(1,0),cex=1.4)
  text(1, 0.15, "Cat. 3", adj=c(1,0),cex=1.2)
  plot(c(0,1),c(0,1),type="l",col="red",xlim=c(0,1),ylim=c(0,1),cex=cexpts,axes=F)
  points(acpres[[3]][1:subsample,1],datapo[[1]][1:subsample,6],cex=cexpts)
  axis(side=2,at=c(0,0.2,0.4,0.6,0.8,1))
  box(which="plot")
  text(1, 0, "AC-PO", adj=c(1,0),cex=1.4)
  text(1, 0.15, "Cat. 1", adj=c(1,0),cex=1.2)
  plot(c(0,1),c(0,1),type="l",col="red",xlim=c(0,1),ylim=c(0,1),cex=cexpts,axes=F)
  points(acpres[[3]][1:subsample,2],datapo[[1]][1:subsample,7],cex=cexpts)
  box(which="plot")
  text(1, 0, "AC-PO", adj=c(1,0),cex=1.4)
  text(1, 0.15, "Cat. 2", adj=c(1,0),cex=1.2)
  plot(c(0,1),c(0,1),type="l",col="red",xlim=c(0,1),ylim=c(0,1),cex=cexpts,axes=F)
  points(acpres[[3]][1:subsample,3],datapo[[1]][1:subsample,8],cex=cexpts)
  box(which="plot")
  text(1, 0, "AC-PO", adj=c(1,0),cex=1.4)
  text(1, 0.15, "Cat. 3", adj=c(1,0),cex=1.2)
  plot(c(0,1),c(0,1),type="l",col="red",xlim=c(0,1),ylim=c(0,1),cex=cexpts,axes=F)
  points(smres[[3]][1:subsample,1],datapo[[1]][1:subsample,6],cex=cexpts)
  axis(side=2,at=c(0,0.2,0.4,0.6,0.8,1))
  axis(side=1,at=c(0,0.2,0.4,0.6,0.8,1))
  box(which="plot")
  text(1, 0, "SLM", adj=c(1,0),cex=1.4)
  text(1, 0.15, "Cat. 1", adj=c(1,0),cex=1.2)
  plot(c(0,1),c(0,1),type="l",col="red",xlim=c(0,1),ylim=c(0,1),cex=cexpts,axes=F)
  points(smres[[3]][1:subsample,2],datapo[[1]][1:subsample,7],cex=cexpts)
  axis(side=1,at=c(0,0.2,0.4,0.6,0.8,1))
  box(which="plot")
  text(1, 0, "SLM", adj=c(1,0),cex=1.4)
  text(1, 0.15, "Cat. 2", adj=c(1,0),cex=1.2)
  plot(c(0,1),c(0,1),type="l",col="red",xlim=c(0,1),ylim=c(0,1),cex=cexpts,axes=F)
  points(smres[[3]][1:subsample,3],datapo[[1]][1:subsample,8],cex=cexpts)
  axis(side=1,at=c(0,0.2,0.4,0.6,0.8,1))
  box(which="plot")
  text(1, 0, "SLM", adj=c(1,0),cex=1.4)
  text(1, 0.15, "Cat. 3", adj=c(1,0),cex=1.2)
  title(xlab = list("Estimated probability",cex=2),
        ylab = list("True probability",cex=2),
        outer = TRUE, line=3)
  
}
plotvstruerisk4 <- function(subsample=1000,cexpts=0.1){
  par(mfrow=c(4,4),
      oma=c(5,5,2,0) + 0.0,
      mar=c(0,0,1,1) + 0.0)
  plot(c(0,1),c(0,1),type="l",col="red",xlim=c(0,1),ylim=c(0,1),cex=cexpts,axes=F)
  points(acnpres[[3]][1:subsample,1],datapo[[1]][1:subsample,5],cex=cexpts)
  axis(side=2,at=c(0,0.2,0.4,0.6,0.8,1))
  box(which="plot")
  text(1, 0, "MLR", adj=c(1,0),cex=1.4)
  text(1, 0.15, "Cat. 1", adj=c(1,0),cex=1.2)
  plot(c(0,1),c(0,1),type="l",col="red",xlim=c(0,1),ylim=c(0,1),cex=cexpts,axes=F)
  points(acnpres[[3]][1:subsample,2],datapo[[1]][1:subsample,6],cex=cexpts)
  box(which="plot")
  text(1, 0, "MLR", adj=c(1,0),cex=1.4)
  text(1, 0.15, "Cat. 2", adj=c(1,0),cex=1.2)
  plot(c(0,1),c(0,1),type="l",col="red",xlim=c(0,1),ylim=c(0,1),cex=cexpts,axes=F)
  points(acnpres[[3]][1:subsample,3],datapo[[1]][1:subsample,7],cex=cexpts)
  box(which="plot")
  text(1, 0, "MLR", adj=c(1,0),cex=1.4)
  text(1, 0.15, "Cat. 3", adj=c(1,0),cex=1.2)
  plot(c(0,1),c(0,1),type="l",col="red",xlim=c(0,1),ylim=c(0,1),cex=cexpts,axes=F)
  points(acnpres[[3]][1:subsample,4],datapo[[1]][1:subsample,8],cex=cexpts)
  box(which="plot")
  text(1, 0, "MLR", adj=c(1,0),cex=1.4)
  text(1, 0.15, "Cat. 4", adj=c(1,0),cex=1.2)
  plot(c(0,1),c(0,1),type="l",col="red",xlim=c(0,1),ylim=c(0,1),cex=cexpts,axes=F)
  points(olrres[[3]][1:subsample,1],datapo[[1]][1:subsample,5],cex=cexpts)
  axis(side=2,at=c(0,0.2,0.4,0.6,0.8,1))
  box(which="plot")
  text(1, 0, "CL-PO", adj=c(1,0),cex=1.4)
  text(1, 0.15, "Cat. 1", adj=c(1,0),cex=1.2)
  plot(c(0,1),c(0,1),type="l",col="red",xlim=c(0,1),ylim=c(0,1),cex=cexpts,axes=F)
  points(olrres[[3]][1:subsample,2],datapo[[1]][1:subsample,6],cex=cexpts)
  box(which="plot")
  text(1, 0, "CL-PO", adj=c(1,0),cex=1.4)
  text(1, 0.15, "Cat. 2", adj=c(1,0),cex=1.2)
  plot(c(0,1),c(0,1),type="l",col="red",xlim=c(0,1),ylim=c(0,1),cex=cexpts,axes=F)
  points(olrres[[3]][1:subsample,3],datapo[[1]][1:subsample,7],cex=cexpts)
  box(which="plot")
  text(1, 0, "CL-PO", adj=c(1,0),cex=1.4)
  text(1, 0.15, "Cat. 3", adj=c(1,0),cex=1.2)
  plot(c(0,1),c(0,1),type="l",col="red",xlim=c(0,1),ylim=c(0,1),cex=cexpts,axes=F)
  points(olrres[[3]][1:subsample,4],datapo[[1]][1:subsample,8],cex=cexpts)
  box(which="plot")
  text(1, 0, "CL-PO", adj=c(1,0),cex=1.4)
  text(1, 0.15, "Cat. 4", adj=c(1,0),cex=1.2)
  plot(c(0,1),c(0,1),type="l",col="red",xlim=c(0,1),ylim=c(0,1),cex=cexpts,axes=F)
  points(acpres[[3]][1:subsample,1],datapo[[1]][1:subsample,5],cex=cexpts)
  axis(side=2,at=c(0,0.2,0.4,0.6,0.8,1))
  box(which="plot")
  text(1, 0, "AC-PO", adj=c(1,0),cex=1.4)
  text(1, 0.15, "Cat. 1", adj=c(1,0),cex=1.2)
  plot(c(0,1),c(0,1),type="l",col="red",xlim=c(0,1),ylim=c(0,1),cex=cexpts,axes=F)
  points(acpres[[3]][1:subsample,2],datapo[[1]][1:subsample,6],cex=cexpts)
  box(which="plot")
  text(1, 0, "AC-PO", adj=c(1,0),cex=1.4)
  text(1, 0.15, "Cat. 2", adj=c(1,0),cex=1.2)
  plot(c(0,1),c(0,1),type="l",col="red",xlim=c(0,1),ylim=c(0,1),cex=cexpts,axes=F)
  points(acpres[[3]][1:subsample,3],datapo[[1]][1:subsample,7],cex=cexpts)
  box(which="plot")
  text(1, 0, "AC-PO", adj=c(1,0),cex=1.4)
  text(1, 0.15, "Cat. 3", adj=c(1,0),cex=1.2)
  plot(c(0,1),c(0,1),type="l",col="red",xlim=c(0,1),ylim=c(0,1),cex=cexpts,axes=F)
  points(acpres[[3]][1:subsample,4],datapo[[1]][1:subsample,8],cex=cexpts)
  box(which="plot")
  text(1, 0, "AC-PO", adj=c(1,0),cex=1.4)
  text(1, 0.15, "Cat. 4", adj=c(1,0),cex=1.2)
  plot(c(0,1),c(0,1),type="l",col="red",xlim=c(0,1),ylim=c(0,1),cex=cexpts,axes=F)
  points(smres[[3]][1:subsample,1],datapo[[1]][1:subsample,5],cex=cexpts)
  axis(side=2,at=c(0,0.2,0.4,0.6,0.8,1))
  axis(side=1,at=c(0,0.2,0.4,0.6,0.8,1))
  box(which="plot")
  text(1, 0, "SLM", adj=c(1,0),cex=1.4)
  text(1, 0.15, "Cat. 1", adj=c(1,0),cex=1.2)
  plot(c(0,1),c(0,1),type="l",col="red",xlim=c(0,1),ylim=c(0,1),cex=cexpts,axes=F)
  points(smres[[3]][1:subsample,2],datapo[[1]][1:subsample,6],cex=cexpts)
  axis(side=1,at=c(0,0.2,0.4,0.6,0.8,1))
  box(which="plot")
  text(1, 0, "SLM", adj=c(1,0),cex=1.4)
  text(1, 0.15, "Cat. 2", adj=c(1,0),cex=1.2)
  plot(c(0,1),c(0,1),type="l",col="red",xlim=c(0,1),ylim=c(0,1),cex=cexpts,axes=F,xlab="Estimated probability class 3")
  points(smres[[3]][1:subsample,3],datapo[[1]][1:subsample,7],cex=cexpts)
  axis(side=1,at=c(0,0.2,0.4,0.6,0.8,1))
  box(which="plot")
  text(1, 0, "SLM", adj=c(1,0),cex=1.4)
  text(1, 0.15, "Cat. 3", adj=c(1,0),cex=1.2)
  plot(c(0,1),c(0,1),type="l",col="red",xlim=c(0,1),ylim=c(0,1),cex=cexpts,axes=F,xlab="Estimated probability class 4")
  points(smres[[3]][1:subsample,4],datapo[[1]][1:subsample,8],cex=cexpts)
  axis(side=1,at=c(0,0.2,0.4,0.6,0.8,1))
  box(which="plot")
  text(1, 0, "SLM", adj=c(1,0),cex=1.4)
  text(1, 0.15, "Cat. 4", adj=c(1,0),cex=1.2)
  title(xlab = list("Estimated probability",cex=2),
        ylab = list("True probability",cex=2),
        outer = TRUE, line=3)
  
}

# Combine results for large dataset, and save
comsave <- function(filep1){
  res1=cbind(acnpres[[1]],olrres[[1]],acpres[[1]],smres[[1]])
  if (nlevels(datapo[[1]]$y)==4){
    rownames(res1) <- c("ci_cat1","ci_cat2","ci_cat3","ci_cat4","cs_cat1","cs_cat2","cs_cat3","cs_cat4","ci_gt1","ci_gt2","ci_gt3","cs_gt1","cs_gt2","cs_gt3","ci_lp1","ci_lp2","ci_lp3","cs_lp1","cs_lp2","cs_lp3","ci_E","cs_E","rMSPE","Brier","ECI","ECIr","ORC")
  }
  if (nlevels(datapo[[1]]$y)==3){
    rownames(res1) <- c("ci_cat1","ci_cat2","ci_cat3","cs_cat1","cs_cat2","cs_cat3","ci_gt1","ci_gt2","cs_gt1","cs_gt2","ci_lp1","ci_lp2","cs_lp1","cs_lp2","ci_E","cs_E","rMSPE","Brier","ECI","ECIr","ORC")
  }
  colnames(res1) <- c("ACNP","OLR","ACP","SMLR")
  
  # coefficients; 4 levels, 3 predictors
  if (nlevels(datapo[[1]]$y)==4){
    coeffs=cbind(c(c(NA,NA,NA),acnpres[[2]]),
                 c(c(NA,NA,NA),olrres[[2]][1:3],olrres[[2]][c(4,4,4)],olrres[[2]][c(5,5,5)],olrres[[2]][c(6,6,6)]),
                 c(c(NA,NA,NA),acpres[[2]][1:3],acpres[[2]][c(4,4,4)],acpres[[2]][c(5,5,5)],acpres[[2]][c(6,6,6)]),
                 c(smres[[2]][1:6],smres[[2]][c(7,7,7)],smres[[2]][c(8,8,8)],smres[[2]][c(9,9,9)]))
    rownames(coeffs) <- c("Phi1","Phi2","Phi3","Int1","Int2","Int3","b1_x1","b2_x1","b3_x1","b1_x2","b2_x2","b3_x2","b1_x3","b2_x3","b3_x3")
    colnames(coeffs) <- c("ACNP","OLR","ACP","SMLR")
  }
  # coefficients; 3 levels, 4 predictors
  if (nlevels(datapo[[1]]$y)==3 & length(acnpres[[2]]) == 10){
    coeffs=cbind(c(c(NA,NA),acnpres[[2]]),
                 c(c(NA,NA),olrres[[2]][1:2],olrres[[2]][c(3,3)],olrres[[2]][c(4,4)],olrres[[2]][c(5,5)],olrres[[2]][c(6,6)]),
                 c(c(NA,NA),acpres[[2]][1:2],acpres[[2]][c(3,3)],acpres[[2]][c(4,4)],acpres[[2]][c(5,5)],acpres[[2]][c(6,6)]),
                 c(smres[[2]][1:4],smres[[2]][c(5,5)],smres[[2]][c(6,6)],smres[[2]][c(7,7)],smres[[2]][c(8,8)]))
    rownames(coeffs) <- c("Phi1","Phi2","Int1","Int2","b1_x1","b2_x1","b1_x2","b2_x2","b1_x3","b2_x3","b1_x4","b2_x4")
    colnames(coeffs) <- c("ACNP","OLR","ACP","SMLR")
  }
  # coefficients; 3 levels, 4 predictors + 4 noise
  if (nlevels(datapo[[1]]$y)==3 & length(acnpres[[2]]) == 18){
     coeffs=cbind(c(c(NA,NA),acnpres[[2]]),
                 c(c(NA,NA),olrres[[2]][1:2],olrres[[2]][c(3,3)],olrres[[2]][c(4,4)],olrres[[2]][c(5,5)],olrres[[2]][c(6,6)],
                   olrres[[2]][c(7,7)],olrres[[2]][c(8,8)],olrres[[2]][c(9,9)],olrres[[2]][c(10,10)]),
                  c(c(NA,NA),acpres[[2]][1:2],acpres[[2]][c(3,3)],acpres[[2]][c(4,4)],acpres[[2]][c(5,5)],acpres[[2]][c(6,6)],
                    acpres[[2]][c(7,7)],acpres[[2]][c(8,8)],acpres[[2]][c(9,9)],acpres[[2]][c(10,10)]),
                 c(smres[[2]][1:4],smres[[2]][c(5,5)],smres[[2]][c(6,6)],smres[[2]][c(7,7)],smres[[2]][c(8,8)],
                   smres[[2]][c(9,9)],smres[[2]][c(10,10)],smres[[2]][c(11,11)],smres[[2]][c(12,12)]))
     rownames(coeffs) <- c("Phi1","Phi2","Int1","Int2","b1_x1","b2_x1","b1_x2","b2_x2","b1_x3","b2_x3","b1_x4","b2_x4",
                           "b1_x5","b2_x5","b1_x6","b2_x6","b1_x7","b2_x7","b1_x8","b2_x8")
     colnames(coeffs) <- c("ACNP","OLR","ACP","SMLR")
    }
  
  settings=datapo[[2]]
  
  save(res1,coeffs,settings,file=filep1)
}
plotsave <- function(filep2,subsample=1000){
  # Do this to avoid huge file size for the figures: randomly select 1000 cases (if total sample size is > 1000)
  if (nlevels(datapo[[1]]$y)==4){
    ref=cbind(as.matrix(seq(0,1,by=0.01)),as.matrix(seq(0,1,by=0.01)),as.matrix(seq(0,1,by=0.01)),as.matrix(seq(0,1,by=0.01)))
    true2=rbind(ref,as.matrix(datapo[[1]][1:subsample,5:8]))
    acnpres2=rbind(ref,as.matrix(acnpres[[3]][1:subsample,1:4]))
    olrres2=rbind(ref,as.matrix(olrres[[3]][1:subsample,1:4]))
    acpres2=rbind(ref,as.matrix(acpres[[3]][1:subsample,1:4]))
    smres2=rbind(ref,as.matrix(smres[[3]][1:subsample,1:4]))
    type=c(rep("red",length(seq(0,1,by=0.01))),rep("black",length(true2)))
    pdf(filep2)
    #par(mar=c(5,5,4,2))
    pairs(~true2[,1]+acnpres2[,1]+olrres2[,1]+acpres2[,1]+smres2[,1],col=type,upper.panel=NULL,xlim=c(0,1),ylim=c(0,1),cex=0.1,labels=c("True","MLR","CL-PO","AC-PO","SLM"),main="Estimated probability: category 1")
    pairs(~true2[,2]+acnpres2[,2]+olrres2[,2]+acpres2[,2]+smres2[,2],col=type,upper.panel=NULL,xlim=c(0,1),ylim=c(0,1),cex=0.1,labels=c("True","MLR","CL-PO","AC-PO","SLM"),main="Estimated probability: category 2")
    pairs(~true2[,3]+acnpres2[,3]+olrres2[,3]+acpres2[,3]+smres2[,3],col=type,upper.panel=NULL,xlim=c(0,1),ylim=c(0,1),cex=0.1,labels=c("True","MLR","CL-PO","AC-PO","SLM"),main="Estimated probability: category 3")
    pairs(~true2[,4]+acnpres2[,4]+olrres2[,4]+acpres2[,4]+smres2[,4],col=type,upper.panel=NULL,xlim=c(0,1),ylim=c(0,1),cex=0.1,labels=c("True","MLR","CL-PO","AC-PO","SLM"),main="Estimated probability: category 4")
    
    par(mfrow=c(2,2))
    
    plotlines4sub(preds=acnpres[[3]],obs=vgamsmps4a,mtxt="MLR",subsample0=subsample)
    plotlines4sub(preds=olrres[[3]],obs=vgamsmps4b,mtxt="CL-PO",subsample0=subsample)
    plotlines4sub(preds=acpres[[3]],obs=vgamsmps4c,mtxt="AC-PO",subsample0=subsample)
    plotlines4sub(preds=smres[[3]],obs=vgamsmps4d,mtxt="SLM",subsample0=subsample)
    
    par(mfrow=c(4,4),
        oma=c(5,5,2,0) + 0.0,
        mar=c(0,0,1,1) + 0.0)
    plot(c(0,1),c(0,1),type="l",col="red",xlim=c(0,1),ylim=c(0,1),cex=0.1,axes=F)
    points(acnpres[[3]][1:subsample,1],datapo[[1]][1:subsample,5],cex=0.1)
    axis(side=2,at=c(0,0.2,0.4,0.6,0.8,1))
    box(which="plot")
    title("Cat. 1                 Cat. 2                 Cat. 3                Cat. 4  ",cex.main=2,outer=T)
    text(1, 0, "MLR", adj=c(1,0),cex=1.2)
    
    plot(c(0,1),c(0,1),type="l",col="red",xlim=c(0,1),ylim=c(0,1),cex=0.1,axes=F)
    points(acnpres[[3]][1:subsample,2],datapo[[1]][1:subsample,6],cex=0.1)
    box(which="plot")
    text(1, 0, "MLR", adj=c(1,0),cex=1.2)
    
    plot(c(0,1),c(0,1),type="l",col="red",xlim=c(0,1),ylim=c(0,1),cex=0.1,axes=F)
    points(acnpres[[3]][1:subsample,3],datapo[[1]][1:subsample,7],cex=0.1)
    box(which="plot")
    text(1, 0, "MLR", adj=c(1,0),cex=1.2)
    
    plot(c(0,1),c(0,1),type="l",col="red",xlim=c(0,1),ylim=c(0,1),cex=0.1,axes=F)
    points(acnpres[[3]][1:subsample,4],datapo[[1]][1:subsample,8],cex=0.1)
    box(which="plot")
    text(1, 0, "MLR", adj=c(1,0),cex=1.2)
    
    plot(c(0,1),c(0,1),type="l",col="red",xlim=c(0,1),ylim=c(0,1),cex=0.1,axes=F)
    points(olrres[[3]][1:subsample,1],datapo[[1]][1:subsample,5],cex=0.1)
    axis(side=2,at=c(0,0.2,0.4,0.6,0.8,1))
    box(which="plot")
    text(1, 0, "CL-PO", adj=c(1,0),cex=1.2)
    
    plot(c(0,1),c(0,1),type="l",col="red",xlim=c(0,1),ylim=c(0,1),cex=0.1,axes=F)
    points(olrres[[3]][1:subsample,2],datapo[[1]][1:subsample,6],cex=0.1)
    box(which="plot")
    text(1, 0, "CL-PO", adj=c(1,0),cex=1.2)
    
    plot(c(0,1),c(0,1),type="l",col="red",xlim=c(0,1),ylim=c(0,1),cex=0.1,axes=F)
    points(olrres[[3]][1:subsample,3],datapo[[1]][1:subsample,7],cex=0.1)
    box(which="plot")
    text(1, 0, "CL-PO", adj=c(1,0),cex=1.2)
    
    plot(c(0,1),c(0,1),type="l",col="red",xlim=c(0,1),ylim=c(0,1),cex=0.1,axes=F)
    points(olrres[[3]][1:subsample,4],datapo[[1]][1:subsample,8],cex=0.1)
    box(which="plot")
    text(1, 0, "CL-PO", adj=c(1,0),cex=1.2)
    
    plot(c(0,1),c(0,1),type="l",col="red",xlim=c(0,1),ylim=c(0,1),cex=0.1,axes=F)
    points(acpres[[3]][1:subsample,1],datapo[[1]][1:subsample,5],cex=0.1)
    axis(side=2,at=c(0,0.2,0.4,0.6,0.8,1))
    box(which="plot")
    text(1, 0, "AC-PO", adj=c(1,0),cex=1.2)
    
    plot(c(0,1),c(0,1),type="l",col="red",xlim=c(0,1),ylim=c(0,1),cex=0.1,axes=F)
    points(acpres[[3]][1:subsample,2],datapo[[1]][1:subsample,6],cex=0.1)
    box(which="plot")
    text(1, 0, "AC-PO", adj=c(1,0),cex=1.2)
    
    plot(c(0,1),c(0,1),type="l",col="red",xlim=c(0,1),ylim=c(0,1),cex=0.1,axes=F)
    points(acpres[[3]][1:subsample,3],datapo[[1]][1:subsample,7],cex=0.1)
    box(which="plot")
    text(1, 0, "AC-PO", adj=c(1,0),cex=1.2)
    
    plot(c(0,1),c(0,1),type="l",col="red",xlim=c(0,1),ylim=c(0,1),cex=0.1,axes=F)
    points(acpres[[3]][1:subsample,4],datapo[[1]][1:subsample,8],cex=0.1)
    box(which="plot")
    text(1, 0, "AC-PO", adj=c(1,0),cex=1.2)
    
    plot(c(0,1),c(0,1),type="l",col="red",xlim=c(0,1),ylim=c(0,1),cex=0.1,axes=F)
    points(smres[[3]][1:subsample,1],datapo[[1]][1:subsample,5],cex=0.1)
    axis(side=2,at=c(0,0.2,0.4,0.6,0.8,1))
    axis(side=1,at=c(0,0.2,0.4,0.6,0.8,1),title="Estimated probability class 1")
    box(which="plot")
    text(1, 0, "SLM", adj=c(1,0),cex=1.2)
    
    plot(c(0,1),c(0,1),type="l",col="red",xlim=c(0,1),ylim=c(0,1),cex=0.1,axes=F)
    points(smres[[3]][1:subsample,2],datapo[[1]][1:subsample,6],cex=0.1)
    axis(side=1,at=c(0,0.2,0.4,0.6,0.8,1),title="Estimated probability class 2")
    box(which="plot")
    text(1, 0, "SLM", adj=c(1,0),cex=1.2)
    
    plot(c(0,1),c(0,1),type="l",col="red",xlim=c(0,1),ylim=c(0,1),cex=0.1,axes=F,xlab="Estimated probability class 3")
    points(smres[[3]][1:subsample,3],datapo[[1]][1:subsample,7],cex=0.1)
    axis(side=1,at=c(0,0.2,0.4,0.6,0.8,1))
    box(which="plot")
    text(1, 0, "SLM", adj=c(1,0),cex=1.2)
    
    plot(c(0,1),c(0,1),type="l",col="red",xlim=c(0,1),ylim=c(0,1),cex=0.1,axes=F,xlab="Estimated probability class 4")
    points(smres[[3]][1:subsample,4],datapo[[1]][1:subsample,8],cex=0.1)
    axis(side=1,at=c(0,0.2,0.4,0.6,0.8,1))
    box(which="plot")
    text(1, 0, "SLM", adj=c(1,0),cex=1.2)
    
    title(xlab = list("Estimated probability",cex=2),
          ylab = list("True probability",cex=2),
          outer = TRUE, line=3)
    
    dev.off()
  }
  if (nlevels(datapo[[1]]$y)==3){
    ref=cbind(as.matrix(seq(0,1,by=0.01)),as.matrix(seq(0,1,by=0.01)),as.matrix(seq(0,1,by=0.01)))
    true2=rbind(ref,as.matrix(datapo[[1]][1:subsample,6:8]))
    acnpres2=rbind(ref,as.matrix(acnpres[[3]][1:subsample,1:3]))
    olrres2=rbind(ref,as.matrix(olrres[[3]][1:subsample,1:3]))
    acpres2=rbind(ref,as.matrix(acpres[[3]][1:subsample,1:3]))
    smres2=rbind(ref,as.matrix(smres[[3]][1:subsample,1:3]))
    type=c(rep("red",length(seq(0,1,by=0.01))),rep("black",length(true2)))
    pdf(filep2)
    #par(mar=c(5,5,4,2))
    pairs(~true2[,1]+acnpres2[,1]+olrres2[,1]+acpres2[,1]+smres2[,1],col=type,upper.panel=NULL,xlim=c(0,1),ylim=c(0,1),cex=0.1,labels=c("True","MLR","CL-PO","AC-PO","SLM"),main="Estimated probability: category 1")
    pairs(~true2[,2]+acnpres2[,2]+olrres2[,2]+acpres2[,2]+smres2[,2],col=type,upper.panel=NULL,xlim=c(0,1),ylim=c(0,1),cex=0.1,labels=c("True","MLR","CL-PO","AC-PO","SLM"),main="Estimated probability: category 2")
    pairs(~true2[,3]+acnpres2[,3]+olrres2[,3]+acpres2[,3]+smres2[,3],col=type,upper.panel=NULL,xlim=c(0,1),ylim=c(0,1),cex=0.1,labels=c("True","MLR","CL-PO","AC-PO","SLM"),main="Estimated probability: category 3")
    
    par(mfrow=c(2,2))
    
    plotlines3sub(preds=acnpres[[3]],obs=vgamsmps4a,mtxt="MLR",subsample0=subsample)
    plotlines3sub(preds=olrres[[3]],obs=vgamsmps4b,mtxt="CL-PO",subsample0=subsample)
    plotlines3sub(preds=acpres[[3]],obs=vgamsmps4c,mtxt="AC-PO",subsample0=subsample)
    plotlines3sub(preds=smres[[3]],obs=vgamsmps4d,mtxt="SLM",subsample0=subsample)
    
    par(mfrow=c(4,3),
        oma=c(5,5,2,0) + 0.0,
        mar=c(0,0,1,1) + 0.0)
    plot(c(0,1),c(0,1),type="l",col="red",xlim=c(0,1),ylim=c(0,1),cex=0.1,axes=F)
    points(acnpres[[3]][1:subsample,1],datapo[[1]][1:subsample,6],cex=0.1)
    axis(side=2,at=c(0,0.2,0.4,0.6,0.8,1))
    box(which="plot")
    title("Category 1               Category 2                Category 3",cex.main=2,outer=T)
    text(1, 0, "MLR", adj=c(1,0),cex=1.2)
    plot(c(0,1),c(0,1),type="l",col="red",xlim=c(0,1),ylim=c(0,1),cex=0.1,axes=F)
    points(acnpres[[3]][1:subsample,2],datapo[[1]][1:subsample,7],cex=0.1)
    box(which="plot")
    text(1, 0, "MLR", adj=c(1,0),cex=1.2)
    plot(c(0,1),c(0,1),type="l",col="red",xlim=c(0,1),ylim=c(0,1),cex=0.1,axes=F)
    points(acnpres[[3]][1:subsample,3],datapo[[1]][1:subsample,8],cex=0.1)
    box(which="plot")
    text(1, 0, "MLR", adj=c(1,0),cex=1.2)
    plot(c(0,1),c(0,1),type="l",col="red",xlim=c(0,1),ylim=c(0,1),cex=0.1,axes=F)
    points(olrres[[3]][1:subsample,1],datapo[[1]][1:subsample,6],cex=0.1)
    axis(side=2,at=c(0,0.2,0.4,0.6,0.8,1))
    box(which="plot")
    text(1, 0, "CL-PO", adj=c(1,0),cex=1.2)
    plot(c(0,1),c(0,1),type="l",col="red",xlim=c(0,1),ylim=c(0,1),cex=0.1,axes=F)
    points(olrres[[3]][1:subsample,2],datapo[[1]][1:subsample,7],cex=0.1)
    box(which="plot")
    text(1, 0, "CL-PO", adj=c(1,0),cex=1.2)
    plot(c(0,1),c(0,1),type="l",col="red",xlim=c(0,1),ylim=c(0,1),cex=0.1,axes=F)
    points(olrres[[3]][1:subsample,3],datapo[[1]][1:subsample,8],cex=0.1)
    box(which="plot")
    text(1, 0, "CL-PO", adj=c(1,0),cex=1.2)
    plot(c(0,1),c(0,1),type="l",col="red",xlim=c(0,1),ylim=c(0,1),cex=0.1,axes=F)
    points(acpres[[3]][1:subsample,1],datapo[[1]][1:subsample,6],cex=0.1)
    axis(side=2,at=c(0,0.2,0.4,0.6,0.8,1))
    box(which="plot")
    text(1, 0, "AC-PO", adj=c(1,0),cex=1.2)
    plot(c(0,1),c(0,1),type="l",col="red",xlim=c(0,1),ylim=c(0,1),cex=0.1,axes=F)
    points(acpres[[3]][1:subsample,2],datapo[[1]][1:subsample,7],cex=0.1)
    box(which="plot")
    text(1, 0, "AC-PO", adj=c(1,0),cex=1.2)
    plot(c(0,1),c(0,1),type="l",col="red",xlim=c(0,1),ylim=c(0,1),cex=0.1,axes=F)
    points(acpres[[3]][1:subsample,3],datapo[[1]][1:subsample,8],cex=0.1)
    box(which="plot")
    text(1, 0, "AC-PO", adj=c(1,0),cex=1.2)
    plot(c(0,1),c(0,1),type="l",col="red",xlim=c(0,1),ylim=c(0,1),cex=0.1,axes=F)
    points(smres[[3]][1:subsample,1],datapo[[1]][1:subsample,6],cex=0.1)
    axis(side=2,at=c(0,0.2,0.4,0.6,0.8,1))
    axis(side=1,at=c(0,0.2,0.4,0.6,0.8,1),title="Estimated probability class 1")
    box(which="plot")
    text(1, 0, "SLM", adj=c(1,0),cex=1.2)
    plot(c(0,1),c(0,1),type="l",col="red",xlim=c(0,1),ylim=c(0,1),cex=0.1,axes=F)
    points(smres[[3]][1:subsample,2],datapo[[1]][1:subsample,7],cex=0.1)
    axis(side=1,at=c(0,0.2,0.4,0.6,0.8,1),title="Estimated probability class 2")
    box(which="plot")
    text(1, 0, "SLM", adj=c(1,0),cex=1.2)
    plot(c(0,1),c(0,1),type="l",col="red",xlim=c(0,1),ylim=c(0,1),cex=0.1,axes=F,xlab="Estimated probability class 3")
    points(smres[[3]][1:subsample,3],datapo[[1]][1:subsample,8],cex=0.1)
    axis(side=1,at=c(0,0.2,0.4,0.6,0.8,1))
    box(which="plot")
    text(1, 0, "SLM", adj=c(1,0),cex=1.2)
    title(#main=list("Class 1                          Class 2                         Class 3",cex=2),
      xlab = list("Estimated probability",cex=2),
      ylab = list("True probability",cex=2),
      outer = TRUE, line=3)
    
    dev.off()
  }
  
  
  # alternative plot to include all data points
  # int=as.data.frame(cbind(datapo[[1]][,5],acnpres[[3]][,1]))
  # ggplot(int,aes(x=V1,y=V2)) +
  #   ggtitle("Plot of 100K Point Dataset") +
  #   xlab("V1") +
  #   ylab("V2") +
  #   stat_bin_hex(colour="white", na.rm=TRUE) +
  #   scale_fill_gradientn(colours=c("gray","black"),
  #                        name = "Frequency",
  #                        na.value=NA)
}
plotsavebin <- function(filep2,subsample=1000){
  # Do this to avoid huge file size for the figures: randomly select 1000 cases (if total sample size is > 1000)
  if (nlevels(datapo[[1]]$y)==4){
    ref=cbind(as.matrix(seq(0,1,by=0.01)),as.matrix(seq(0,1,by=0.01)),as.matrix(seq(0,1,by=0.01)),as.matrix(seq(0,1,by=0.01)))
    true2=rbind(ref,as.matrix(datapo[[1]][1:subsample,5:8]))
    acnpres2=rbind(ref,as.matrix(acnpres[[3]][1:subsample,1:4]))
    olrres2=rbind(ref,as.matrix(olrres[[3]][1:subsample,1:4]))
    acpres2=rbind(ref,as.matrix(acpres[[3]][1:subsample,1:4]))
    smres2=rbind(ref,as.matrix(smres[[3]][1:subsample,1:4]))
    type=c(rep("red",length(seq(0,1,by=0.01))),rep("black",length(true2)))
    pdf(filep2)
    #par(mar=c(5,5,4,2))
    pairs(~true2[,1]+acnpres2[,1]+olrres2[,1]+acpres2[,1]+smres2[,1],col=type,upper.panel=NULL,xlim=c(0,1),ylim=c(0,1),cex=0.1,labels=c("True","MLR","CL-PO","AC-PO","SLM"),main="Estimated probability: category 1")
    pairs(~true2[,2]+acnpres2[,2]+olrres2[,2]+acpres2[,2]+smres2[,2],col=type,upper.panel=NULL,xlim=c(0,1),ylim=c(0,1),cex=0.1,labels=c("True","MLR","CL-PO","AC-PO","SLM"),main="Estimated probability: category 2")
    pairs(~true2[,3]+acnpres2[,3]+olrres2[,3]+acpres2[,3]+smres2[,3],col=type,upper.panel=NULL,xlim=c(0,1),ylim=c(0,1),cex=0.1,labels=c("True","MLR","CL-PO","AC-PO","SLM"),main="Estimated probability: category 3")
    pairs(~true2[,4]+acnpres2[,4]+olrres2[,4]+acpres2[,4]+smres2[,4],col=type,upper.panel=NULL,xlim=c(0,1),ylim=c(0,1),cex=0.1,labels=c("True","MLR","CL-PO","AC-PO","SLM"),main="Estimated probability: category 4")
    
    par(mfrow=c(2,2))
    
    plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4a,mtxt="MLR",subsample0=subsample)
    plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4b,mtxt="CL-PO",subsample0=subsample)
    plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4c,mtxt="AC-PO",subsample0=subsample)
    plotscatter4sub(preds=smres[[3]],obs=vgamsmps4d,mtxt="SLM",subsample0=subsample)
    
    par(mfrow=c(4,4),
        oma=c(5,5,2,0) + 0.0,
        mar=c(0,0,1,1) + 0.0)
    plot(c(0,1),c(0,1),type="l",col="red",xlim=c(0,1),ylim=c(0,1),cex=0.1,axes=F)
    points(acnpres[[3]][1:subsample,1],datapo[[1]][1:subsample,5],cex=0.1)
    axis(side=2,at=c(0,0.2,0.4,0.6,0.8,1))
    box(which="plot")
    title("Category 1             Category 2             Category 3            Category 4  ",cex.main=2,outer=T)
    text(1, 0, "MLR", adj=c(1,0),cex=1.2)
    
    plot(c(0,1),c(0,1),type="l",col="red",xlim=c(0,1),ylim=c(0,1),cex=0.1,axes=F)
    points(acnpres[[3]][1:subsample,2],datapo[[1]][1:subsample,6],cex=0.1)
    box(which="plot")
    text(1, 0, "MLR", adj=c(1,0),cex=1.2)
    
    plot(c(0,1),c(0,1),type="l",col="red",xlim=c(0,1),ylim=c(0,1),cex=0.1,axes=F)
    points(acnpres[[3]][1:subsample,3],datapo[[1]][1:subsample,7],cex=0.1)
    box(which="plot")
    text(1, 0, "MLR", adj=c(1,0),cex=1.2)
    
    plot(c(0,1),c(0,1),type="l",col="red",xlim=c(0,1),ylim=c(0,1),cex=0.1,axes=F)
    points(acnpres[[3]][1:subsample,4],datapo[[1]][1:subsample,8],cex=0.1)
    box(which="plot")
    text(1, 0, "MLR", adj=c(1,0),cex=1.2)
    
    plot(c(0,1),c(0,1),type="l",col="red",xlim=c(0,1),ylim=c(0,1),cex=0.1,axes=F)
    points(olrres[[3]][1:subsample,1],datapo[[1]][1:subsample,5],cex=0.1)
    axis(side=2,at=c(0,0.2,0.4,0.6,0.8,1))
    box(which="plot")
    text(1, 0, "CL-PO", adj=c(1,0),cex=1.2)
    
    plot(c(0,1),c(0,1),type="l",col="red",xlim=c(0,1),ylim=c(0,1),cex=0.1,axes=F)
    points(olrres[[3]][1:subsample,2],datapo[[1]][1:subsample,6],cex=0.1)
    box(which="plot")
    text(1, 0, "CL-PO", adj=c(1,0),cex=1.2)
    
    plot(c(0,1),c(0,1),type="l",col="red",xlim=c(0,1),ylim=c(0,1),cex=0.1,axes=F)
    points(olrres[[3]][1:subsample,3],datapo[[1]][1:subsample,7],cex=0.1)
    box(which="plot")
    text(1, 0, "CL-PO", adj=c(1,0),cex=1.2)
    
    plot(c(0,1),c(0,1),type="l",col="red",xlim=c(0,1),ylim=c(0,1),cex=0.1,axes=F)
    points(olrres[[3]][1:subsample,4],datapo[[1]][1:subsample,8],cex=0.1)
    box(which="plot")
    text(1, 0, "CL-PO", adj=c(1,0),cex=1.2)
    
    plot(c(0,1),c(0,1),type="l",col="red",xlim=c(0,1),ylim=c(0,1),cex=0.1,axes=F)
    points(acpres[[3]][1:subsample,1],datapo[[1]][1:subsample,5],cex=0.1)
    axis(side=2,at=c(0,0.2,0.4,0.6,0.8,1))
    box(which="plot")
    text(1, 0, "AC-PO", adj=c(1,0),cex=1.2)
    
    plot(c(0,1),c(0,1),type="l",col="red",xlim=c(0,1),ylim=c(0,1),cex=0.1,axes=F)
    points(acpres[[3]][1:subsample,2],datapo[[1]][1:subsample,6],cex=0.1)
    box(which="plot")
    text(1, 0, "AC-PO", adj=c(1,0),cex=1.2)
    
    plot(c(0,1),c(0,1),type="l",col="red",xlim=c(0,1),ylim=c(0,1),cex=0.1,axes=F)
    points(acpres[[3]][1:subsample,3],datapo[[1]][1:subsample,7],cex=0.1)
    box(which="plot")
    text(1, 0, "AC-PO", adj=c(1,0),cex=1.2)
    
    plot(c(0,1),c(0,1),type="l",col="red",xlim=c(0,1),ylim=c(0,1),cex=0.1,axes=F)
    points(acpres[[3]][1:subsample,4],datapo[[1]][1:subsample,8],cex=0.1)
    box(which="plot")
    text(1, 0, "AC-PO", adj=c(1,0),cex=1.2)
    
    plot(c(0,1),c(0,1),type="l",col="red",xlim=c(0,1),ylim=c(0,1),cex=0.1,axes=F)
    points(smres[[3]][1:subsample,1],datapo[[1]][1:subsample,5],cex=0.1)
    axis(side=2,at=c(0,0.2,0.4,0.6,0.8,1))
    axis(side=1,at=c(0,0.2,0.4,0.6,0.8,1),title="Estimated probability class 1")
    box(which="plot")
    text(1, 0, "SLM", adj=c(1,0),cex=1.2)
    
    plot(c(0,1),c(0,1),type="l",col="red",xlim=c(0,1),ylim=c(0,1),cex=0.1,axes=F)
    points(smres[[3]][1:subsample,2],datapo[[1]][1:subsample,6],cex=0.1)
    axis(side=1,at=c(0,0.2,0.4,0.6,0.8,1),title="Estimated probability class 2")
    box(which="plot")
    text(1, 0, "SLM", adj=c(1,0),cex=1.2)
    
    plot(c(0,1),c(0,1),type="l",col="red",xlim=c(0,1),ylim=c(0,1),cex=0.1,axes=F,xlab="Estimated probability class 3")
    points(smres[[3]][1:subsample,3],datapo[[1]][1:subsample,7],cex=0.1)
    axis(side=1,at=c(0,0.2,0.4,0.6,0.8,1))
    box(which="plot")
    text(1, 0, "SLM", adj=c(1,0),cex=1.2)
    
    plot(c(0,1),c(0,1),type="l",col="red",xlim=c(0,1),ylim=c(0,1),cex=0.1,axes=F,xlab="Estimated probability class 3")
    points(smres[[3]][1:subsample,4],datapo[[1]][1:subsample,8],cex=0.1)
    axis(side=1,at=c(0,0.2,0.4,0.6,0.8,1))
    box(which="plot")
    text(1, 0, "SLM", adj=c(1,0),cex=1.2)
    
    title(xlab = list("Estimated probability",cex=2),
          ylab = list("True probability",cex=2),
          outer = TRUE, line=3)
    
    dev.off()
  }
  if (nlevels(datapo[[1]]$y)==3){
    ref=cbind(as.matrix(seq(0,1,by=0.01)),as.matrix(seq(0,1,by=0.01)),as.matrix(seq(0,1,by=0.01)))
    true2=rbind(ref,as.matrix(datapo[[1]][1:subsample,6:8]))
    acnpres2=rbind(ref,as.matrix(acnpres[[3]][1:subsample,1:3]))
    olrres2=rbind(ref,as.matrix(olrres[[3]][1:subsample,1:3]))
    acpres2=rbind(ref,as.matrix(acpres[[3]][1:subsample,1:3]))
    smres2=rbind(ref,as.matrix(smres[[3]][1:subsample,1:3]))
    type=c(rep("red",length(seq(0,1,by=0.01))),rep("black",length(true2)))
    pdf(filep2)
    #par(mar=c(5,5,4,2))
    pairs(~true2[,1]+acnpres2[,1]+olrres2[,1]+acpres2[,1]+smres2[,1],col=type,upper.panel=NULL,xlim=c(0,1),ylim=c(0,1),cex=0.1,labels=c("True","MLR","CL-PO","AC-PO","SLM"),main="Estimated probability: category 1")
    pairs(~true2[,2]+acnpres2[,2]+olrres2[,2]+acpres2[,2]+smres2[,2],col=type,upper.panel=NULL,xlim=c(0,1),ylim=c(0,1),cex=0.1,labels=c("True","MLR","CL-PO","AC-PO","SLM"),main="Estimated probability: category 2")
    pairs(~true2[,3]+acnpres2[,3]+olrres2[,3]+acpres2[,3]+smres2[,3],col=type,upper.panel=NULL,xlim=c(0,1),ylim=c(0,1),cex=0.1,labels=c("True","MLR","CL-PO","AC-PO","SLM"),main="Estimated probability: category 3")
    
    par(mfrow=c(2,2))
    
    plotscatter3sub(preds=acnpres[[3]],obs=vgamsmps4a,mtxt="MLR",subsample0=subsample)
    plotscatter3sub(preds=olrres[[3]],obs=vgamsmps4b,mtxt="CL-PO",subsample0=subsample)
    plotscatter3sub(preds=acpres[[3]],obs=vgamsmps4c,mtxt="AC-PO",subsample0=subsample)
    plotscatter3sub(preds=smres[[3]],obs=vgamsmps4d,mtxt="SLM",subsample0=subsample)
    
    par(mfrow=c(4,3),
        oma=c(5,5,2,0) + 0.0,
        mar=c(0,0,1,1) + 0.0)
    plot(c(0,1),c(0,1),type="l",col="red",xlim=c(0,1),ylim=c(0,1),cex=0.1,axes=F)
    points(acnpres[[3]][1:subsample,1],datapo[[1]][1:subsample,6],cex=0.1)
    axis(side=2,at=c(0,0.2,0.4,0.6,0.8,1))
    box(which="plot")
    title("Category 1                     Category 2                    Category 3   ",cex.main=2,outer=T)
    text(1, 0, "MLR", adj=c(1,0),cex=1.2)
    plot(c(0,1),c(0,1),type="l",col="red",xlim=c(0,1),ylim=c(0,1),cex=0.1,axes=F)
    points(acnpres[[3]][1:subsample,2],datapo[[1]][1:subsample,7],cex=0.1)
    box(which="plot")
    text(1, 0, "MLR", adj=c(1,0),cex=1.2)
    plot(c(0,1),c(0,1),type="l",col="red",xlim=c(0,1),ylim=c(0,1),cex=0.1,axes=F)
    points(acnpres[[3]][1:subsample,3],datapo[[1]][1:subsample,8],cex=0.1)
    box(which="plot")
    text(1, 0, "MLR", adj=c(1,0),cex=1.2)
    plot(c(0,1),c(0,1),type="l",col="red",xlim=c(0,1),ylim=c(0,1),cex=0.1,axes=F)
    points(olrres[[3]][1:subsample,1],datapo[[1]][1:subsample,6],cex=0.1)
    axis(side=2,at=c(0,0.2,0.4,0.6,0.8,1))
    box(which="plot")
    text(1, 0, "CL-PO", adj=c(1,0),cex=1.2)
    plot(c(0,1),c(0,1),type="l",col="red",xlim=c(0,1),ylim=c(0,1),cex=0.1,axes=F)
    points(olrres[[3]][1:subsample,2],datapo[[1]][1:subsample,7],cex=0.1)
    box(which="plot")
    text(1, 0, "CL-PO", adj=c(1,0),cex=1.2)
    plot(c(0,1),c(0,1),type="l",col="red",xlim=c(0,1),ylim=c(0,1),cex=0.1,axes=F)
    points(olrres[[3]][1:subsample,3],datapo[[1]][1:subsample,8],cex=0.1)
    box(which="plot")
    text(1, 0, "CL-PO", adj=c(1,0),cex=1.2)
    plot(c(0,1),c(0,1),type="l",col="red",xlim=c(0,1),ylim=c(0,1),cex=0.1,axes=F)
    points(acpres[[3]][1:subsample,1],datapo[[1]][1:subsample,6],cex=0.1)
    axis(side=2,at=c(0,0.2,0.4,0.6,0.8,1))
    box(which="plot")
    text(1, 0, "AC-PO", adj=c(1,0),cex=1.2)
    plot(c(0,1),c(0,1),type="l",col="red",xlim=c(0,1),ylim=c(0,1),cex=0.1,axes=F)
    points(acpres[[3]][1:subsample,2],datapo[[1]][1:subsample,7],cex=0.1)
    box(which="plot")
    text(1, 0, "AC-PO", adj=c(1,0),cex=1.2)
    plot(c(0,1),c(0,1),type="l",col="red",xlim=c(0,1),ylim=c(0,1),cex=0.1,axes=F)
    points(acpres[[3]][1:subsample,3],datapo[[1]][1:subsample,8],cex=0.1)
    box(which="plot")
    text(1, 0, "AC-PO", adj=c(1,0),cex=1.2)
    plot(c(0,1),c(0,1),type="l",col="red",xlim=c(0,1),ylim=c(0,1),cex=0.1,axes=F)
    points(smres[[3]][1:subsample,1],datapo[[1]][1:subsample,6],cex=0.1)
    axis(side=2,at=c(0,0.2,0.4,0.6,0.8,1))
    axis(side=1,at=c(0,0.2,0.4,0.6,0.8,1),title="Estimated probability class 1")
    box(which="plot")
    text(1, 0, "SLM", adj=c(1,0),cex=1.2)
    plot(c(0,1),c(0,1),type="l",col="red",xlim=c(0,1),ylim=c(0,1),cex=0.1,axes=F)
    points(smres[[3]][1:subsample,2],datapo[[1]][1:subsample,7],cex=0.1)
    axis(side=1,at=c(0,0.2,0.4,0.6,0.8,1),title="Estimated probability class 2")
    box(which="plot")
    text(1, 0, "SLM", adj=c(1,0),cex=1.2)
    plot(c(0,1),c(0,1),type="l",col="red",xlim=c(0,1),ylim=c(0,1),cex=0.1,axes=F,xlab="Estimated probability class 3")
    points(smres[[3]][1:subsample,3],datapo[[1]][1:subsample,8],cex=0.1)
    axis(side=1,at=c(0,0.2,0.4,0.6,0.8,1))
    box(which="plot")
    text(1, 0, "SLM", adj=c(1,0),cex=1.2)
    title(#main=list("Class 1                          Class 2                         Class 3",cex=2),
      xlab = list("Estimated probability",cex=2),
      ylab = list("True probability",cex=2),
      outer = TRUE, line=3)
    
    dev.off()
  }
  
  
  # alternative plot to include all data points
  # int=as.data.frame(cbind(datapo[[1]][,5],acnpres[[3]][,1]))
  # ggplot(int,aes(x=V1,y=V2)) +
  #   ggtitle("Plot of 100K Point Dataset") +
  #   xlab("V1") +
  #   ylab("V2") +
  #   stat_bin_hex(colour="white", na.rm=TRUE) +
  #   scale_fill_gradientn(colours=c("gray","black"),
  #                        name = "Frequency",
  #                        na.value=NA)
}
plotcalsave <- function(filep2x,subsample=1000){ # compare full calibration model with separate binary plots
  if (nlevels(datapo[[1]]$y)==4){
    pdf(filep2x)
    
    par(mfrow=c(2,2))
    plotlines4sub(preds=acnpres[[3]],obs=vgamsmps4a,mtxt="MLR",subsample0=subsample)
    plotlines4sub(preds=olrres[[3]],obs=vgamsmps4b,mtxt="CL-PO",subsample0=subsample)
    plotlines4sub(preds=acpres[[3]],obs=vgamsmps4c,mtxt="AC-PO",subsample0=subsample)
    plotlines4sub(preds=smres[[3]],obs=vgamsmps4d,mtxt="SLM",subsample0=subsample)
    
    par(mfrow=c(2,2))
    plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4a,mtxt="MLR",subsample0=subsample)
    plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4b,mtxt="CL-PO",subsample0=subsample)
    plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4c,mtxt="AC-PO",subsample0=subsample)
    plotscatter4sub(preds=smres[[3]],obs=vgamsmps4d,mtxt="SLM",subsample0=subsample)
    
    par(mfrow=c(2,2))
    plotlines4sub(preds=acnpres[[3]],obs=vgamsmps4rcsa,mtxt="MLR",subsample0=subsample)
    plotlines4sub(preds=olrres[[3]],obs=vgamsmps4rcsb,mtxt="CL-PO",subsample0=subsample)
    plotlines4sub(preds=acpres[[3]],obs=vgamsmps4rcsc,mtxt="AC-PO",subsample0=subsample)
    plotlines4sub(preds=smres[[3]],obs=vgamsmps4rcsd,mtxt="SLM",subsample0=subsample)
    
    par(mfrow=c(2,2))
    plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4rcsa,mtxt="MLR",subsample0=subsample)
    plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4rcsb,mtxt="CL-PO",subsample0=subsample)
    plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4rcsc,mtxt="AC-PO",subsample0=subsample)
    plotscatter4sub(preds=smres[[3]],obs=vgamsmps4rcsd,mtxt="SLM",subsample0=subsample)
    
    par(mfrow=c(2,2))
    plotlines4sub_dich(preds=acnpres[[3]],obs1=dichcala1,obs2=dichcala2,obs3=dichcala3,obs4=dichcala4,mtxt="MLR",subsample0=subsample)
    plotlines4sub_dich(preds=olrres[[3]],obs1=dichcalb1,obs2=dichcalb2,obs3=dichcalb3,obs4=dichcalb4,mtxt="CL-PO",subsample0=subsample)
    plotlines4sub_dich(preds=acpres[[3]],obs1=dichcalc1,obs2=dichcalc2,obs3=dichcalc3,obs4=dichcalc4,mtxt="AC-PO",subsample0=subsample)
    plotlines4sub_dich(preds=smres[[3]],obs1=dichcald1,obs2=dichcald2,obs3=dichcald3,obs4=dichcald4,mtxt="SLM",subsample0=subsample)
    
    par(mfrow=c(2,2))
    plotscatter4sub_dich(preds=acnpres[[3]],obs1=dichcala1,obs2=dichcala2,obs3=dichcala3,obs4=dichcala4,mtxt="MLR",subsample0=subsample)
    plotscatter4sub_dich(preds=olrres[[3]],obs1=dichcalb1,obs2=dichcalb2,obs3=dichcalb3,obs4=dichcalb4,mtxt="CL-PO",subsample0=subsample)
    plotscatter4sub_dich(preds=acpres[[3]],obs1=dichcalc1,obs2=dichcalc2,obs3=dichcalc3,obs4=dichcalc4,mtxt="AC-PO",subsample0=subsample)
    plotscatter4sub_dich(preds=smres[[3]],obs1=dichcald1,obs2=dichcald2,obs3=dichcald3,obs4=dichcald4,mtxt="SLM",subsample0=subsample)
    
    dev.off()
  }
  if (nlevels(datapo[[1]]$y)==3){
    pdf(filep2x)
    
    par(mfrow=c(2,2))
    plotlines3sub(preds=acnpres[[3]],obs=vgamsmps4a,mtxt="MLR",subsample0=subsample)
    plotlines3sub(preds=olrres[[3]],obs=vgamsmps4b,mtxt="CL-PO",subsample0=subsample)
    plotlines3sub(preds=acpres[[3]],obs=vgamsmps4c,mtxt="AC-PO",subsample0=subsample)
    plotlines3sub(preds=smres[[3]],obs=vgamsmps4d,mtxt="SLM",subsample0=subsample)
    
    par(mfrow=c(2,2))
    plotscatter3sub(preds=acnpres[[3]],obs=vgamsmps4a,mtxt="MLR",subsample0=subsample)
    plotscatter3sub(preds=olrres[[3]],obs=vgamsmps4b,mtxt="CL-PO",subsample0=subsample)
    plotscatter3sub(preds=acpres[[3]],obs=vgamsmps4c,mtxt="AC-PO",subsample0=subsample)
    plotscatter3sub(preds=smres[[3]],obs=vgamsmps4d,mtxt="SLM",subsample0=subsample)
    
    par(mfrow=c(2,2))
    plotlines3sub(preds=acnpres[[3]],obs=vgamsmps4rcsa,mtxt="MLR",subsample0=subsample)
    plotlines3sub(preds=olrres[[3]],obs=vgamsmps4rcsb,mtxt="CL-PO",subsample0=subsample)
    plotlines3sub(preds=acpres[[3]],obs=vgamsmps4rcsc,mtxt="AC-PO",subsample0=subsample)
    plotlines3sub(preds=smres[[3]],obs=vgamsmps4rcsd,mtxt="SLM",subsample0=subsample)
    
    par(mfrow=c(2,2))
    plotscatter3sub(preds=acnpres[[3]],obs=vgamsmps4rcsa,mtxt="MLR",subsample0=subsample)
    plotscatter3sub(preds=olrres[[3]],obs=vgamsmps4rcsb,mtxt="CL-PO",subsample0=subsample)
    plotscatter3sub(preds=acpres[[3]],obs=vgamsmps4rcsc,mtxt="AC-PO",subsample0=subsample)
    plotscatter3sub(preds=smres[[3]],obs=vgamsmps4rcsd,mtxt="SLM",subsample0=subsample)
    
    par(mfrow=c(2,2))
    plotlines3sub_dich(preds=acnpres[[3]],obs1=dichcala1,obs2=dichcala2,obs3=dichcala3,mtxt="MLR",subsample0=subsample)
    plotlines3sub_dich(preds=olrres[[3]],obs1=dichcalb1,obs2=dichcalb2,obs3=dichcalb3,mtxt="CL-PO",subsample0=subsample)
    plotlines3sub_dich(preds=acpres[[3]],obs1=dichcalc1,obs2=dichcalc2,obs3=dichcalc3,mtxt="AC-PO",subsample0=subsample)
    plotlines3sub_dich(preds=smres[[3]],obs1=dichcald1,obs2=dichcald2,obs3=dichcald3,mtxt="SLM",subsample0=subsample)
    
    par(mfrow=c(2,2))
    plotscatter3sub_dich(preds=acnpres[[3]],obs1=dichcala1,obs2=dichcala2,obs3=dichcala3,mtxt="MLR",subsample0=subsample)
    plotscatter3sub_dich(preds=olrres[[3]],obs1=dichcalb1,obs2=dichcalb2,obs3=dichcalb3,mtxt="CL-PO",subsample0=subsample)
    plotscatter3sub_dich(preds=acpres[[3]],obs1=dichcalc1,obs2=dichcalc2,obs3=dichcalc3,mtxt="AC-PO",subsample0=subsample)
    plotscatter3sub_dich(preds=smres[[3]],obs1=dichcald1,obs2=dichcald2,obs3=dichcald3,mtxt="SLM",subsample0=subsample)
    
    dev.off()
  }
  
}
# Combine small sample results, and save
comsave_smalln <- function(filep3){
  if (nlevels(datapo[[1]]$y)==4){
    res=cbind(colMeans(mlrsim[,1:27]),colMeans(olrsim[,1:27]),colMeans(acpsim[,1:27]),colMeans(smsim[,1:27]),
              colMedians(mlrsim[,1:27]),colMedians(olrsim[,1:27]),colMedians(acpsim[,1:27]),colMedians(smsim[,1:27]))
    rownames(res) <- c("ci_cat1","ci_cat2","ci_cat3","ci_cat4","cs_cat1","cs_cat2","cs_cat3","cs_cat4","ci_gt1","ci_gt2","ci_gt3","cs_gt1","cs_gt2","cs_gt3","ci_lp1","ci_lp2","ci_lp3","cs_lp1","cs_lp2","cs_lp3","ci_E","cs_E","rMSPE","Brier","ECI","ECIr","ORC")
    colnames(res) <- c("ACNPmean","OLRmean","ACPmean","SMLRmean","ACNPmedian","OLRmedian","ACPmedian","SMLRmedian")
  }
  if (nlevels(datapo[[1]]$y)==3){
    res=cbind(colMeans(mlrsim[,1:21]),colMeans(olrsim[,1:21]),colMeans(acpsim[,1:21]),colMeans(smsim[,1:21]),
              colMedians(mlrsim[,1:21]),colMedians(olrsim[,1:21]),colMedians(acpsim[,1:21]),colMedians(smsim[,1:21]))
    rownames(res) <- c("ci_cat1","ci_cat2","ci_cat3","cs_cat1","cs_cat2","cs_cat3","ci_gt1","ci_gt2","cs_gt1","cs_gt2","ci_lp1","ci_lp2","cs_lp1","cs_lp2","ci_E","cs_E","rMSPE","Brier","ECI","ECIr","ORC")
    colnames(res) <- c("ACNPmean","OLRmean","ACPmean","SMLRmean","ACNPmedian","OLRmedian","ACPmedian","SMLRmedian")
  }
  
  save(mlrsim,olrsim,acpsim,smsim,res,file=filep3)
}

# Compare flexible recalibration approaches
flexrecalk4 <- function(filenm){
  vgamsmps4azz = vgam(datapo[[1]]$y ~ sm.ps(log(acnpres[[3]][,2]/acnpres[[3]][,1]),df=4) + sm.ps(log(acnpres[[3]][,3]/acnpres[[3]][,1]),df=4) + sm.ps(log(acnpres[[3]][,4]/acnpres[[3]][,1]),df=4),family=multinomial(refLevel = 1))
  vgamsmps4bzz = vgam(datapo[[1]]$y ~ sm.ps(log(olrres[[3]][,2]/olrres[[3]][,1]),df=4) + sm.ps(log(olrres[[3]][,3]/olrres[[3]][,1]),df=4) + sm.ps(log(olrres[[3]][,4]/olrres[[3]][,1]),df=4),family=multinomial(refLevel = 1))
  vgamsmps4czz = vgam(datapo[[1]]$y ~ sm.ps(log(acpres[[3]][,2]/acpres[[3]][,1]),df=4) + sm.ps(log(acpres[[3]][,3]/acpres[[3]][,1]),df=4) + sm.ps(log(acpres[[3]][,4]/acpres[[3]][,1]),df=4),family=multinomial(refLevel = 1))
  vgamsmps4dzz = vgam(datapo[[1]]$y ~ sm.ps(log(smres[[3]][,2]/smres[[3]][,1]),df=4) + sm.ps(log(smres[[3]][,3]/smres[[3]][,1]),df=4) + sm.ps(log(smres[[3]][,4]/smres[[3]][,1]),df=4),family=multinomial(refLevel = 1))
  vgamsmps4axx = vgam(datapo[[1]]$y ~ sm.ps(log(acnpres[[3]][,1]/(acnpres[[3]][,2] + acnpres[[3]][,3] + acnpres[[3]][,4])),df=4) + 
                        sm.ps(log((acnpres[[3]][,1] + acnpres[[3]][,2])/(acnpres[[3]][,3] + acnpres[[3]][,4])),df=4) + 
                        sm.ps(log((acnpres[[3]][,1] + acnpres[[3]][,2] + acnpres[[3]][,3])/acnpres[[3]][,4]),df=4),family=multinomial(refLevel = 1))
  vgamsmps4bxx = vgam(datapo[[1]]$y ~ sm.ps(log(olrres[[3]][,1]/(olrres[[3]][,2] + olrres[[3]][,3] + olrres[[3]][,4])),df=4) + 
                        sm.ps(log((olrres[[3]][,1] + olrres[[3]][,2])/(olrres[[3]][,3] + olrres[[3]][,4])),df=4) + 
                        sm.ps(log((olrres[[3]][,1] + olrres[[3]][,2] + olrres[[3]][,3])/olrres[[3]][,4]),df=4),family=multinomial(refLevel = 1))
  vgamsmps4cxx = vgam(datapo[[1]]$y ~ sm.ps(log(acpres[[3]][,1]/(acpres[[3]][,2] + acpres[[3]][,3] + acpres[[3]][,4])),df=4) + 
                        sm.ps(log((acpres[[3]][,1] + acpres[[3]][,2])/(acpres[[3]][,3] + acpres[[3]][,4])),df=4) + 
                        sm.ps(log((acpres[[3]][,1] + acpres[[3]][,2] + acpres[[3]][,3])/acpres[[3]][,4]),df=4),family=multinomial(refLevel = 1))
  vgamsmps4dxx = vgam(datapo[[1]]$y ~ sm.ps(log(smres[[3]][,1]/(smres[[3]][,2] + smres[[3]][,3] + smres[[3]][,4])),df=4) + 
                        sm.ps(log((smres[[3]][,1] + smres[[3]][,2])/(smres[[3]][,3] + smres[[3]][,4])),df=4) + 
                        sm.ps(log((smres[[3]][,1] + smres[[3]][,2] + smres[[3]][,3])/smres[[3]][,4]),df=4),family=multinomial(refLevel = 1))
  vgamsmps4ayy = vgam(datapo[[1]]$y ~ sm.ps(log(acnpres[[3]][,1]/(1 - acnpres[[3]][,1])),df=4) + 
                        sm.ps(log(acnpres[[3]][,2]/(1 - acnpres[[3]][,2])),df=4) + 
                        sm.ps(log(acnpres[[3]][,3]/(1 - acnpres[[3]][,3])),df=4), family=multinomial(refLevel = 1))
  vgamsmps4byy = vgam(datapo[[1]]$y ~ sm.ps(log(olrres[[3]][,1]/(1 - olrres[[3]][,1])),df=4) + 
                        sm.ps(log(olrres[[3]][,2]/(1 - olrres[[3]][,2])),df=4) + 
                        sm.ps(log(olrres[[3]][,3]/(1 - olrres[[3]][,3])),df=4), family=multinomial(refLevel = 1))
  vgamsmps4cyy = vgam(datapo[[1]]$y ~ sm.ps(log(acpres[[3]][,1]/(1 - acpres[[3]][,1])),df=4) + 
                        sm.ps(log(acpres[[3]][,2]/(1 - acpres[[3]][,2])),df=4) + 
                        sm.ps(log(acpres[[3]][,3]/(1 - acpres[[3]][,3])),df=4), family=multinomial(refLevel = 1))
  vgamsmps4dyy = vgam(datapo[[1]]$y ~ sm.ps(log(smres[[3]][,1]/(1 - smres[[3]][,1])),df=4) + 
                        sm.ps(log(smres[[3]][,2]/(1 - smres[[3]][,2])),df=4) + 
                        sm.ps(log(smres[[3]][,3]/(1 - smres[[3]][,3])),df=4), family=multinomial(refLevel = 1))
  
  vgamsmps4azzz = vgam(datapo[[1]]$y ~ sm.ps(log(acnpres[[3]][,2]/acnpres[[3]][,1]),df=4) + sm.ps(log(acnpres[[3]][,3]/acnpres[[3]][,1]),df=4) + sm.ps(log(acnpres[[3]][,4]/acnpres[[3]][,1]),df=4),family=cratio(parallel = F))
  vgamsmps4bzzz = vgam(datapo[[1]]$y ~ sm.ps(log(olrres[[3]][,2]/olrres[[3]][,1]),df=4) + sm.ps(log(olrres[[3]][,3]/olrres[[3]][,1]),df=4) + sm.ps(log(olrres[[3]][,4]/olrres[[3]][,1]),df=4),family=cratio(parallel = F))
  vgamsmps4czzz = vgam(datapo[[1]]$y ~ sm.ps(log(acpres[[3]][,2]/acpres[[3]][,1]),df=4) + sm.ps(log(acpres[[3]][,3]/acpres[[3]][,1]),df=4) + sm.ps(log(acpres[[3]][,4]/acpres[[3]][,1]),df=4),family=cratio(parallel = F))
  vgamsmps4dzzz = vgam(datapo[[1]]$y ~ sm.ps(log(smres[[3]][,2]/smres[[3]][,1]),df=4) + sm.ps(log(smres[[3]][,3]/smres[[3]][,1]),df=4) + sm.ps(log(smres[[3]][,4]/smres[[3]][,1]),df=4),family=cratio(parallel = F))
  vgamsmps4axxx = vgam(datapo[[1]]$y ~ sm.ps(log(acnpres[[3]][,1]/(acnpres[[3]][,2] + acnpres[[3]][,3] + acnpres[[3]][,4])),df=4) + 
                         sm.ps(log((acnpres[[3]][,1] + acnpres[[3]][,2])/(acnpres[[3]][,3] + acnpres[[3]][,4])),df=4) + 
                         sm.ps(log((acnpres[[3]][,1] + acnpres[[3]][,2] + acnpres[[3]][,3])/acnpres[[3]][,4]),df=4),family=cratio(parallel = F))
  vgamsmps4bxxx = vgam(datapo[[1]]$y ~ sm.ps(log(olrres[[3]][,1]/(olrres[[3]][,2] + olrres[[3]][,3] + olrres[[3]][,4])),df=4) + 
                         sm.ps(log((olrres[[3]][,1] + olrres[[3]][,2])/(olrres[[3]][,3] + olrres[[3]][,4])),df=4) + 
                         sm.ps(log((olrres[[3]][,1] + olrres[[3]][,2] + olrres[[3]][,3])/olrres[[3]][,4]),df=4),family=cratio(parallel = F))
  vgamsmps4cxxx = vgam(datapo[[1]]$y ~ sm.ps(log(acpres[[3]][,1]/(acpres[[3]][,2] + acpres[[3]][,3] + acpres[[3]][,4])),df=4) + 
                         sm.ps(log((acpres[[3]][,1] + acpres[[3]][,2])/(acpres[[3]][,3] + acpres[[3]][,4])),df=4) + 
                         sm.ps(log((acpres[[3]][,1] + acpres[[3]][,2] + acpres[[3]][,3])/acpres[[3]][,4]),df=4),family=cratio(parallel = F))
  vgamsmps4dxxx = vgam(datapo[[1]]$y ~ sm.ps(log(smres[[3]][,1]/(smres[[3]][,2] + smres[[3]][,3] + smres[[3]][,4])),df=4) + 
                         sm.ps(log((smres[[3]][,1] + smres[[3]][,2])/(smres[[3]][,3] + smres[[3]][,4])),df=4) + 
                         sm.ps(log((smres[[3]][,1] + smres[[3]][,2] + smres[[3]][,3])/smres[[3]][,4]),df=4),family=cratio(parallel = F))
  vgamsmps4ayyy = vgam(datapo[[1]]$y ~ sm.ps(log(acnpres[[3]][,1]/(1 - acnpres[[3]][,1])),df=4) + 
                         sm.ps(log(acnpres[[3]][,2]/(1 - acnpres[[3]][,2])),df=4) + 
                         sm.ps(log(acnpres[[3]][,3]/(1 - acnpres[[3]][,3])),df=4), family=cratio(parallel = F))
  vgamsmps4byyy = vgam(datapo[[1]]$y ~ sm.ps(log(olrres[[3]][,1]/(1 - olrres[[3]][,1])),df=4) + 
                         sm.ps(log(olrres[[3]][,2]/(1 - olrres[[3]][,2])),df=4) + 
                         sm.ps(log(olrres[[3]][,3]/(1 - olrres[[3]][,3])),df=4), family=cratio(parallel = F))
  vgamsmps4cyyy = vgam(datapo[[1]]$y ~ sm.ps(log(acpres[[3]][,1]/(1 - acpres[[3]][,1])),df=4) + 
                         sm.ps(log(acpres[[3]][,2]/(1 - acpres[[3]][,2])),df=4) + 
                         sm.ps(log(acpres[[3]][,3]/(1 - acpres[[3]][,3])),df=4), family=cratio(parallel = F))
  vgamsmps4dyyy = vgam(datapo[[1]]$y ~ sm.ps(log(smres[[3]][,1]/(1 - smres[[3]][,1])),df=4) + 
                         sm.ps(log(smres[[3]][,2]/(1 - smres[[3]][,2])),df=4) + 
                         sm.ps(log(smres[[3]][,3]/(1 - smres[[3]][,3])),df=4), family=cratio(parallel = F))
  
  ecifrm = c(eci_rel(calout=vgamsmps4azz,preds=acnpres[[3]],k=4,outc=datapo[[1]]$y),
             eci_rel(calout=vgamsmps4bzz,preds=olrres[[3]],k=4,outc=datapo[[1]]$y),
             eci_rel(calout=vgamsmps4czz,preds=acpres[[3]],k=4,outc=datapo[[1]]$y),
             eci_rel(calout=vgamsmps4dzz,preds=smres[[3]],k=4,outc=datapo[[1]]$y),
             eci_rel(calout=vgamsmps4azzz,preds=acnpres[[3]],k=4,outc=datapo[[1]]$y),
             eci_rel(calout=vgamsmps4bzzz,preds=olrres[[3]],k=4,outc=datapo[[1]]$y),
             eci_rel(calout=vgamsmps4czzz,preds=acpres[[3]],k=4,outc=datapo[[1]]$y),
             eci_rel(calout=vgamsmps4dzzz,preds=smres[[3]],k=4,outc=datapo[[1]]$y),
             eci_rel(calout=vgamsmps4axx,preds=acnpres[[3]],k=4,outc=datapo[[1]]$y),
             eci_rel(calout=vgamsmps4bxx,preds=olrres[[3]],k=4,outc=datapo[[1]]$y),
             eci_rel(calout=vgamsmps4cxx,preds=acpres[[3]],k=4,outc=datapo[[1]]$y),
             eci_rel(calout=vgamsmps4dxx,preds=smres[[3]],k=4,outc=datapo[[1]]$y),
             eci_rel(calout=vgamsmps4axxx,preds=acnpres[[3]],k=4,outc=datapo[[1]]$y),
             eci_rel(calout=vgamsmps4bxxx,preds=olrres[[3]],k=4,outc=datapo[[1]]$y),
             eci_rel(calout=vgamsmps4cxxx,preds=acpres[[3]],k=4,outc=datapo[[1]]$y),
             eci_rel(calout=vgamsmps4dxxx,preds=smres[[3]],k=4,outc=datapo[[1]]$y),
             eci_rel(calout=vgamsmps4ayy,preds=acnpres[[3]],k=4,outc=datapo[[1]]$y),
             eci_rel(calout=vgamsmps4byy,preds=olrres[[3]],k=4,outc=datapo[[1]]$y),
             eci_rel(calout=vgamsmps4cyy,preds=acpres[[3]],k=4,outc=datapo[[1]]$y),
             eci_rel(calout=vgamsmps4dyy,preds=smres[[3]],k=4,outc=datapo[[1]]$y),
             eci_rel(calout=vgamsmps4ayyy,preds=acnpres[[3]],k=4,outc=datapo[[1]]$y),
             eci_rel(calout=vgamsmps4byyy,preds=olrres[[3]],k=4,outc=datapo[[1]]$y),
             eci_rel(calout=vgamsmps4cyyy,preds=acpres[[3]],k=4,outc=datapo[[1]]$y),
             eci_rel(calout=vgamsmps4dyyy,preds=smres[[3]],k=4,outc=datapo[[1]]$y))
  
  save(acnpres, olrres, acpres, smres, vgamsmps4azz, vgamsmps4bzz, vgamsmps4czz, vgamsmps4dzz, vgamsmps4azzz, vgamsmps4bzzz, vgamsmps4czzz, vgamsmps4dzzz,
       vgamsmps4axx, vgamsmps4bxx, vgamsmps4cxx, vgamsmps4dxx, vgamsmps4axxx, vgamsmps4bxxx, vgamsmps4cxxx, vgamsmps4dxxx,
       vgamsmps4ayy, vgamsmps4byy, vgamsmps4cyy, vgamsmps4dyy, vgamsmps4ayyy, vgamsmps4byyy, vgamsmps4cyyy, vgamsmps4dyyy,
       ecifrm, file=filenm)
  
}
flexrecalk3 <- function(filenm){
  vgamsmps4azz = vgam(datapo[[1]]$y ~ sm.ps(log(acnpres[[3]][,2]/acnpres[[3]][,1]),df=4) + sm.ps(log(acnpres[[3]][,3]/acnpres[[3]][,1]),df=4),family=multinomial(refLevel = 1))
  vgamsmps4bzz = vgam(datapo[[1]]$y ~ sm.ps(log(olrres[[3]][,2]/olrres[[3]][,1]),df=4) + sm.ps(log(olrres[[3]][,3]/olrres[[3]][,1]),df=4),family=multinomial(refLevel = 1))
  vgamsmps4czz = vgam(datapo[[1]]$y ~ sm.ps(log(acpres[[3]][,2]/acpres[[3]][,1]),df=4) + sm.ps(log(acpres[[3]][,3]/acpres[[3]][,1]),df=4),family=multinomial(refLevel = 1))
  vgamsmps4dzz = vgam(datapo[[1]]$y ~ sm.ps(log(smres[[3]][,2]/smres[[3]][,1]),df=4) + sm.ps(log(smres[[3]][,3]/smres[[3]][,1]),df=4),family=multinomial(refLevel = 1))
  vgamsmps4axx = vgam(datapo[[1]]$y ~ sm.ps(log(acnpres[[3]][,1]/(acnpres[[3]][,2] + acnpres[[3]][,3])),df=4) + 
                        sm.ps(log((acnpres[[3]][,1] + acnpres[[3]][,2])/acnpres[[3]][,3]),df=4),family=multinomial(refLevel = 1))
  vgamsmps4bxx = vgam(datapo[[1]]$y ~ sm.ps(log(olrres[[3]][,1]/(olrres[[3]][,2] + olrres[[3]][,3])),df=4) + 
                        sm.ps(log((olrres[[3]][,1] + olrres[[3]][,2])/olrres[[3]][,3]),df=4),family=multinomial(refLevel = 1))
  vgamsmps4cxx = vgam(datapo[[1]]$y ~ sm.ps(log(acpres[[3]][,1]/(acpres[[3]][,2] + acpres[[3]][,3])),df=4) + 
                        sm.ps(log((acpres[[3]][,1] + acpres[[3]][,2])/acpres[[3]][,3]),df=4),family=multinomial(refLevel = 1))
  vgamsmps4dxx = vgam(datapo[[1]]$y ~ sm.ps(log(smres[[3]][,1]/(smres[[3]][,2] + smres[[3]][,3])),df=4) + 
                        sm.ps(log((smres[[3]][,1] + smres[[3]][,2])/smres[[3]][,3]),df=4),family=multinomial(refLevel = 1))
  vgamsmps4ayy = vgam(datapo[[1]]$y ~ sm.ps(log(acnpres[[3]][,1]/(1 - acnpres[[3]][,1])),df=4) + 
                        sm.ps(log(acnpres[[3]][,2]/(1 - acnpres[[3]][,2])),df=4), family=multinomial(refLevel = 1))
  vgamsmps4byy = vgam(datapo[[1]]$y ~ sm.ps(log(olrres[[3]][,1]/(1 - olrres[[3]][,1])),df=4) + 
                        sm.ps(log(olrres[[3]][,2]/(1 - olrres[[3]][,2])),df=4), family=multinomial(refLevel = 1))
  vgamsmps4cyy = vgam(datapo[[1]]$y ~ sm.ps(log(acpres[[3]][,1]/(1 - acpres[[3]][,1])),df=4) + 
                        sm.ps(log(acpres[[3]][,2]/(1 - acpres[[3]][,2])),df=4), family=multinomial(refLevel = 1))
  vgamsmps4dyy = vgam(datapo[[1]]$y ~ sm.ps(log(smres[[3]][,1]/(1 - smres[[3]][,1])),df=4) + 
                        sm.ps(log(smres[[3]][,2]/(1 - smres[[3]][,2])),df=4), family=multinomial(refLevel = 1))
  
  vgamsmps4azzz = vgam(datapo[[1]]$y ~ sm.ps(log(acnpres[[3]][,2]/acnpres[[3]][,1]),df=4) + sm.ps(log(acnpres[[3]][,3]/acnpres[[3]][,1]),df=4),family=cratio(parallel = F))
  vgamsmps4bzzz = vgam(datapo[[1]]$y ~ sm.ps(log(olrres[[3]][,2]/olrres[[3]][,1]),df=4) + sm.ps(log(olrres[[3]][,3]/olrres[[3]][,1]),df=4),family=cratio(parallel = F))
  vgamsmps4czzz = vgam(datapo[[1]]$y ~ sm.ps(log(acpres[[3]][,2]/acpres[[3]][,1]),df=4) + sm.ps(log(acpres[[3]][,3]/acpres[[3]][,1]),df=4),family=cratio(parallel = F))
  vgamsmps4dzzz = vgam(datapo[[1]]$y ~ sm.ps(log(smres[[3]][,2]/smres[[3]][,1]),df=4) + sm.ps(log(smres[[3]][,3]/smres[[3]][,1]),df=4),family=cratio(parallel = F))
  vgamsmps4axxx = vgam(datapo[[1]]$y ~ sm.ps(log(acnpres[[3]][,1]/(acnpres[[3]][,2] + acnpres[[3]][,3])),df=4) + 
                         sm.ps(log((acnpres[[3]][,1] + acnpres[[3]][,2])/acnpres[[3]][,3]),df=4),family=cratio(parallel = F))
  vgamsmps4bxxx = vgam(datapo[[1]]$y ~ sm.ps(log(olrres[[3]][,1]/(olrres[[3]][,2] + olrres[[3]][,3])),df=4) + 
                         sm.ps(log((olrres[[3]][,1] + olrres[[3]][,2])/olrres[[3]][,3]),df=4),family=cratio(parallel = F))
  vgamsmps4cxxx = vgam(datapo[[1]]$y ~ sm.ps(log(acpres[[3]][,1]/(acpres[[3]][,2] + acpres[[3]][,3])),df=4) + 
                         sm.ps(log((acpres[[3]][,1] + acpres[[3]][,2])/acpres[[3]][,3]),df=4),family=cratio(parallel = F))
  vgamsmps4dxxx = vgam(datapo[[1]]$y ~ sm.ps(log(smres[[3]][,1]/(smres[[3]][,2] + smres[[3]][,3])),df=4) + 
                         sm.ps(log((smres[[3]][,1] + smres[[3]][,2])/smres[[3]][,3]),df=4),family=cratio(parallel = F))
  vgamsmps4ayyy = vgam(datapo[[1]]$y ~ sm.ps(log(acnpres[[3]][,1]/(1 - acnpres[[3]][,1])),df=4) + 
                         sm.ps(log(acnpres[[3]][,2]/(1 - acnpres[[3]][,2])),df=4), family=cratio(parallel = F))
  vgamsmps4byyy = vgam(datapo[[1]]$y ~ sm.ps(log(olrres[[3]][,1]/(1 - olrres[[3]][,1])),df=4) + 
                         sm.ps(log(olrres[[3]][,2]/(1 - olrres[[3]][,2])),df=4), family=cratio(parallel = F))
  vgamsmps4cyyy = vgam(datapo[[1]]$y ~ sm.ps(log(acpres[[3]][,1]/(1 - acpres[[3]][,1])),df=4) + 
                         sm.ps(log(acpres[[3]][,2]/(1 - acpres[[3]][,2])),df=4), family=cratio(parallel = F))
  vgamsmps4dyyy = vgam(datapo[[1]]$y ~ sm.ps(log(smres[[3]][,1]/(1 - smres[[3]][,1])),df=4) + 
                         sm.ps(log(smres[[3]][,2]/(1 - smres[[3]][,2])),df=4), family=cratio(parallel = F))
  
  ecifrm = c(eci_rel(calout=vgamsmps4azz,preds=acnpres[[3]],k=3,outc=datapo[[1]]$y),
             eci_rel(calout=vgamsmps4bzz,preds=olrres[[3]],k=3,outc=datapo[[1]]$y),
             eci_rel(calout=vgamsmps4czz,preds=acpres[[3]],k=3,outc=datapo[[1]]$y),
             eci_rel(calout=vgamsmps4dzz,preds=smres[[3]],k=3,outc=datapo[[1]]$y),
             eci_rel(calout=vgamsmps4azzz,preds=acnpres[[3]],k=3,outc=datapo[[1]]$y),
             eci_rel(calout=vgamsmps4bzzz,preds=olrres[[3]],k=3,outc=datapo[[1]]$y),
             eci_rel(calout=vgamsmps4czzz,preds=acpres[[3]],k=3,outc=datapo[[1]]$y),
             eci_rel(calout=vgamsmps4dzzz,preds=smres[[3]],k=3,outc=datapo[[1]]$y),
             eci_rel(calout=vgamsmps4axx,preds=acnpres[[3]],k=3,outc=datapo[[1]]$y),
             eci_rel(calout=vgamsmps4bxx,preds=olrres[[3]],k=3,outc=datapo[[1]]$y),
             eci_rel(calout=vgamsmps4cxx,preds=acpres[[3]],k=3,outc=datapo[[1]]$y),
             eci_rel(calout=vgamsmps4dxx,preds=smres[[3]],k=3,outc=datapo[[1]]$y),
             eci_rel(calout=vgamsmps4axxx,preds=acnpres[[3]],k=3,outc=datapo[[1]]$y),
             eci_rel(calout=vgamsmps4bxxx,preds=olrres[[3]],k=3,outc=datapo[[1]]$y),
             eci_rel(calout=vgamsmps4cxxx,preds=acpres[[3]],k=3,outc=datapo[[1]]$y),
             eci_rel(calout=vgamsmps4dxxx,preds=smres[[3]],k=3,outc=datapo[[1]]$y),
             eci_rel(calout=vgamsmps4ayy,preds=acnpres[[3]],k=3,outc=datapo[[1]]$y),
             eci_rel(calout=vgamsmps4byy,preds=olrres[[3]],k=3,outc=datapo[[1]]$y),
             eci_rel(calout=vgamsmps4cyy,preds=acpres[[3]],k=3,outc=datapo[[1]]$y),
             eci_rel(calout=vgamsmps4dyy,preds=smres[[3]],k=3,outc=datapo[[1]]$y),
             eci_rel(calout=vgamsmps4ayyy,preds=acnpres[[3]],k=3,outc=datapo[[1]]$y),
             eci_rel(calout=vgamsmps4byyy,preds=olrres[[3]],k=3,outc=datapo[[1]]$y),
             eci_rel(calout=vgamsmps4cyyy,preds=acpres[[3]],k=3,outc=datapo[[1]]$y),
             eci_rel(calout=vgamsmps4dyyy,preds=smres[[3]],k=3,outc=datapo[[1]]$y))
  
  save(acnpres, olrres, acpres, smres, vgamsmps4azz, vgamsmps4bzz, vgamsmps4czz, vgamsmps4dzz, vgamsmps4azzz, vgamsmps4bzzz, vgamsmps4czzz, vgamsmps4dzzz,
       vgamsmps4axx, vgamsmps4bxx, vgamsmps4cxx, vgamsmps4dxx, vgamsmps4axxx, vgamsmps4bxxx, vgamsmps4cxxx, vgamsmps4dxxx,
       vgamsmps4ayy, vgamsmps4byy, vgamsmps4cyy, vgamsmps4dyy, vgamsmps4ayyy, vgamsmps4byyy, vgamsmps4cyyy, vgamsmps4dyyy,
       ecifrm, file=filenm)
  
}
