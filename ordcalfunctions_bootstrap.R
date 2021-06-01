# Bootstrap for MLR
bootmlr <- function(dataset,outc,modelformula,seedval,bootnum,appres,ncat){
  # for some reason, function below does not work without this line, haven't figured out why
bootr=list()
set.seed(seedval) # first set the random number seed, then you can redo the bootstrapping and still get the same results
for (i in 1:bootnum){  
  train_data <- dataset[sample(row.names(dataset),replace=TRUE),]
  # do the modeling
  mb <- vglm(modelformula, family=multinomial(refLevel = "1"), data=train_data)
  
  # predict the values on the bootstrap data
  mpredb <- predictvglm(mb,newdata=train_data,type="response")
  mlpredb <- predictvglm(mb,newdata=train_data,type="link")
  
  # predict the values on the original, unresampled data
  mpredo <- predictvglm(mb,newdata=dataset,type="response")
  mlpredo <- predictvglm(mb,newdata=dataset,type="link")
  
  mcaloutb = calout(out=train_data[,outc],preds=mpredb,k=ncat)
  mcalouto = calout(out=dataset[,outc],preds=mpredo,k=ncat)
  
  mcaldoutb = caldout(out=train_data[,outc],preds=mpredb,k=ncat)
  mcaldouto = caldout(out=dataset[,outc],preds=mpredo,k=ncat)
  
  if(ncat==3){
  mrecalib <- coefficients(vglm(train_data[,outc] ~ 1, offset = mlpredb[,1:2], family=multinomial(refLevel = "1")))[c(1,2)]
  mrecalio <- coefficients(vglm(dataset[,outc] ~ 1, offset = mlpredo[,1:2], family=multinomial(refLevel = "1")))[c(1,2)]
  mrecalsb <- coefficients(vglm(train_data[,outc] ~ mlpredb[,1] + mlpredb[,2], constraints=list("(Intercept)"=diag(2),"mlpredb[, 1]"=rbind(1,0),"mlpredb[, 2]"=rbind(0,1)), family=multinomial(refLevel = "1")))[c(3,4)]
  mrecalso <- coefficients(vglm(dataset[,outc] ~ mlpredo[,1] + mlpredo[,2], constraints=list("(Intercept)"=diag(2),"mlpredo[, 1]"=rbind(1,0),"mlpredo[, 2]"=rbind(0,1)), family=multinomial(refLevel = "1")))[c(3,4)]
  } else {
    if(ncat==4) {
      mrecalib <- coefficients(vglm(train_data[,outc] ~ 1, offset = mlpredb[,1:3], family=multinomial(refLevel = "1")))[c(1,2,3)]
      mrecalio <- coefficients(vglm(dataset[,outc] ~ 1, offset = mlpredo[,1:3], family=multinomial(refLevel = "1")))[c(1,2,3)]
      mrecalsb <- coefficients(vglm(train_data[,outc] ~ mlpredb[,1] + mlpredb[,2] + mlpredb[,3], constraints=list("(Intercept)"=diag(3),"mlpredb[, 1]"=rbind(1,0,0),"mlpredb[, 2]"=rbind(0,1,0),"mlpredb[, 3]"=rbind(0,0,1)), family=multinomial(refLevel = "1")))[c(4:6)]
      mrecalso <- coefficients(vglm(dataset[,outc] ~ mlpredo[,1] + mlpredo[,2] + mlpredo[,3], constraints=list("(Intercept)"=diag(3),"mlpredo[, 1]"=rbind(1,0,0),"mlpredo[, 2]"=rbind(0,1,0),"mlpredo[, 3]"=rbind(0,0,1)), family=multinomial(refLevel = "1")))[c(4:6)]
    } else {
      if(ncat==5) {
        mrecalib <- coefficients(vglm(train_data[,outc] ~ 1, offset = mlpredb[,1:4], family=multinomial(refLevel = "1")))[c(1:4)]
        mrecalio <- coefficients(vglm(dataset[,outc] ~ 1, offset = mlpredo[,1:4], family=multinomial(refLevel = "1")))[c(1:4)]
        mrecalsb <- coefficients(vglm(train_data[,outc] ~ mlpredb[,1] + mlpredb[,2] + mlpredb[,3] + mlpredb[,4], constraints=list("(Intercept)"=diag(4),"mlpredb[, 1]"=rbind(1,0,0,0),"mlpredb[, 2]"=rbind(0,1,0,0),"mlpredb[, 3]"=rbind(0,0,1,0),"mlpredb[, 4]"=rbind(0,0,0,1)), family=multinomial(refLevel = "1")))[c(5:8)]
        mrecalso <- coefficients(vglm(dataset[,outc] ~ mlpredo[,1] + mlpredo[,2] + mlpredo[,3] + mlpredo[,4], constraints=list("(Intercept)"=diag(4),"mlpredo[, 1]"=rbind(1,0,0,0),"mlpredo[, 2]"=rbind(0,1,0,0),"mlpredo[, 3]"=rbind(0,0,1,0),"mlpredo[, 4]"=rbind(0,0,0,1)), family=multinomial(refLevel = "1")))[c(5:8)]
      }
    }
  }

  mcaleb = calE(out=train_data[,outc],preds=mpredb,k=ncat)
  mcaleo = calE(out=dataset[,outc],preds=mpredo,k=ncat)

  # mvgamsmps4bb = vglm(train_data[,outc] ~ rcs(log(mpredb[,1]/(1-mpredb[,1])),3) + rcs(log(mpredb[,2]/(1-mpredb[,2])),3) + rcs(log(mpredb[,3]/(1-mpredb[,3])),3) + rcs(log(mpredb[,4]/(1-mpredb[,4])),3),family=multinomial(refLevel = "1"))
  # mvgamsmps4bo = vglm(dataset[,outc] ~ rcs(log(mpredo[,1]/(1-mpredo[,1])),3) + rcs(log(mpredo[,2]/(1-mpredo[,2])),3) + rcs(log(mpredo[,3]/(1-mpredo[,3])),3) + rcs(log(mpredo[,4]/(1-mpredo[,4])),3),family=multinomial(refLevel = "1"))
  # mECIrbb = eci_rel(calout=mvgamsmps4bb,preds=mpredb,k=5,outc=train_data[,outc])
  # mECIrbo = eci_rel(calout=mvgamsmps4bo,preds=mpredo,k=5,outc=dataset[,outc])
  
  mcb = orc(out=train_data[,outc],preds=mpredb,k=ncat)
  mco = orc(out=dataset[,outc],preds=mpredo,k=ncat)

  bootr[[i]] <- c(mcaloutb[,1]-mcalouto[,1],mcaloutb[,2]-mcalouto[,2],mcaldoutb[,1]-mcaldouto[,1],mcaldoutb[,2]-mcaldouto[,2],mrecalib-mrecalio,mrecalsb-mrecalso,mcaleb-mcaleo,0,mcb-mco)#mECIrbb-mECIrbo,mcb-mco)

  print(i)
  
  }

cbind(appres,appres-colMeans(matrix(unlist(bootr),nrow=bootnum,ncol=length(unlist(bootr))/bootnum,byrow=T)), colMeans(matrix(unlist(bootr),nrow=bootnum,ncol=length(unlist(bootr))/bootnum,byrow=T)))
  
}

#corrected=bootmlr(ncat=3,dataset=iota3,outc="ordinal3O",bootnum=10,appres=results_app,modelformula=ordinal3O ~ oncocenter + Age + lesdmax + propsolid + papnr + Shadows + Ascites,seedval=3436)


# Bootstrap for OLR (cumulative logit with proportional odds)
bootolr <- function(dataset,outc,modelformula,seedval,bootnum,appres,ncat){
  # for some reason, function below does not work without this line, haven't figured out why
  bootr=list()
  set.seed(seedval) # first set the random number seed, then you can redo the bootstrapping and still get the same results
  for (i in 1:bootnum){  
    train_data <- dataset[sample(row.names(dataset),replace=TRUE),]
    # do the modeling
    mb <- vglm(modelformula, family=cumulative(parallel=T), data=train_data)
    
    # predict the values on the bootstrap data
    mpredb <- predictvglm(mb,newdata=train_data,type="response")
    mlpredb <- predictvglm(mb,newdata=train_data,type="link")
    
    # predict the values on the original, unresampled data
    mpredo <- predictvglm(mb,newdata=dataset,type="response")
    mlpredo <- predictvglm(mb,newdata=dataset,type="link")
    
    mcaloutb = calout(out=train_data[,outc],preds=mpredb,k=ncat)
    mcalouto = calout(out=dataset[,outc],preds=mpredo,k=ncat)
    
    mcaldoutb = caldout(out=train_data[,outc],preds=mpredb,k=ncat)
    mcaldouto = caldout(out=dataset[,outc],preds=mpredo,k=ncat)
    
    if(ncat==3){
      mrecalib <- coefficients(vglm(train_data[,outc] ~ 1, offset = mlpredb[,1], family=cumulative(parallel=T)))[1]
      mrecalio <- coefficients(vglm(dataset[,outc] ~ 1, offset = mlpredo[,1], family=cumulative(parallel=T)))[1]
      mrecalsb <- coefficients(vglm(train_data[,outc] ~ mlpredb[,1], family=cumulative(parallel=T)))[c(3)]
      mrecalso <- coefficients(vglm(dataset[,outc] ~ mlpredo[,1], family=cumulative(parallel=T)))[c(3)]
      mrecalib=c(mrecalib, coefficients(vglm(train_data[,outc] ~ 1, offset = mlpredb[,2], family=cumulative(parallel=T)))[2])
      mrecalio=c(mrecalio, coefficients(vglm(dataset[,outc] ~ 1, offset = mlpredo[,2], family=cumulative(parallel=T)))[2])
      mrecalsb=c(mrecalsb, coefficients(vglm(train_data[,outc] ~ mlpredb[,2], family=cumulative(parallel=T)))[c(3)])
      mrecalso=c(mrecalso, coefficients(vglm(dataset[,outc] ~ mlpredo[,2], family=cumulative(parallel=T)))[c(3)])
    }
    if(ncat==4){
      mrecalib <- coefficients(vglm(train_data[,outc] ~ 1, offset = mlpredb[,1], family=cumulative(parallel=T)))[1]
      mrecalio <- coefficients(vglm(dataset[,outc] ~ 1, offset = mlpredo[,1], family=cumulative(parallel=T)))[1]
      mrecalsb <- coefficients(vglm(train_data[,outc] ~ mlpredb[,1], family=cumulative(parallel=T)))[c(4)]
      mrecalso <- coefficients(vglm(dataset[,outc] ~ mlpredo[,1], family=cumulative(parallel=T)))[c(4)]
      mrecalib=c(mrecalib, coefficients(vglm(train_data[,outc] ~ 1, offset = mlpredb[,2], family=cumulative(parallel=T)))[2])
      mrecalio=c(mrecalio, coefficients(vglm(dataset[,outc] ~ 1, offset = mlpredo[,2], family=cumulative(parallel=T)))[2])
      mrecalsb=c(mrecalsb, coefficients(vglm(train_data[,outc] ~ mlpredb[,2], family=cumulative(parallel=T)))[c(4)])
      mrecalso=c(mrecalso, coefficients(vglm(dataset[,outc] ~ mlpredo[,2], family=cumulative(parallel=T)))[c(4)])
      mrecalib=c(mrecalib, coefficients(vglm(train_data[,outc] ~ 1, offset = mlpredb[,3], family=cumulative(parallel=T)))[3])
      mrecalio=c(mrecalio, coefficients(vglm(dataset[,outc] ~ 1, offset = mlpredo[,3], family=cumulative(parallel=T)))[3])
      mrecalsb=c(mrecalsb, coefficients(vglm(train_data[,outc] ~ mlpredb[,3], family=cumulative(parallel=T)))[c(4)])
      mrecalso=c(mrecalso, coefficients(vglm(dataset[,outc] ~ mlpredo[,3], family=cumulative(parallel=T)))[c(4)])
    }
    if(ncat==5){
      mrecalib <- coefficients(vglm(train_data[,outc] ~ 1, offset = mlpredb[,1], family=cumulative(parallel=T)))[1]
      mrecalio <- coefficients(vglm(dataset[,outc] ~ 1, offset = mlpredo[,1], family=cumulative(parallel=T)))[1]
      mrecalsb <- coefficients(vglm(train_data[,outc] ~ mlpredb[,1], family=cumulative(parallel=T)))[c(5)]
      mrecalso <- coefficients(vglm(dataset[,outc] ~ mlpredo[,1], family=cumulative(parallel=T)))[c(5)]
      mrecalib=c(mrecalib, coefficients(vglm(train_data[,outc] ~ 1, offset = mlpredb[,2], family=cumulative(parallel=T)))[2])
      mrecalio=c(mrecalio, coefficients(vglm(dataset[,outc] ~ 1, offset = mlpredo[,2], family=cumulative(parallel=T)))[2])
      mrecalsb=c(mrecalsb, coefficients(vglm(train_data[,outc] ~ mlpredb[,2], family=cumulative(parallel=T)))[c(5)])
      mrecalso=c(mrecalso, coefficients(vglm(dataset[,outc] ~ mlpredo[,2], family=cumulative(parallel=T)))[c(5)])
      mrecalib=c(mrecalib, coefficients(vglm(train_data[,outc] ~ 1, offset = mlpredb[,3], family=cumulative(parallel=T)))[3])
      mrecalio=c(mrecalio, coefficients(vglm(dataset[,outc] ~ 1, offset = mlpredo[,3], family=cumulative(parallel=T)))[3])
      mrecalsb=c(mrecalsb, coefficients(vglm(train_data[,outc] ~ mlpredb[,3], family=cumulative(parallel=T)))[c(5)])
      mrecalso=c(mrecalso, coefficients(vglm(dataset[,outc] ~ mlpredo[,3], family=cumulative(parallel=T)))[c(5)])
      mrecalib=c(mrecalib, coefficients(vglm(train_data[,outc] ~ 1, offset = mlpredb[,4], family=cumulative(parallel=T)))[4])
      mrecalio=c(mrecalio, coefficients(vglm(dataset[,outc] ~ 1, offset = mlpredo[,4], family=cumulative(parallel=T)))[4])
      mrecalsb=c(mrecalsb, coefficients(vglm(train_data[,outc] ~ mlpredb[,4], family=cumulative(parallel=T)))[c(5)])
      mrecalso=c(mrecalso, coefficients(vglm(dataset[,outc] ~ mlpredo[,4], family=cumulative(parallel=T)))[c(5)])
    }

    # Overall recalibration based on E
    mcaleb = calE(out=train_data[,outc],preds=mpredb,k=ncat)
    mcaleo = calE(out=dataset[,outc],preds=mpredo,k=ncat)
    
    # mvgamsmps4bb = vglm(train_data[,outc] ~ rcs(log(mpredb[,1]/(1-mpredb[,1])),3) + rcs(log(mpredb[,2]/(1-mpredb[,2])),3) + rcs(log(mpredb[,3]/(1-mpredb[,3])),3) + rcs(log(mpredb[,4]/(1-mpredb[,4])),3),family=multinomial(refLevel = "1"))
    # mvgamsmps4bo = vglm(dataset[,outc] ~ rcs(log(mpredo[,1]/(1-mpredo[,1])),3) + rcs(log(mpredo[,2]/(1-mpredo[,2])),3) + rcs(log(mpredo[,3]/(1-mpredo[,3])),3) + rcs(log(mpredo[,4]/(1-mpredo[,4])),3),family=multinomial(refLevel = "1"))
    # mECIrbb = eci_rel(calout=mvgamsmps4bb,preds=mpredb,k=5,outc=train_data[,outc])
    # mECIrbo = eci_rel(calout=mvgamsmps4bo,preds=mpredo,k=5,outc=dataset[,outc])
    
    mcb = orc(out=train_data[,outc],preds=mpredb,k=ncat)
    mco = orc(out=dataset[,outc],preds=mpredo,k=ncat)
    
    bootr[[i]] <- c(mcaloutb[,1]-mcalouto[,1],mcaloutb[,2]-mcalouto[,2],mcaldoutb[,1]-mcaldouto[,1],mcaldoutb[,2]-mcaldouto[,2],mrecalib-mrecalio,mrecalsb-mrecalso,mcaleb-mcaleo,0,mcb-mco)#mECIrbb-mECIrbo,mcb-mco)
    
    print(i)
    
  }
  
  cbind(appres,appres-colMeans(matrix(unlist(bootr),nrow=bootnum,ncol=length(unlist(bootr))/bootnum,byrow=T)), colMeans(matrix(unlist(bootr),nrow=bootnum,ncol=length(unlist(bootr))/bootnum,byrow=T)))
  
}


# Bootstrap for adjacent category with PO
bootacp <- function(dataset,outc,modelformula,seedval,bootnum,appres,ncat){
  # for some reason, function below does not work without this line, haven't figured out why
  bootr=list()
  set.seed(seedval) # first set the random number seed, then you can redo the bootstrapping and still get the same results
  for (i in 1:bootnum){  
    train_data <- dataset[sample(row.names(dataset),replace=TRUE),]
    # do the modeling
    mb <- vglm(modelformula, family=acat(parallel=T), data=train_data)
    
    # predict the values on the bootstrap data
    mpredb <- predictvglm(mb,newdata=train_data,type="response")
    mlpredb <- predictvglm(mb,newdata=train_data,type="link")
    
    # predict the values on the original, unresampled data
    mpredo <- predictvglm(mb,newdata=dataset,type="response")
    mlpredo <- predictvglm(mb,newdata=dataset,type="link")
    
    mcaloutb = calout(out=train_data[,outc],preds=mpredb,k=ncat)
    mcalouto = calout(out=dataset[,outc],preds=mpredo,k=ncat)
    
    mcaldoutb = caldout(out=train_data[,outc],preds=mpredb,k=ncat)
    mcaldouto = caldout(out=dataset[,outc],preds=mpredo,k=ncat)
    
    if(ncat==3){
      mrecalib <- coefficients(vglm(train_data[,outc] ~ 1, offset = mlpredb[,1], family=acat(parallel=T)))[1]
      mrecalio <- coefficients(vglm(dataset[,outc] ~ 1, offset = mlpredo[,1], family=acat(parallel=T)))[1]
      mrecalsb <- coefficients(vglm(train_data[,outc] ~ mlpredb[,1], family=acat(parallel=T)))[c(3)]
      mrecalso <- coefficients(vglm(dataset[,outc] ~ mlpredo[,1], family=acat(parallel=T)))[c(3)]
      mrecalib=c(mrecalib, coefficients(vglm(train_data[,outc] ~ 1, offset = mlpredb[,2], family=acat(parallel=T)))[2])
      mrecalio=c(mrecalio, coefficients(vglm(dataset[,outc] ~ 1, offset = mlpredo[,2], family=acat(parallel=T)))[2])
      mrecalsb=c(mrecalsb, coefficients(vglm(train_data[,outc] ~ mlpredb[,2], family=acat(parallel=T)))[c(3)])
      mrecalso=c(mrecalso, coefficients(vglm(dataset[,outc] ~ mlpredo[,2], family=acat(parallel=T)))[c(3)])
    }
    if(ncat==4){
      mrecalib <- coefficients(vglm(train_data[,outc] ~ 1, offset = mlpredb[,1], family=acat(parallel=T)))[1]
      mrecalio <- coefficients(vglm(dataset[,outc] ~ 1, offset = mlpredo[,1], family=acat(parallel=T)))[1]
      mrecalsb <- coefficients(vglm(train_data[,outc] ~ mlpredb[,1], family=acat(parallel=T)))[c(4)]
      mrecalso <- coefficients(vglm(dataset[,outc] ~ mlpredo[,1], family=acat(parallel=T)))[c(4)]
      mrecalib=c(mrecalib, coefficients(vglm(train_data[,outc] ~ 1, offset = mlpredb[,2], family=acat(parallel=T)))[2])
      mrecalio=c(mrecalio, coefficients(vglm(dataset[,outc] ~ 1, offset = mlpredo[,2], family=acat(parallel=T)))[2])
      mrecalsb=c(mrecalsb, coefficients(vglm(train_data[,outc] ~ mlpredb[,2], family=acat(parallel=T)))[c(4)])
      mrecalso=c(mrecalso, coefficients(vglm(dataset[,outc] ~ mlpredo[,2], family=acat(parallel=T)))[c(4)])
      mrecalib=c(mrecalib, coefficients(vglm(train_data[,outc] ~ 1, offset = mlpredb[,3], family=acat(parallel=T)))[3])
      mrecalio=c(mrecalio, coefficients(vglm(dataset[,outc] ~ 1, offset = mlpredo[,3], family=acat(parallel=T)))[3])
      mrecalsb=c(mrecalsb, coefficients(vglm(train_data[,outc] ~ mlpredb[,3], family=acat(parallel=T)))[c(4)])
      mrecalso=c(mrecalso, coefficients(vglm(dataset[,outc] ~ mlpredo[,3], family=acat(parallel=T)))[c(4)])
    }
    if(ncat==5){
      mrecalib <- coefficients(vglm(train_data[,outc] ~ 1, offset = mlpredb[,1], family=acat(parallel=T)))[1]
      mrecalio <- coefficients(vglm(dataset[,outc] ~ 1, offset = mlpredo[,1], family=acat(parallel=T)))[1]
      mrecalsb <- coefficients(vglm(train_data[,outc] ~ mlpredb[,1], family=acat(parallel=T)))[c(5)]
      mrecalso <- coefficients(vglm(dataset[,outc] ~ mlpredo[,1], family=acat(parallel=T)))[c(5)]
      mrecalib=c(mrecalib, coefficients(vglm(train_data[,outc] ~ 1, offset = mlpredb[,2], family=acat(parallel=T)))[2])
      mrecalio=c(mrecalio, coefficients(vglm(dataset[,outc] ~ 1, offset = mlpredo[,2], family=acat(parallel=T)))[2])
      mrecalsb=c(mrecalsb, coefficients(vglm(train_data[,outc] ~ mlpredb[,2], family=acat(parallel=T)))[c(5)])
      mrecalso=c(mrecalso, coefficients(vglm(dataset[,outc] ~ mlpredo[,2], family=acat(parallel=T)))[c(5)])
      mrecalib=c(mrecalib, coefficients(vglm(train_data[,outc] ~ 1, offset = mlpredb[,3], family=acat(parallel=T)))[3])
      mrecalio=c(mrecalio, coefficients(vglm(dataset[,outc] ~ 1, offset = mlpredo[,3], family=acat(parallel=T)))[3])
      mrecalsb=c(mrecalsb, coefficients(vglm(train_data[,outc] ~ mlpredb[,3], family=acat(parallel=T)))[c(5)])
      mrecalso=c(mrecalso, coefficients(vglm(dataset[,outc] ~ mlpredo[,3], family=acat(parallel=T)))[c(5)])
      mrecalib=c(mrecalib, coefficients(vglm(train_data[,outc] ~ 1, offset = mlpredb[,4], family=acat(parallel=T)))[4])
      mrecalio=c(mrecalio, coefficients(vglm(dataset[,outc] ~ 1, offset = mlpredo[,4], family=acat(parallel=T)))[4])
      mrecalsb=c(mrecalsb, coefficients(vglm(train_data[,outc] ~ mlpredb[,4], family=acat(parallel=T)))[c(5)])
      mrecalso=c(mrecalso, coefficients(vglm(dataset[,outc] ~ mlpredo[,4], family=acat(parallel=T)))[c(5)])
    }
    
    # Overall recalibration based on E
    mcaleb = calE(out=train_data[,outc],preds=mpredb,k=ncat)
    mcaleo = calE(out=dataset[,outc],preds=mpredo,k=ncat)
    
    # mvgamsmps4bb = vglm(train_data[,outc] ~ rcs(log(mpredb[,1]/(1-mpredb[,1])),3) + rcs(log(mpredb[,2]/(1-mpredb[,2])),3) + rcs(log(mpredb[,3]/(1-mpredb[,3])),3) + rcs(log(mpredb[,4]/(1-mpredb[,4])),3),family=multinomial(refLevel = "1"))
    # mvgamsmps4bo = vglm(dataset[,outc] ~ rcs(log(mpredo[,1]/(1-mpredo[,1])),3) + rcs(log(mpredo[,2]/(1-mpredo[,2])),3) + rcs(log(mpredo[,3]/(1-mpredo[,3])),3) + rcs(log(mpredo[,4]/(1-mpredo[,4])),3),family=multinomial(refLevel = "1"))
    # mECIrbb = eci_rel(calout=mvgamsmps4bb,preds=mpredb,k=5,outc=train_data[,outc])
    # mECIrbo = eci_rel(calout=mvgamsmps4bo,preds=mpredo,k=5,outc=dataset[,outc])
    
    mcb = orc(out=train_data[,outc],preds=mpredb,k=ncat)
    mco = orc(out=dataset[,outc],preds=mpredo,k=ncat)
    
    bootr[[i]] <- c(mcaloutb[,1]-mcalouto[,1],mcaloutb[,2]-mcalouto[,2],mcaldoutb[,1]-mcaldouto[,1],mcaldoutb[,2]-mcaldouto[,2],mrecalib-mrecalio,mrecalsb-mrecalso,mcaleb-mcaleo,0,mcb-mco)#mECIrbb-mECIrbo,mcb-mco)
    
    print(i)
    
  }
  
  cbind(appres,appres-colMeans(matrix(unlist(bootr),nrow=bootnum,ncol=length(unlist(bootr))/bootnum,byrow=T)), colMeans(matrix(unlist(bootr),nrow=bootnum,ncol=length(unlist(bootr))/bootnum,byrow=T)))
  
}


# Bootstrap for continuation ratio with proportional odds
bootcrp <- function(dataset,outc,modelformula,seedval,bootnum,appres,ncat){
  # for some reason, function below does not work without this line, haven't figured out why
  bootr=list()
  set.seed(seedval) # first set the random number seed, then you can redo the bootstrapping and still get the same results
  for (i in 1:bootnum){  
    train_data <- dataset[sample(row.names(dataset),replace=TRUE),]
    # do the modeling
    mb <- vglm(modelformula, family=cratio(parallel=T), data=train_data)
    
    # predict the values on the bootstrap data
    mpredb <- predictvglm(mb,newdata=train_data,type="response")
    mlpredb <- predictvglm(mb,newdata=train_data,type="link")
    
    # predict the values on the original, unresampled data
    mpredo <- predictvglm(mb,newdata=dataset,type="response")
    mlpredo <- predictvglm(mb,newdata=dataset,type="link")
    
    mcaloutb = calout(out=train_data[,outc],preds=mpredb,k=ncat)
    mcalouto = calout(out=dataset[,outc],preds=mpredo,k=ncat)
    
    mcaldoutb = caldout(out=train_data[,outc],preds=mpredb,k=ncat)
    mcaldouto = caldout(out=dataset[,outc],preds=mpredo,k=ncat)
    
    if(ncat==3){
      mrecalib <- coefficients(vglm(train_data[,outc] ~ 1, offset = mlpredb[,1], family=cratio(parallel=T)))[1]
      mrecalio <- coefficients(vglm(dataset[,outc] ~ 1, offset = mlpredo[,1], family=cratio(parallel=T)))[1]
      mrecalsb <- coefficients(vglm(train_data[,outc] ~ mlpredb[,1], family=cratio(parallel=T)))[c(3)]
      mrecalso <- coefficients(vglm(dataset[,outc] ~ mlpredo[,1], family=cratio(parallel=T)))[c(3)]
      mrecalib=c(mrecalib, coefficients(vglm(train_data[,outc] ~ 1, offset = mlpredb[,2], family=cratio(parallel=T)))[2])
      mrecalio=c(mrecalio, coefficients(vglm(dataset[,outc] ~ 1, offset = mlpredo[,2], family=cratio(parallel=T)))[2])
      mrecalsb=c(mrecalsb, coefficients(vglm(train_data[,outc] ~ mlpredb[,2], family=cratio(parallel=T)))[c(3)])
      mrecalso=c(mrecalso, coefficients(vglm(dataset[,outc] ~ mlpredo[,2], family=cratio(parallel=T)))[c(3)])
    }
    if(ncat==4){
      mrecalib <- coefficients(vglm(train_data[,outc] ~ 1, offset = mlpredb[,1], family=cratio(parallel=T)))[1]
      mrecalio <- coefficients(vglm(dataset[,outc] ~ 1, offset = mlpredo[,1], family=cratio(parallel=T)))[1]
      mrecalsb <- coefficients(vglm(train_data[,outc] ~ mlpredb[,1], family=cratio(parallel=T)))[c(4)]
      mrecalso <- coefficients(vglm(dataset[,outc] ~ mlpredo[,1], family=cratio(parallel=T)))[c(4)]
      mrecalib=c(mrecalib, coefficients(vglm(train_data[,outc] ~ 1, offset = mlpredb[,2], family=cratio(parallel=T)))[2])
      mrecalio=c(mrecalio, coefficients(vglm(dataset[,outc] ~ 1, offset = mlpredo[,2], family=cratio(parallel=T)))[2])
      mrecalsb=c(mrecalsb, coefficients(vglm(train_data[,outc] ~ mlpredb[,2], family=cratio(parallel=T)))[c(4)])
      mrecalso=c(mrecalso, coefficients(vglm(dataset[,outc] ~ mlpredo[,2], family=cratio(parallel=T)))[c(4)])
      mrecalib=c(mrecalib, coefficients(vglm(train_data[,outc] ~ 1, offset = mlpredb[,3], family=cratio(parallel=T)))[3])
      mrecalio=c(mrecalio, coefficients(vglm(dataset[,outc] ~ 1, offset = mlpredo[,3], family=cratio(parallel=T)))[3])
      mrecalsb=c(mrecalsb, coefficients(vglm(train_data[,outc] ~ mlpredb[,3], family=cratio(parallel=T)))[c(4)])
      mrecalso=c(mrecalso, coefficients(vglm(dataset[,outc] ~ mlpredo[,3], family=cratio(parallel=T)))[c(4)])
    }
    if(ncat==5){
      mrecalib <- coefficients(vglm(train_data[,outc] ~ 1, offset = mlpredb[,1], family=cratio(parallel=T)))[1]
      mrecalio <- coefficients(vglm(dataset[,outc] ~ 1, offset = mlpredo[,1], family=cratio(parallel=T)))[1]
      mrecalsb <- coefficients(vglm(train_data[,outc] ~ mlpredb[,1], family=cratio(parallel=T)))[c(5)]
      mrecalso <- coefficients(vglm(dataset[,outc] ~ mlpredo[,1], family=cratio(parallel=T)))[c(5)]
      mrecalib=c(mrecalib, coefficients(vglm(train_data[,outc] ~ 1, offset = mlpredb[,2], family=cratio(parallel=T)))[2])
      mrecalio=c(mrecalio, coefficients(vglm(dataset[,outc] ~ 1, offset = mlpredo[,2], family=cratio(parallel=T)))[2])
      mrecalsb=c(mrecalsb, coefficients(vglm(train_data[,outc] ~ mlpredb[,2], family=cratio(parallel=T)))[c(5)])
      mrecalso=c(mrecalso, coefficients(vglm(dataset[,outc] ~ mlpredo[,2], family=cratio(parallel=T)))[c(5)])
      mrecalib=c(mrecalib, coefficients(vglm(train_data[,outc] ~ 1, offset = mlpredb[,3], family=cratio(parallel=T)))[3])
      mrecalio=c(mrecalio, coefficients(vglm(dataset[,outc] ~ 1, offset = mlpredo[,3], family=cratio(parallel=T)))[3])
      mrecalsb=c(mrecalsb, coefficients(vglm(train_data[,outc] ~ mlpredb[,3], family=cratio(parallel=T)))[c(5)])
      mrecalso=c(mrecalso, coefficients(vglm(dataset[,outc] ~ mlpredo[,3], family=cratio(parallel=T)))[c(5)])
      mrecalib=c(mrecalib, coefficients(vglm(train_data[,outc] ~ 1, offset = mlpredb[,4], family=cratio(parallel=T)))[4])
      mrecalio=c(mrecalio, coefficients(vglm(dataset[,outc] ~ 1, offset = mlpredo[,4], family=cratio(parallel=T)))[4])
      mrecalsb=c(mrecalsb, coefficients(vglm(train_data[,outc] ~ mlpredb[,4], family=cratio(parallel=T)))[c(5)])
      mrecalso=c(mrecalso, coefficients(vglm(dataset[,outc] ~ mlpredo[,4], family=cratio(parallel=T)))[c(5)])
    }
    
    # Overall recalibration based on E
    mcaleb = calE(out=train_data[,outc],preds=mpredb,k=ncat)
    mcaleo = calE(out=dataset[,outc],preds=mpredo,k=ncat)
    
    # mvgamsmps4bb = vglm(train_data[,outc] ~ rcs(log(mpredb[,1]/(1-mpredb[,1])),3) + rcs(log(mpredb[,2]/(1-mpredb[,2])),3) + rcs(log(mpredb[,3]/(1-mpredb[,3])),3) + rcs(log(mpredb[,4]/(1-mpredb[,4])),3),family=multinomial(refLevel = "1"))
    # mvgamsmps4bo = vglm(dataset[,outc] ~ rcs(log(mpredo[,1]/(1-mpredo[,1])),3) + rcs(log(mpredo[,2]/(1-mpredo[,2])),3) + rcs(log(mpredo[,3]/(1-mpredo[,3])),3) + rcs(log(mpredo[,4]/(1-mpredo[,4])),3),family=multinomial(refLevel = "1"))
    # mECIrbb = eci_rel(calout=mvgamsmps4bb,preds=mpredb,k=5,outc=train_data[,outc])
    # mECIrbo = eci_rel(calout=mvgamsmps4bo,preds=mpredo,k=5,outc=dataset[,outc])
    
    mcb = orc(out=train_data[,outc],preds=mpredb,k=ncat)
    mco = orc(out=dataset[,outc],preds=mpredo,k=ncat)
    
    bootr[[i]] <- c(mcaloutb[,1]-mcalouto[,1],mcaloutb[,2]-mcalouto[,2],mcaldoutb[,1]-mcaldouto[,1],mcaldoutb[,2]-mcaldouto[,2],mrecalib-mrecalio,mrecalsb-mrecalso,mcaleb-mcaleo,0,mcb-mco)#mECIrbb-mECIrbo,mcb-mco)
    
    print(i)
    
  }
  
  cbind(appres,appres-colMeans(matrix(unlist(bootr),nrow=bootnum,ncol=length(unlist(bootr))/bootnum,byrow=T)), colMeans(matrix(unlist(bootr),nrow=bootnum,ncol=length(unlist(bootr))/bootnum,byrow=T)))
  
}


# Bootstrap for continuation ratio without proportional odds
bootcrnp <- function(dataset,outc,modelformula,seedval,bootnum,appres,ncat){
  # for some reason, function below does not work without this line, haven't figured out why
  bootr=list()
  set.seed(seedval) # first set the random number seed, then you can redo the bootstrapping and still get the same results
  for (i in 1:bootnum){  
    train_data <- dataset[sample(row.names(dataset),replace=TRUE),]
    # do the modeling
    mb <- vglm(modelformula, family=cratio(parallel = F), data=train_data)
    
    # predict the values on the bootstrap data
    mpredb <- predictvglm(mb,newdata=train_data,type="response")
    mlpredb <- predictvglm(mb,newdata=train_data,type="link")
    
    # predict the values on the original, unresampled data
    mpredo <- predictvglm(mb,newdata=dataset,type="response")
    mlpredo <- predictvglm(mb,newdata=dataset,type="link")
    
    mcaloutb = calout(out=train_data[,outc],preds=mpredb,k=ncat)
    mcalouto = calout(out=dataset[,outc],preds=mpredo,k=ncat)
    
    mcaldoutb = caldout(out=train_data[,outc],preds=mpredb,k=ncat)
    mcaldouto = caldout(out=dataset[,outc],preds=mpredo,k=ncat)
    
    if(ncat==3){
      mrecalib <- coefficients(vglm(train_data[,outc] ~ 1, offset = mlpredb[,1:2], family=cratio(parallel=F)))[c(1,2)]
      mrecalio <- coefficients(vglm(dataset[,outc] ~ 1, offset = mlpredo[,1:2], family=cratio(parallel=F)))[c(1,2)]
      mrecalsb <- coefficients(vglm(train_data[,outc] ~ mlpredb[,1] + mlpredb[,2], constraints=list("(Intercept)"=diag(2),"mlpredb[, 1]"=rbind(1,0),"mlpredb[, 2]"=rbind(0,1)), family=cratio(parallel=F)))[c(3,4)]
      mrecalso <- coefficients(vglm(dataset[,outc] ~ mlpredo[,1] + mlpredo[,2], constraints=list("(Intercept)"=diag(2),"mlpredo[, 1]"=rbind(1,0),"mlpredo[, 2]"=rbind(0,1)), family=cratio(parallel=F)))[c(3,4)]
    } else {
      if(ncat==4) {
        mrecalib <- coefficients(vglm(train_data[,outc] ~ 1, offset = mlpredb[,1:3], family=cratio(parallel=F)))[c(1,2,3)]
        mrecalio <- coefficients(vglm(dataset[,outc] ~ 1, offset = mlpredo[,1:3], family=cratio(parallel=F)))[c(1,2,3)]
        mrecalsb <- coefficients(vglm(train_data[,outc] ~ mlpredb[,1] + mlpredb[,2] + mlpredb[,3], constraints=list("(Intercept)"=diag(3),"mlpredb[, 1]"=rbind(1,0,0),"mlpredb[, 2]"=rbind(0,1,0),"mlpredb[, 3]"=rbind(0,0,1)), family=cratio(parallel=F)))[c(4:6)]
        mrecalso <- coefficients(vglm(dataset[,outc] ~ mlpredo[,1] + mlpredo[,2] + mlpredo[,3], constraints=list("(Intercept)"=diag(3),"mlpredo[, 1]"=rbind(1,0,0),"mlpredo[, 2]"=rbind(0,1,0),"mlpredo[, 3]"=rbind(0,0,1)), family=cratio(parallel=F)))[c(4:6)]
      } else {
        if(ncat==5) {
          mrecalib <- coefficients(vglm(train_data[,outc] ~ 1, offset = mlpredb[,1:4], family=cratio(parallel=F)))[c(1:4)]
          mrecalio <- coefficients(vglm(dataset[,outc] ~ 1, offset = mlpredo[,1:4], family=cratio(parallel=F)))[c(1:4)]
          mrecalsb <- coefficients(vglm(train_data[,outc] ~ mlpredb[,1] + mlpredb[,2] + mlpredb[,3] + mlpredb[,4], constraints=list("(Intercept)"=diag(4),"mlpredb[, 1]"=rbind(1,0,0,0),"mlpredb[, 2]"=rbind(0,1,0,0),"mlpredb[, 3]"=rbind(0,0,1,0),"mlpredb[, 4]"=rbind(0,0,0,1)), family=cratio(parallel=F)))[c(5:8)]
          mrecalso <- coefficients(vglm(dataset[,outc] ~ mlpredo[,1] + mlpredo[,2] + mlpredo[,3] + mlpredo[,4], constraints=list("(Intercept)"=diag(4),"mlpredo[, 1]"=rbind(1,0,0,0),"mlpredo[, 2]"=rbind(0,1,0,0),"mlpredo[, 3]"=rbind(0,0,1,0),"mlpredo[, 4]"=rbind(0,0,0,1)), family=cratio(parallel=F)))[c(5:8)]
        }
      }
    }
    
    # Overall recalibration based on E
    mcaleb = calE(out=train_data[,outc],preds=mpredb,k=ncat)
    mcaleo = calE(out=dataset[,outc],preds=mpredo,k=ncat)
    
    # mvgamsmps4bb = vglm(train_data[,outc] ~ rcs(log(mpredb[,1]/(1-mpredb[,1])),3) + rcs(log(mpredb[,2]/(1-mpredb[,2])),3) + rcs(log(mpredb[,3]/(1-mpredb[,3])),3) + rcs(log(mpredb[,4]/(1-mpredb[,4])),3),family=multinomial(refLevel = "1"))
    # mvgamsmps4bo = vglm(dataset[,outc] ~ rcs(log(mpredo[,1]/(1-mpredo[,1])),3) + rcs(log(mpredo[,2]/(1-mpredo[,2])),3) + rcs(log(mpredo[,3]/(1-mpredo[,3])),3) + rcs(log(mpredo[,4]/(1-mpredo[,4])),3),family=multinomial(refLevel = "1"))
    # mECIrbb = eci_rel(calout=mvgamsmps4bb,preds=mpredb,k=5,outc=train_data[,outc])
    # mECIrbo = eci_rel(calout=mvgamsmps4bo,preds=mpredo,k=5,outc=dataset[,outc])
    
    mcb = orc(out=train_data[,outc],preds=mpredb,k=ncat)
    mco = orc(out=dataset[,outc],preds=mpredo,k=ncat)
    
    bootr[[i]] <- c(mcaloutb[,1]-mcalouto[,1],mcaloutb[,2]-mcalouto[,2],mcaldoutb[,1]-mcaldouto[,1],mcaldoutb[,2]-mcaldouto[,2],mrecalib-mrecalio,mrecalsb-mrecalso,mcaleb-mcaleo,0,mcb-mco)#mECIrbb-mECIrbo,mcb-mco)
    
    print(i)
    
  }
  
  cbind(appres,appres-colMeans(matrix(unlist(bootr),nrow=bootnum,ncol=length(unlist(bootr))/bootnum,byrow=T)), colMeans(matrix(unlist(bootr),nrow=bootnum,ncol=length(unlist(bootr))/bootnum,byrow=T)))
  
}


# Bootstrap for stereotype model
bootsm <- function(dataset,outc,modelformula,seedval,bootnum,appres,ncat){
  # for some reason, function below does not work without this line, haven't figured out why
  bootr=list()
  set.seed(seedval) # first set the random number seed, then you can redo the bootstrapping and still get the same results
  for (i in 1:bootnum){  
    train_data <- dataset[sample(row.names(dataset),replace=TRUE),]
    # do the modeling
    mb <- rrvglm(modelformula, family=multinomial(refLevel = "1"), data=train_data)
    
    # predict the values on the bootstrap data
    mpredb <- predictvglm(mb,newdata=train_data,type="response")
    mlpredb <- predictvglm(mb,newdata=train_data,type="link")
    
    # predict the values on the original, unresampled data
    mpredo <- predictvglm(mb,newdata=dataset,type="response")
    mlpredo <- predictvglm(mb,newdata=dataset,type="link")
    
    mcaloutb = calout(out=train_data[,outc],preds=mpredb,k=ncat)
    mcalouto = calout(out=dataset[,outc],preds=mpredo,k=ncat)
    
    mcaldoutb = caldout(out=train_data[,outc],preds=mpredb,k=ncat)
    mcaldouto = caldout(out=dataset[,outc],preds=mpredo,k=ncat)
    
    if(ncat==3){
      mrecalib <- coefficients(vglm(train_data[,outc] ~ 1, offset = mlpredb[,1:2], family=multinomial(refLevel = "1")))[c(1,2)]
      mrecalio <- coefficients(vglm(dataset[,outc] ~ 1, offset = mlpredo[,1:2], family=multinomial(refLevel = "1")))[c(1,2)]
      mrecalsb <- coefficients(vglm(train_data[,outc] ~ mlpredb[,1] + mlpredb[,2], constraints=list("(Intercept)"=diag(2),"mlpredb[, 1]"=rbind(1,0),"mlpredb[, 2]"=rbind(0,1)), family=multinomial(refLevel = "1")))[c(3,4)]
      mrecalso <- coefficients(vglm(dataset[,outc] ~ mlpredo[,1] + mlpredo[,2], constraints=list("(Intercept)"=diag(2),"mlpredo[, 1]"=rbind(1,0),"mlpredo[, 2]"=rbind(0,1)), family=multinomial(refLevel = "1")))[c(3,4)]
    } else {
      if(ncat==4) {
        mrecalib <- coefficients(vglm(train_data[,outc] ~ 1, offset = mlpredb[,1:3], family=multinomial(refLevel = "1")))[c(1,2,3)]
        mrecalio <- coefficients(vglm(dataset[,outc] ~ 1, offset = mlpredo[,1:3], family=multinomial(refLevel = "1")))[c(1,2,3)]
        mrecalsb <- coefficients(vglm(train_data[,outc] ~ mlpredb[,1] + mlpredb[,2] + mlpredb[,3], constraints=list("(Intercept)"=diag(3),"mlpredb[, 1]"=rbind(1,0,0),"mlpredb[, 2]"=rbind(0,1,0),"mlpredb[, 3]"=rbind(0,0,1)), family=multinomial(refLevel = "1")))[c(4:6)]
        mrecalso <- coefficients(vglm(dataset[,outc] ~ mlpredo[,1] + mlpredo[,2] + mlpredo[,3], constraints=list("(Intercept)"=diag(3),"mlpredo[, 1]"=rbind(1,0,0),"mlpredo[, 2]"=rbind(0,1,0),"mlpredo[, 3]"=rbind(0,0,1)), family=multinomial(refLevel = "1")))[c(4:6)]
      } else {
        if(ncat==5) {
          mrecalib <- coefficients(vglm(train_data[,outc] ~ 1, offset = mlpredb[,1:4], family=multinomial(refLevel = "1")))[c(1:4)]
          mrecalio <- coefficients(vglm(dataset[,outc] ~ 1, offset = mlpredo[,1:4], family=multinomial(refLevel = "1")))[c(1:4)]
          mrecalsb <- coefficients(vglm(train_data[,outc] ~ mlpredb[,1] + mlpredb[,2] + mlpredb[,3] + mlpredb[,4], constraints=list("(Intercept)"=diag(4),"mlpredb[, 1]"=rbind(1,0,0,0),"mlpredb[, 2]"=rbind(0,1,0,0),"mlpredb[, 3]"=rbind(0,0,1,0),"mlpredb[, 4]"=rbind(0,0,0,1)), family=multinomial(refLevel = "1")))[c(5:8)]
          mrecalso <- coefficients(vglm(dataset[,outc] ~ mlpredo[,1] + mlpredo[,2] + mlpredo[,3] + mlpredo[,4], constraints=list("(Intercept)"=diag(4),"mlpredo[, 1]"=rbind(1,0,0,0),"mlpredo[, 2]"=rbind(0,1,0,0),"mlpredo[, 3]"=rbind(0,0,1,0),"mlpredo[, 4]"=rbind(0,0,0,1)), family=multinomial(refLevel = "1")))[c(5:8)]
        }
      }
    }
    
    # Overall recalibration based on E
    mcaleb = calE(out=train_data[,outc],preds=mpredb,k=ncat)
    mcaleo = calE(out=dataset[,outc],preds=mpredo,k=ncat)
    
    # mvgamsmps4bb = vglm(train_data[,outc] ~ rcs(log(mpredb[,1]/(1-mpredb[,1])),3) + rcs(log(mpredb[,2]/(1-mpredb[,2])),3) + rcs(log(mpredb[,3]/(1-mpredb[,3])),3) + rcs(log(mpredb[,4]/(1-mpredb[,4])),3),family=multinomial(refLevel = "1"))
    # mvgamsmps4bo = vglm(dataset[,outc] ~ rcs(log(mpredo[,1]/(1-mpredo[,1])),3) + rcs(log(mpredo[,2]/(1-mpredo[,2])),3) + rcs(log(mpredo[,3]/(1-mpredo[,3])),3) + rcs(log(mpredo[,4]/(1-mpredo[,4])),3),family=multinomial(refLevel = "1"))
    # mECIrbb = eci_rel(calout=mvgamsmps4bb,preds=mpredb,k=5,outc=train_data[,outc])
    # mECIrbo = eci_rel(calout=mvgamsmps4bo,preds=mpredo,k=5,outc=dataset[,outc])
    
    mcb = orc(out=train_data[,outc],preds=mpredb,k=ncat)
    mco = orc(out=dataset[,outc],preds=mpredo,k=ncat)
    
    bootr[[i]] <- c(mcaloutb[,1]-mcalouto[,1],mcaloutb[,2]-mcalouto[,2],mcaldoutb[,1]-mcaldouto[,1],mcaldoutb[,2]-mcaldouto[,2],mrecalib-mrecalio,mrecalsb-mrecalso,mcaleb-mcaleo,0,mcb-mco)#mECIrbb-mECIrbo,mcb-mco)
    
    print(i)
    
  }
  
  cbind(appres,appres-colMeans(matrix(unlist(bootr),nrow=bootnum,ncol=length(unlist(bootr))/bootnum,byrow=T)), colMeans(matrix(unlist(bootr),nrow=bootnum,ncol=length(unlist(bootr))/bootnum,byrow=T)))
  
}
