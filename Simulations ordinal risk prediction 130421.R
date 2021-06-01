####################################################
####################################################
# ORDINAL RISK PREDICTION MODELS: SIMULATION STUDY #
####################################################
####################################################

# Disclaimer: written by Ben Van Calster. The efficiency/speed of this code was "not necessarily optimized",
# and may therefore not qualify for the R coding Champions League. Apologies to those who may be disappointed.
# But it did what it was supposed to do.


#############
# LIBRARIES #
#############

install.packages("pacman")
pacman::p_load(nnet,gmodels,psych,MASS,VGAM,Hmisc,pROC,clusterPower,ordinal, DescTools,R.utils,mvtnorm,robustbase,rms)

# Set working directory to save files
setwd("C:\\Ben\\Methodological research\\Ordinal prediction models")


#############
# FUNCTIONS #
#############

# Load ordinal functions; adjusted path as necessary!
source("ordcalfunctions.R")
source("ordcalfunctions_simulationstudy.R")


###############################
# SIMULATIONS UNDER MLR TRUTH #
###############################

# SCENARIO 1

{
# Get large dataset
set.seed(34634) # the seed for the large datasets
datapo=getdata_k3p4(mu_x1=c(0,0.4,0.8),
                    mu_x2=c(0,0.3,0.6),
                    mu_x3=c(0,0.4,0.8),
                    mu_x4=c(0,0.3,0.6),
                    yprev=c(33.34,33.33,33.33),#yprev=c(55,30,15),
                    mult=2000)

# Run models for large dataset
mlrfit=vglm(y ~ x1 + x2 + x3 + x4, family=multinomial(refLevel = "1"), data=datapo[[1]])
acnpres=mlrfunc(datafit=datapo[[1]],dataval=datapo[[1]],ECI="nonpar",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
olrres=olrfunc(datafit=datapo[[1]],dataval=datapo[[1]],ECI="nonpar",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
acpres=acpfunc(datafit=datapo[[1]],dataval=datapo[[1]],ECI="nonpar",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
smres=smfunc(datafit=datapo[[1]],dataval=datapo[[1]],ECI="nonpar",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
comsave(filep1="Ordinal prediction model simulations\\ordsimk3p4_equi_equal_moderateorc.RData")
vgamsmps4a = vgam(datapo[[1]]$y ~ sm.ps(log(acnpres[[3]][,2]/acnpres[[3]][,1]),df=4) + sm.ps(log(acnpres[[3]][,3]/acnpres[[3]][,1]),df=4),family=multinomial(refLevel = 1))
vgamsmps4b = vgam(datapo[[1]]$y ~ sm.ps(log(olrres[[3]][,2]/olrres[[3]][,1]),df=4) + sm.ps(log(olrres[[3]][,3]/olrres[[3]][,1]),df=4),family=multinomial(refLevel = 1))
vgamsmps4c = vgam(datapo[[1]]$y ~ sm.ps(log(acpres[[3]][,2]/acpres[[3]][,1]),df=4) + sm.ps(log(acpres[[3]][,3]/acpres[[3]][,1]),df=4),family=multinomial(refLevel = 1))
vgamsmps4d = vgam(datapo[[1]]$y ~ sm.ps(log(smres[[3]][,2]/smres[[3]][,1]),df=4) + sm.ps(log(smres[[3]][,3]/smres[[3]][,1]),df=4),family=multinomial(refLevel = 1))
plotsave(filep2="Ordinal prediction model simulations\\ordsimk3p4_equi_equal_moderateorc.pdf")

# Smooth calibration plots for paper
par(mfrow=c(2,2))
plotlines3sub(preds=acnpres[[3]],obs=vgamsmps4a,mtxt=NULL,subsample0=1000,cexlab=1.4,showleg=0,showtext=1,texttoshow="MLR",cextext=1.3)
plotlines3sub(preds=olrres[[3]],obs=vgamsmps4b,mtxt=NULL,subsample0=1000,cexlab=1.4,showleg=0,showtext=1,texttoshow="CL-PO",cextext=1.3)
plotlines3sub(preds=acpres[[3]],obs=vgamsmps4c,mtxt=NULL,subsample0=1000,cexlab=1.4,showleg=0,showtext=1,texttoshow="AC-PO",cextext=1.3)
plotlines3sub(preds=smres[[3]],obs=vgamsmps4d,mtxt=NULL,subsample0=1000,cexlab=1.4,showleg=0,showtext=1,texttoshow="SLM",cextext=1.3)
# Scatter plot vs true risk for paper
plotvstruerisk3(subsample=1000)

# dev.off()
# olrx <- orm(y ~ x1 + x2 + x3 +x4, data=datapo[[1]], x=T)
# par(mfrow=c(2,2))
# resid(olrx, 'partial',pl=T,label.curves=F) 

# small sample simulations to check overfitting (N=100)
set.seed(4682)
mlrsim=matrix(data=NA,nrow=200,ncol=27)
olrsim=matrix(data=NA,nrow=200,ncol=27)
acpsim=matrix(data=NA,nrow=200,ncol=27)
smsim=matrix(data=NA,nrow=200,ncol=27)
for(i in 1:200) {
  datatr=getdata_k3p4(mu_x1=c(0,0.4,0.8),
                      mu_x2=c(0,0.3,0.6),
                      mu_x3=c(0,0.4,0.8),
                      mu_x4=c(0,0.3,0.6),
                      yprev=c(33.34,33.33,33.33),
                      mult=1)
  mlrrun=mlrfunc(datafit=datatr[[1]],dataval=datapo[[1]],ECI="rcs3",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
  mlrsim[i,1:length(mlrrun[[1]])]=mlrrun[[1]]
  olrrun=olrfunc(datafit=datatr[[1]],dataval=datapo[[1]],ECI="rcs3",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
  olrsim[i,1:length(olrrun[[1]])]=olrrun[[1]]
  acprun=acpfunc(datafit=datatr[[1]],dataval=datapo[[1]],ECI="rcs3",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
  acpsim[i,1:length(acprun[[1]])]=acprun[[1]]
  smrun=smfunc(datafit=datatr[[1]],dataval=datapo[[1]],ECI="rcs3",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
  smsim[i,1:length(smrun[[1]])]=smrun[[1]]
  print(i)
}
comsave_smalln(filep3="Ordinal prediction model simulations\\ordsimk3p4_equi_equal_moderateorc_n100.RData")

# small sample simulations to check overfitting (N=500)
set.seed(4682)
mlrsim=matrix(data=NA,nrow=200,ncol=27)
olrsim=matrix(data=NA,nrow=200,ncol=27)
acpsim=matrix(data=NA,nrow=200,ncol=27)
smsim=matrix(data=NA,nrow=200,ncol=27)
for(i in 1:200) {
  datatr=getdata_k3p4(mu_x1=c(0,0.4,0.8),
                      mu_x2=c(0,0.3,0.6),
                      mu_x3=c(0,0.4,0.8),
                      mu_x4=c(0,0.3,0.6),
                      yprev=c(33.34,33.33,33.33),
                      mult=5)
  mlrrun=mlrfunc(datafit=datatr[[1]],dataval=datapo[[1]],ECI="rcs3",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
  mlrsim[i,1:length(mlrrun[[1]])]=mlrrun[[1]]
  olrrun=olrfunc(datafit=datatr[[1]],dataval=datapo[[1]],ECI="rcs3",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
  olrsim[i,1:length(olrrun[[1]])]=olrrun[[1]]
  acprun=acpfunc(datafit=datatr[[1]],dataval=datapo[[1]],ECI="rcs3",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
  acpsim[i,1:length(acprun[[1]])]=acprun[[1]]
  smrun=smfunc(datafit=datatr[[1]],dataval=datapo[[1]],ECI="rcs3",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
  smsim[i,1:length(smrun[[1]])]=smrun[[1]]
  print(i)
}
comsave_smalln(filep3="Ordinal prediction model simulations\\ordsimk3p4_equi_equal_moderateorc_n500.RData")

}


# SCENARIO 2

{
# Get large dataset
set.seed(34634) # the seed for the large datasets
datapo=getdata_k3p4(mu_x1=c(0,0.4,0.8),
                    mu_x2=c(0,0.3,0.6),
                    mu_x3=c(0,0.4,0.8),
                    mu_x4=c(0,0.3,0.6),
                    yprev=c(55,30,15),#yprev=c(55,30,15),
                    mult=2000)

# Run models for large dataset
mlrfit=vglm(y ~ x1 + x2 + x3 + x4, family=multinomial(refLevel = "1"), data=datapo[[1]])
acnpres=mlrfunc(datafit=datapo[[1]],dataval=datapo[[1]],ECI="nonpar",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
olrres=olrfunc(datafit=datapo[[1]],dataval=datapo[[1]],ECI="nonpar",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
acpres=acpfunc(datafit=datapo[[1]],dataval=datapo[[1]],ECI="nonpar",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
smres=smfunc(datafit=datapo[[1]],dataval=datapo[[1]],ECI="nonpar",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
comsave(filep="Ordinal prediction model simulations\\ordsimk3p4_equi_unequal_moderateorc.RData")
vgamsmps4a = vgam(datapo[[1]]$y ~ sm.ps(log(acnpres[[3]][,2]/acnpres[[3]][,1]),df=4) + sm.ps(log(acnpres[[3]][,3]/acnpres[[3]][,1]),df=4),family=multinomial(refLevel = 1))
vgamsmps4b = vgam(datapo[[1]]$y ~ sm.ps(log(olrres[[3]][,2]/olrres[[3]][,1]),df=4) + sm.ps(log(olrres[[3]][,3]/olrres[[3]][,1]),df=4),family=multinomial(refLevel = 1))
vgamsmps4c = vgam(datapo[[1]]$y ~ sm.ps(log(acpres[[3]][,2]/acpres[[3]][,1]),df=4) + sm.ps(log(acpres[[3]][,3]/acpres[[3]][,1]),df=4),family=multinomial(refLevel = 1))
vgamsmps4d = vgam(datapo[[1]]$y ~ sm.ps(log(smres[[3]][,2]/smres[[3]][,1]),df=4) + sm.ps(log(smres[[3]][,3]/smres[[3]][,1]),df=4),family=multinomial(refLevel = 1))
plotsave(filep2="Ordinal prediction model simulations\\ordsimk3p4_equi_unequal_moderateorc.pdf")

# Smooth calibration plots for paper
par(mfrow=c(2,2))
plotlines3sub(preds=acnpres[[3]],obs=vgamsmps4a,mtxt=NULL,subsample0=1000,cexlab=1.4,showleg=0,showtext=1,texttoshow="MLR",cextext=1.3)
plotlines3sub(preds=olrres[[3]],obs=vgamsmps4b,mtxt=NULL,subsample0=1000,cexlab=1.4,showleg=0,showtext=1,texttoshow="CL-PO",cextext=1.3)
plotlines3sub(preds=acpres[[3]],obs=vgamsmps4c,mtxt=NULL,subsample0=1000,cexlab=1.4,showleg=0,showtext=1,texttoshow="AC-PO",cextext=1.3)
plotlines3sub(preds=smres[[3]],obs=vgamsmps4d,mtxt=NULL,subsample0=1000,cexlab=1.4,showleg=0,showtext=1,texttoshow="SLM",cextext=1.3)
# Scatter plot vs true risk for paper
plotvstruerisk3(subsample=1000)

# dev.off()
# olrx <- orm(y ~ x1 + x2 + x3 +x4, data=datapo[[1]], x=T)
# par(mfrow=c(2,2))
# resid(olrx, 'partial',pl=T,label.curves=F) 

# small sample simulations to check overfitting (N=100)
set.seed(4682)
mlrsim=matrix(data=NA,nrow=200,ncol=27)
olrsim=matrix(data=NA,nrow=200,ncol=27)
acpsim=matrix(data=NA,nrow=200,ncol=27)
smsim=matrix(data=NA,nrow=200,ncol=27)
for(i in 1:200) {
  datatr=getdata_k3p4(mu_x1=c(0,0.4,0.8),
                      mu_x2=c(0,0.3,0.6),
                      mu_x3=c(0,0.4,0.8),
                      mu_x4=c(0,0.3,0.6),
                      yprev=c(55,30,15),
                      mult=1)
  mlrrun=mlrfunc(datafit=datatr[[1]],dataval=datapo[[1]],ECI="rcs3",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
  mlrsim[i,1:length(mlrrun[[1]])]=mlrrun[[1]]
  olrrun=olrfunc(datafit=datatr[[1]],dataval=datapo[[1]],ECI="rcs3",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
  olrsim[i,1:length(olrrun[[1]])]=olrrun[[1]]
  acprun=acpfunc(datafit=datatr[[1]],dataval=datapo[[1]],ECI="rcs3",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
  acpsim[i,1:length(acprun[[1]])]=acprun[[1]]
  smrun=smfunc(datafit=datatr[[1]],dataval=datapo[[1]],ECI="rcs3",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
  smsim[i,1:length(smrun[[1]])]=smrun[[1]]
  print(i)
}
comsave_smalln(filep3="Ordinal prediction model simulations\\ordsimk3p4_equi_unequal_moderateorc_n100.RData")

# small sample simulations to check overfitting (N=500)
set.seed(4682)
mlrsim=matrix(data=NA,nrow=200,ncol=27)
olrsim=matrix(data=NA,nrow=200,ncol=27)
acpsim=matrix(data=NA,nrow=200,ncol=27)
smsim=matrix(data=NA,nrow=200,ncol=27)
for(i in 1:200) {
  datatr=getdata_k3p4(mu_x1=c(0,0.4,0.8),
                      mu_x2=c(0,0.3,0.6),
                      mu_x3=c(0,0.4,0.8),
                      mu_x4=c(0,0.3,0.6),
                      yprev=c(55,30,15),
                      mult=5)
  mlrrun=mlrfunc(datafit=datatr[[1]],dataval=datapo[[1]],ECI="rcs3",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
  mlrsim[i,1:length(mlrrun[[1]])]=mlrrun[[1]]
  olrrun=olrfunc(datafit=datatr[[1]],dataval=datapo[[1]],ECI="rcs3",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
  olrsim[i,1:length(olrrun[[1]])]=olrrun[[1]]
  acprun=acpfunc(datafit=datatr[[1]],dataval=datapo[[1]],ECI="rcs3",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
  acpsim[i,1:length(acprun[[1]])]=acprun[[1]]
  smrun=smfunc(datafit=datatr[[1]],dataval=datapo[[1]],ECI="rcs3",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
  smsim[i,1:length(smrun[[1]])]=smrun[[1]]
  print(i)
}
comsave_smalln(filep3="Ordinal prediction model simulations\\ordsimk3p4_equi_unequal_moderateorc_n500.RData")

}


# SCENARIO 3

{
# Get large dataset
set.seed(34634) # the seed for the large datasets
datapo=getdata_k3p4(mu_x1=c(0,0.7,0.8),
                    mu_x2=c(0,0.6,0.6),
                    mu_x3=c(0,0.5,0.8),
                    mu_x4=c(0,0.1,0.6),
                    yprev=c(33.34,33.33,33.33),#yprev=c(55,30,15),
                    mult=2000)

# Run models for large dataset
mlrfit=vglm(y ~ x1 + x2 + x3 + x4, family=multinomial(refLevel = "1"), data=datapo[[1]])
acnpres=mlrfunc(datafit=datapo[[1]],dataval=datapo[[1]],ECI="nonpar",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
olrres=olrfunc(datafit=datapo[[1]],dataval=datapo[[1]],ECI="nonpar",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
acpres=acpfunc(datafit=datapo[[1]],dataval=datapo[[1]],ECI="nonpar",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
smres=smfunc(datafit=datapo[[1]],dataval=datapo[[1]],ECI="nonpar",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
comsave(filep="Ordinal prediction model simulations\\ordsimk3p4_nonequi_equal_moderateorc.RData")
vgamsmps4a = vgam(datapo[[1]]$y ~ sm.ps(log(acnpres[[3]][,2]/acnpres[[3]][,1]),df=4) + sm.ps(log(acnpres[[3]][,3]/acnpres[[3]][,1]),df=4),family=multinomial(refLevel = 1))
vgamsmps4b = vgam(datapo[[1]]$y ~ sm.ps(log(olrres[[3]][,2]/olrres[[3]][,1]),df=4) + sm.ps(log(olrres[[3]][,3]/olrres[[3]][,1]),df=4),family=multinomial(refLevel = 1))
vgamsmps4c = vgam(datapo[[1]]$y ~ sm.ps(log(acpres[[3]][,2]/acpres[[3]][,1]),df=4) + sm.ps(log(acpres[[3]][,3]/acpres[[3]][,1]),df=4),family=multinomial(refLevel = 1))
vgamsmps4d = vgam(datapo[[1]]$y ~ sm.ps(log(smres[[3]][,2]/smres[[3]][,1]),df=4) + sm.ps(log(smres[[3]][,3]/smres[[3]][,1]),df=4),family=multinomial(refLevel = 1))
plotsave(filep2="Ordinal prediction model simulations\\ordsimk3p4_nonequi_equal_moderateorc.pdf")

# Smooth calibration plots for paper
par(mfrow=c(2,2))
plotlines3sub(preds=acnpres[[3]],obs=vgamsmps4a,mtxt=NULL,subsample0=1000,cexlab=1.4,showleg=0,showtext=1,texttoshow="MLR",cextext=1.3)
plotlines3sub(preds=olrres[[3]],obs=vgamsmps4b,mtxt=NULL,subsample0=1000,spar2=1.2,cexlab=1.4,showleg=0,showtext=1,texttoshow="CL-PO",cextext=1.3)
plotlines3sub(preds=acpres[[3]],obs=vgamsmps4c,mtxt=NULL,subsample0=1000,spar2=1.4,cexlab=1.4,showleg=0,showtext=1,texttoshow="AC-PO",cextext=1.3)
plotlines3sub(preds=smres[[3]],obs=vgamsmps4d,mtxt=NULL,subsample0=1000,cexlab=1.4,showleg=0,showtext=1,texttoshow="SLM",cextext=1.3)
# Scatter plot vs true risk for paper
plotvstruerisk3(subsample=1000)

# dev.off()
# olrx <- orm(y ~ x1 + x2 + x3 +x4, data=datapo[[1]], x=T)
# par(mfrow=c(2,2))
# resid(olrx, 'partial',pl=T,label.curves=F) 

# small sample simulations to check overfitting (N=100)
set.seed(4682)
mlrsim=matrix(data=NA,nrow=200,ncol=27)
olrsim=matrix(data=NA,nrow=200,ncol=27)
acpsim=matrix(data=NA,nrow=200,ncol=27)
smsim=matrix(data=NA,nrow=200,ncol=27)
for(i in 1:200) {
  datatr=getdata_k3p4(mu_x1=c(0,0.7,0.8),
                      mu_x2=c(0,0.6,0.6),
                      mu_x3=c(0,0.5,0.8),
                      mu_x4=c(0,0.1,0.6),
                      yprev=c(33.34,33.33,33.33),
                      mult=1)
  mlrrun=mlrfunc(datafit=datatr[[1]],dataval=datapo[[1]],ECI="rcs3",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
  mlrsim[i,1:length(mlrrun[[1]])]=mlrrun[[1]]
  olrrun=olrfunc(datafit=datatr[[1]],dataval=datapo[[1]],ECI="rcs3",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
  olrsim[i,1:length(olrrun[[1]])]=olrrun[[1]]
  acprun=acpfunc(datafit=datatr[[1]],dataval=datapo[[1]],ECI="rcs3",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
  acpsim[i,1:length(acprun[[1]])]=acprun[[1]]
  smrun=smfunc(datafit=datatr[[1]],dataval=datapo[[1]],ECI="rcs3",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
  smsim[i,1:length(smrun[[1]])]=smrun[[1]]
  print(i)
}
comsave_smalln(filep3="Ordinal prediction model simulations\\ordsimk3p4_nonequi_equal_moderateorc_n100.RData")

# small sample simulations to check overfitting (N=500)
set.seed(4682)
mlrsim=matrix(data=NA,nrow=200,ncol=27)
olrsim=matrix(data=NA,nrow=200,ncol=27)
acpsim=matrix(data=NA,nrow=200,ncol=27)
smsim=matrix(data=NA,nrow=200,ncol=27)
for(i in 1:200) {
  datatr=getdata_k3p4(mu_x1=c(0,0.7,0.8),
                      mu_x2=c(0,0.6,0.6),
                      mu_x3=c(0,0.5,0.8),
                      mu_x4=c(0,0.1,0.6),
                      yprev=c(33.34,33.33,33.33),
                      mult=5)
  mlrrun=mlrfunc(datafit=datatr[[1]],dataval=datapo[[1]],ECI="rcs3",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
  mlrsim[i,1:length(mlrrun[[1]])]=mlrrun[[1]]
  olrrun=olrfunc(datafit=datatr[[1]],dataval=datapo[[1]],ECI="rcs3",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
  olrsim[i,1:length(olrrun[[1]])]=olrrun[[1]]
  acprun=acpfunc(datafit=datatr[[1]],dataval=datapo[[1]],ECI="rcs3",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
  acpsim[i,1:length(acprun[[1]])]=acprun[[1]]
  smrun=smfunc(datafit=datatr[[1]],dataval=datapo[[1]],ECI="rcs3",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
  smsim[i,1:length(smrun[[1]])]=smrun[[1]]
  print(i)
}
comsave_smalln(filep3="Ordinal prediction model simulations\\ordsimk3p4_nonequi_equal_moderateorc_n500.RData")

}


# SCENARIO 4

{
# Get large dataset
set.seed(34634) # the seed for the large datasets
datapo=getdata_k3p4(mu_x1=c(0,0.7,0.8),
                    mu_x2=c(0,0.6,0.6),
                    mu_x3=c(0,0.5,0.8),
                    mu_x4=c(0,0.1,0.6),
                    yprev=c(55,30,15),#yprev=c(55,30,15),
                    mult=2000)

# Run models for large dataset
mlrfit=vglm(y ~ x1 + x2 + x3 + x4, family=multinomial(refLevel = "1"), data=datapo[[1]])
acnpres=mlrfunc(datafit=datapo[[1]],dataval=datapo[[1]],ECI="nonpar",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
olrres=olrfunc(datafit=datapo[[1]],dataval=datapo[[1]],ECI="nonpar",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
acpres=acpfunc(datafit=datapo[[1]],dataval=datapo[[1]],ECI="nonpar",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
smres=smfunc(datafit=datapo[[1]],dataval=datapo[[1]],ECI="nonpar",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
comsave(filep="Ordinal prediction model simulations\\ordsimk3p4_nonequi_unequal_moderateorc.RData")
vgamsmps4a = vgam(datapo[[1]]$y ~ sm.ps(log(acnpres[[3]][,2]/acnpres[[3]][,1]),df=4) + sm.ps(log(acnpres[[3]][,3]/acnpres[[3]][,1]),df=4),family=multinomial(refLevel = 1))
vgamsmps4b = vgam(datapo[[1]]$y ~ sm.ps(log(olrres[[3]][,2]/olrres[[3]][,1]),df=4) + sm.ps(log(olrres[[3]][,3]/olrres[[3]][,1]),df=4),family=multinomial(refLevel = 1))
vgamsmps4c = vgam(datapo[[1]]$y ~ sm.ps(log(acpres[[3]][,2]/acpres[[3]][,1]),df=4) + sm.ps(log(acpres[[3]][,3]/acpres[[3]][,1]),df=4),family=multinomial(refLevel = 1))
vgamsmps4d = vgam(datapo[[1]]$y ~ sm.ps(log(smres[[3]][,2]/smres[[3]][,1]),df=4) + sm.ps(log(smres[[3]][,3]/smres[[3]][,1]),df=4),family=multinomial(refLevel = 1))
plotsave(filep2="Ordinal prediction model simulations\\ordsimk3p4_nonequi_unequal_moderateorc.pdf")

# Smooth calibration plots for paper
par(mfrow=c(2,2))
plotlines3sub(preds=acnpres[[3]],obs=vgamsmps4a,mtxt=NULL,subsample0=1000,cexlab=1.4,showleg=0,showtext=1,texttoshow="MLR",cextext=1.3)
plotlines3sub(preds=olrres[[3]],obs=vgamsmps4b,mtxt=NULL,subsample0=1000,cexlab=1.4,showleg=0,showtext=1,texttoshow="CL-PO",cextext=1.3)
plotlines3sub(preds=acpres[[3]],obs=vgamsmps4c,mtxt=NULL,subsample0=1000,cexlab=1.4,showleg=0,showtext=1,texttoshow="AC-PO",cextext=1.3)
plotlines3sub(preds=smres[[3]],obs=vgamsmps4d,mtxt=NULL,subsample0=1000,cexlab=1.4,showleg=0,showtext=1,texttoshow="SLM",cextext=1.3)
# Scatter plot vs true risk for paper
plotvstruerisk3(subsample=1000)

# dev.off()
# olrx <- orm(y ~ x1 + x2 + x3 +x4, data=datapo[[1]], x=T)
# par(mfrow=c(2,2))
# resid(olrx, 'partial',pl=T,label.curves=F) 

# small sample simulations to check overfitting (N=100)
set.seed(4682)
mlrsim=matrix(data=NA,nrow=200,ncol=27)
olrsim=matrix(data=NA,nrow=200,ncol=27)
acpsim=matrix(data=NA,nrow=200,ncol=27)
smsim=matrix(data=NA,nrow=200,ncol=27)
for(i in 1:200) {
  datatr=getdata_k3p4(mu_x1=c(0,0.7,0.8),
                      mu_x2=c(0,0.6,0.6),
                      mu_x3=c(0,0.5,0.8),
                      mu_x4=c(0,0.1,0.6),
                      yprev=c(55,30,15),
                      mult=1)
  mlrrun=mlrfunc(datafit=datatr[[1]],dataval=datapo[[1]],ECI="rcs3",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
  mlrsim[i,1:length(mlrrun[[1]])]=mlrrun[[1]]
  olrrun=olrfunc(datafit=datatr[[1]],dataval=datapo[[1]],ECI="rcs3",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
  olrsim[i,1:length(olrrun[[1]])]=olrrun[[1]]
  acprun=acpfunc(datafit=datatr[[1]],dataval=datapo[[1]],ECI="rcs3",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
  acpsim[i,1:length(acprun[[1]])]=acprun[[1]]
  smrun=smfunc(datafit=datatr[[1]],dataval=datapo[[1]],ECI="rcs3",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
  smsim[i,1:length(smrun[[1]])]=smrun[[1]]
  print(i)
}
comsave_smalln(filep3="Ordinal prediction model simulations\\ordsimk3p4_nonequi_unequal_moderateorc_n100.RData")

# small sample simulations to check overfitting (N=500)
set.seed(4682)
mlrsim=matrix(data=NA,nrow=200,ncol=27)
olrsim=matrix(data=NA,nrow=200,ncol=27)
acpsim=matrix(data=NA,nrow=200,ncol=27)
smsim=matrix(data=NA,nrow=200,ncol=27)
for(i in 1:200) {
  datatr=getdata_k3p4(mu_x1=c(0,0.7,0.8),
                      mu_x2=c(0,0.6,0.6),
                      mu_x3=c(0,0.5,0.8),
                      mu_x4=c(0,0.1,0.6),
                      yprev=c(55,30,15),
                      mult=5)
  mlrrun=mlrfunc(datafit=datatr[[1]],dataval=datapo[[1]],ECI="rcs3",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
  mlrsim[i,1:length(mlrrun[[1]])]=mlrrun[[1]]
  olrrun=olrfunc(datafit=datatr[[1]],dataval=datapo[[1]],ECI="rcs3",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
  olrsim[i,1:length(olrrun[[1]])]=olrrun[[1]]
  acprun=acpfunc(datafit=datatr[[1]],dataval=datapo[[1]],ECI="rcs3",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
  acpsim[i,1:length(acprun[[1]])]=acprun[[1]]
  smrun=smfunc(datafit=datatr[[1]],dataval=datapo[[1]],ECI="rcs3",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
  smsim[i,1:length(smrun[[1]])]=smrun[[1]]
  print(i)
}
comsave_smalln(filep3="Ordinal prediction model simulations\\ordsimk3p4_nonequi_unequal_moderateorc_n500.RData")

}


# SCENARIO 5

{
# Get large dataset
set.seed(34634) # the seed for the large datasets
datapo=getdata_k3p4(mu_x1=c(0,0.7,0.8),
                    mu_x2=c(0,0.7,0.6),
                    mu_x3=c(0,0,1),
                    mu_x4=c(0.3,0,0.3),
                    yprev=c(55,30,15),#yprev=c(55,30,15),
                    mult=2000)

# Run models for large dataset
mlrfit=vglm(y ~ x1 + x2 + x3 + x4, family=multinomial(refLevel = "1"), data=datapo[[1]])
acnpres=mlrfunc(datafit=datapo[[1]],dataval=datapo[[1]],ECI="nonpar",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
olrres=olrfunc(datafit=datapo[[1]],dataval=datapo[[1]],ECI="nonpar",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
acpres=acpfunc(datafit=datapo[[1]],dataval=datapo[[1]],ECI="nonpar",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
smres=smfunc(datafit=datapo[[1]],dataval=datapo[[1]],ECI="nonpar",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
comsave(filep="Ordinal prediction model simulations\\ordsimk3p4_xnonequi_unequal_moderateorc.RData")
vgamsmps4a = vgam(datapo[[1]]$y ~ sm.ps(log(acnpres[[3]][,2]/acnpres[[3]][,1]),df=4) + sm.ps(log(acnpres[[3]][,3]/acnpres[[3]][,1]),df=4),family=multinomial(refLevel = 1))
vgamsmps4b = vgam(datapo[[1]]$y ~ sm.ps(log(olrres[[3]][,2]/olrres[[3]][,1]),df=4) + sm.ps(log(olrres[[3]][,3]/olrres[[3]][,1]),df=4),family=multinomial(refLevel = 1))
vgamsmps4c = vgam(datapo[[1]]$y ~ sm.ps(log(acpres[[3]][,2]/acpres[[3]][,1]),df=4) + sm.ps(log(acpres[[3]][,3]/acpres[[3]][,1]),df=4),family=multinomial(refLevel = 1))
vgamsmps4d = vgam(datapo[[1]]$y ~ sm.ps(log(smres[[3]][,2]/smres[[3]][,1]),df=4) + sm.ps(log(smres[[3]][,3]/smres[[3]][,1]),df=4),family=multinomial(refLevel = 1))
plotsave(filep2="Ordinal prediction model simulations\\ordsimk3p4_xnonequi_unequal_moderateorc.pdf")

# Smooth calibration plots for paper
par(mfrow=c(2,2))
plotlines3sub(preds=acnpres[[3]],obs=vgamsmps4a,mtxt=NULL,subsample0=1000,cexlab=1.4,showleg=0,showtext=1,texttoshow="MLR",cextext=1.3)
plotlines3sub(preds=olrres[[3]],obs=vgamsmps4b,mtxt=NULL,subsample0=1000,cexlab=1.4,showleg=0,showtext=1,texttoshow="CL-PO",cextext=1.3)
plotlines3sub(preds=acpres[[3]],obs=vgamsmps4c,mtxt=NULL,subsample0=1000,cexlab=1.4,showleg=0,showtext=1,texttoshow="AC-PO",cextext=1.3)
plotlines3sub(preds=smres[[3]],obs=vgamsmps4d,mtxt=NULL,subsample0=1000,cexlab=1.4,showleg=0,showtext=1,texttoshow="SLM",cextext=1.3)
# Scatter plot vs true risk for paper
plotvstruerisk3(subsample=1000)

# dev.off()
# olrx <- orm(y ~ x1 + x2 + x3 +x4, data=datapo[[1]], x=T)
# par(mfrow=c(2,2))
# resid(olrx, 'partial',pl=T,label.curves=F) 

}


# SCENARIO 6

{
# Get large dataset
set.seed(34634) # the seed for the large datasets
datapo=getdata_k4p3(mu_x1=c(0,0,1,1),
                    mu_x2=c(0,0.8,0.8,0.9),
                    mu_x3=c(0.2,0,0.9,1.0),
                    yprev=c(40,25,20,15),#yprev=c(55,30,15),
                    mult=2000)

# Run models for large dataset
mlrfit=vglm(y ~ x1 + x2 + x3, family=multinomial(refLevel = "1"), data=datapo[[1]])
acnpres=mlrfunc(datafit=datapo[[1]],dataval=datapo[[1]],ECI="nonpar",levels=4,modelformula=y ~ x1 + x2 + x3)
olrres=olrfunc(datafit=datapo[[1]],dataval=datapo[[1]],ECI="nonpar",levels=4,modelformula=y ~ x1 + x2 + x3)
acpres=acpfunc(datafit=datapo[[1]],dataval=datapo[[1]],ECI="nonpar",levels=4,modelformula=y ~ x1 + x2 + x3)
smres=smfunc(datafit=datapo[[1]],dataval=datapo[[1]],ECI="nonpar",levels=4,modelformula=y ~ x1 + x2 + x3)
comsave(filep="Ordinal prediction model simulations\\ordsimk4p3_xnonequi_unequal_moderateorc.RData")
vgamsmps4a = vgam(datapo[[1]]$y ~ sm.ps(log(acnpres[[3]][,2]/acnpres[[3]][,1]),df=4) + sm.ps(log(acnpres[[3]][,3]/acnpres[[3]][,1]),df=4) + sm.ps(log(acnpres[[3]][,4]/acnpres[[3]][,1]),df=4),family=multinomial(refLevel = 1))
vgamsmps4b = vgam(datapo[[1]]$y ~ sm.ps(log(olrres[[3]][,2]/olrres[[3]][,1]),df=4) + sm.ps(log(olrres[[3]][,3]/olrres[[3]][,1]),df=4) + sm.ps(log(olrres[[3]][,4]/olrres[[3]][,1]),df=4),family=multinomial(refLevel = 1))
vgamsmps4c = vgam(datapo[[1]]$y ~ sm.ps(log(acpres[[3]][,2]/acpres[[3]][,1]),df=4) + sm.ps(log(acpres[[3]][,3]/acpres[[3]][,1]),df=4) + sm.ps(log(acpres[[3]][,4]/acpres[[3]][,1]),df=4),family=multinomial(refLevel = 1))
vgamsmps4d = vgam(datapo[[1]]$y ~ sm.ps(log(smres[[3]][,2]/smres[[3]][,1]),df=4) + sm.ps(log(smres[[3]][,3]/smres[[3]][,1]),df=4) + sm.ps(log(smres[[3]][,4]/smres[[3]][,1]),df=4),family=multinomial(refLevel = 1))
plotsave(filep2="Ordinal prediction model simulations\\ordsimk4p3_xnonequi_unequal_moderateorc.pdf")

# Smooth calibration plots for paper
par(mfrow=c(2,2))
plotlines4sub(preds=acnpres[[3]],obs=vgamsmps4a,mtxt=NULL,subsample0=1000,cexlab=1.4,showleg=0,showtext=1,texttoshow="MLR",cextext=1.3)
plotlines4sub(preds=olrres[[3]],obs=vgamsmps4b,mtxt=NULL,subsample0=1000,spar2=1.4,cexlab=1.4,showleg=0,showtext=1,texttoshow="CL-PO",cextext=1.3)
plotlines4sub(preds=acpres[[3]],obs=vgamsmps4c,mtxt=NULL,subsample0=1000,spar2=1.2,cexlab=1.4,showleg=0,showtext=1,texttoshow="AC-PO",cextext=1.3)
plotlines4sub(preds=smres[[3]],obs=vgamsmps4d,mtxt=NULL,subsample0=1000,cexlab=1.4,showleg=0,showtext=1,texttoshow="SLM",cextext=1.3)
# Scatter plot vs true risk for paper
plotvstruerisk4(subsample=1000)

# dev.off()
# olrx <- orm(y ~ x1 + x2 + x3, data=datapo[[1]], x=T)
# par(mfrow=c(2,2))
# resid(olrx, 'partial',pl=T,label.curves=F) 

}


# SCENARIO 7

{
# Get large dataset
set.seed(34634) # the seed for the large datasets
datapo=getdata_k4p3(mu_x1=c(0,0,0.6,0.6),
                    mu_x2=c(0,0.4,0.4,0.5),
                    mu_x3=c(0.1,0,0.6,0.7),
                    yprev=c(40,25,20,15),#yprev=c(55,30,15),
                    mult=2000)

# Run models for large dataset
mlrfit=vglm(y ~ x1 + x2 + x3, family=multinomial(refLevel = "1"), data=datapo[[1]])
acnpres=mlrfunc(datafit=datapo[[1]],dataval=datapo[[1]],ECI="nonpar",levels=4,modelformula=y ~ x1 + x2 + x3)
olrres=olrfunc(datafit=datapo[[1]],dataval=datapo[[1]],ECI="nonpar",levels=4,modelformula=y ~ x1 + x2 + x3)
acpres=acpfunc(datafit=datapo[[1]],dataval=datapo[[1]],ECI="nonpar",levels=4,modelformula=y ~ x1 + x2 + x3)
smres=smfunc(datafit=datapo[[1]],dataval=datapo[[1]],ECI="nonpar",levels=4,modelformula=y ~ x1 + x2 + x3)
comsave(filep="Ordinal prediction model simulations\\ordsimk4p3_xnonequi_unequal_loworc.RData")
vgamsmps4a = vgam(datapo[[1]]$y ~ sm.ps(log(acnpres[[3]][,2]/acnpres[[3]][,1]),df=4) + sm.ps(log(acnpres[[3]][,3]/acnpres[[3]][,1]),df=4) + sm.ps(log(acnpres[[3]][,4]/acnpres[[3]][,1]),df=4),family=multinomial(refLevel = 1))
vgamsmps4b = vgam(datapo[[1]]$y ~ sm.ps(log(olrres[[3]][,2]/olrres[[3]][,1]),df=4) + sm.ps(log(olrres[[3]][,3]/olrres[[3]][,1]),df=4) + sm.ps(log(olrres[[3]][,4]/olrres[[3]][,1]),df=4),family=multinomial(refLevel = 1))
vgamsmps4c = vgam(datapo[[1]]$y ~ sm.ps(log(acpres[[3]][,2]/acpres[[3]][,1]),df=4) + sm.ps(log(acpres[[3]][,3]/acpres[[3]][,1]),df=4) + sm.ps(log(acpres[[3]][,4]/acpres[[3]][,1]),df=4),family=multinomial(refLevel = 1))
vgamsmps4d = vgam(datapo[[1]]$y ~ sm.ps(log(smres[[3]][,2]/smres[[3]][,1]),df=4) + sm.ps(log(smres[[3]][,3]/smres[[3]][,1]),df=4) + sm.ps(log(smres[[3]][,4]/smres[[3]][,1]),df=4),family=multinomial(refLevel = 1))
plotsave(filep2="Ordinal prediction model simulations\\ordsimk4p3_xnonequi_unequal_loworc.pdf")

# Smooth calibration plots for paper
par(mfrow=c(2,2))
plotlines4sub(preds=acnpres[[3]],obs=vgamsmps4a,mtxt=NULL,subsample0=1000,cexlab=1.4,showleg=0,showtext=1,texttoshow="MLR",cextext=1.3)
plotlines4sub(preds=olrres[[3]],obs=vgamsmps4b,mtxt=NULL,subsample0=1000,spar2=1.5,cexlab=1.4,showleg=0,showtext=1,texttoshow="CL-PO",cextext=1.3)
plotlines4sub(preds=acpres[[3]],obs=vgamsmps4c,mtxt=NULL,subsample0=1000,spar2=1.5,cexlab=1.4,showleg=0,showtext=1,texttoshow="AC-PO",cextext=1.3)
plotlines4sub(preds=smres[[3]],obs=vgamsmps4d,mtxt=NULL,subsample0=1000,cexlab=1.4,showleg=0,showtext=1,texttoshow="SLM",cextext=1.3)
# Scatter plot vs true risk for paper
plotvstruerisk4(subsample=1000)

# dev.off()
# olrx <- orm(y ~ x1 + x2 + x3, data=datapo[[1]], x=T)
# par(mfrow=c(2,2))
# resid(olrx, 'partial',pl=T,label.curves=F) 

}


# SCENARIO 8

{
# Get large dataset
set.seed(34634) # the seed for the large datasets
datapo=getdata_k4p3(mu_x1=c(0,0,0.6,0.6),
                    mu_x2=c(0,0.4,0.4,0.5),
                    mu_x3=c(0.1,0,0.6,0.7),
                    yprev=c(45,30,20,5),#yprev=c(55,30,15),
                    mult=2000)

# Run models for large dataset
mlrfit=vglm(y ~ x1 + x2 + x3, family=multinomial(refLevel = "1"), data=datapo[[1]])
acnpres=mlrfunc(datafit=datapo[[1]],dataval=datapo[[1]],ECI="nonpar",levels=4,modelformula=y ~ x1 + x2 + x3)
olrres=olrfunc(datafit=datapo[[1]],dataval=datapo[[1]],ECI="nonpar",levels=4,modelformula=y ~ x1 + x2 + x3)
acpres=acpfunc(datafit=datapo[[1]],dataval=datapo[[1]],ECI="nonpar",levels=4,modelformula=y ~ x1 + x2 + x3)
smres=smfunc(datafit=datapo[[1]],dataval=datapo[[1]],ECI="nonpar",levels=4,modelformula=y ~ x1 + x2 + x3)
comsave(filep="Ordinal prediction model simulations\\ordsimk4p3_xnonequi_xunequal_loworc.RData")
vgamsmps4a = vgam(datapo[[1]]$y ~ sm.ps(log(acnpres[[3]][,2]/acnpres[[3]][,1]),df=4) + sm.ps(log(acnpres[[3]][,3]/acnpres[[3]][,1]),df=4) + sm.ps(log(acnpres[[3]][,4]/acnpres[[3]][,1]),df=4),family=multinomial(refLevel = 1))
vgamsmps4b = vgam(datapo[[1]]$y ~ sm.ps(log(olrres[[3]][,2]/olrres[[3]][,1]),df=4) + sm.ps(log(olrres[[3]][,3]/olrres[[3]][,1]),df=4) + sm.ps(log(olrres[[3]][,4]/olrres[[3]][,1]),df=4),family=multinomial(refLevel = 1))
vgamsmps4c = vgam(datapo[[1]]$y ~ sm.ps(log(acpres[[3]][,2]/acpres[[3]][,1]),df=4) + sm.ps(log(acpres[[3]][,3]/acpres[[3]][,1]),df=4) + sm.ps(log(acpres[[3]][,4]/acpres[[3]][,1]),df=4),family=multinomial(refLevel = 1))
vgamsmps4d = vgam(datapo[[1]]$y ~ sm.ps(log(smres[[3]][,2]/smres[[3]][,1]),df=4) + sm.ps(log(smres[[3]][,3]/smres[[3]][,1]),df=4) + sm.ps(log(smres[[3]][,4]/smres[[3]][,1]),df=4),family=multinomial(refLevel = 1))
plotsave(filep2="Ordinal prediction model simulations\\ordsimk4p3_xnonequi_xunequal_loworc.pdf")

# Smooth calibration plots for paper
par(mfrow=c(2,2))
plotlines4sub(preds=acnpres[[3]],obs=vgamsmps4a,mtxt=NULL,subsample0=1000,cexlab=1.4,showleg=0,showtext=1,texttoshow="MLR",cextext=1.3)
plotlines4sub(preds=olrres[[3]],obs=vgamsmps4b,mtxt=NULL,subsample0=1000,cexlab=1.4,showleg=0,showtext=1,texttoshow="CL-PO",cextext=1.3)
plotlines4sub(preds=acpres[[3]],obs=vgamsmps4c,mtxt=NULL,subsample0=1000,spar2=1.5,cexlab=1.4,showleg=0,showtext=1,texttoshow="AC-PO",cextext=1.3)
plotlines4sub(preds=smres[[3]],obs=vgamsmps4d,mtxt=NULL,subsample0=1000,cexlab=1.4,showleg=0,showtext=1,texttoshow="SLM",cextext=1.3)
# Scatter plot vs true risk for paper
plotvstruerisk4(subsample=1000)

# dev.off()
# olrx <- orm(y ~ x1 + x2 + x3, data=datapo[[1]], x=T)
# olry <- vglm(y ~ x1 + x2 + x3, family=cumulative(parallel=T), data=datapo[[1]])
# par(mfrow=c(2,2))
# resid(olrx, 'partial',pl=T,label.curves=F) 

}


# SCENARIO 9

{
# Get large dataset
set.seed(34634) # the seed for the large datasets
datapo=getdata_k3p4b(px1=c(0.2,0.55,0.58),
                    px2=c(0.2,0.5,0.5),
                    px3=c(0.2,0.45,0.58),
                    px4=c(0.2,0.25,0.5),
                    yprev=c(55,30,15),#yprev=c(55,30,15),
                    mult=2000)
CrossTable(datapo[[1]]$y)

# Run models for large dataset
mlrfit=vglm(y ~ x1 + x2 + x3 + x4, family=multinomial(refLevel = "1"), data=datapo[[1]])
acnpres=mlrfunc(datafit=datapo[[1]],dataval=datapo[[1]],ECI="nonpar",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
olrres=olrfunc(datafit=datapo[[1]],dataval=datapo[[1]],ECI="nonpar",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
acpres=acpfunc(datafit=datapo[[1]],dataval=datapo[[1]],ECI="nonpar",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
smres=smfunc(datafit=datapo[[1]],dataval=datapo[[1]],ECI="nonpar",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
comsave(filep="Ordinal prediction model simulations\\ordsimk3p4bin_nonequi_unequal_moderateorc.RData")
vgamsmps4a = vgam(datapo[[1]]$y ~ sm.ps(log(acnpres[[3]][,2]/acnpres[[3]][,1]),df=4) + sm.ps(log(acnpres[[3]][,3]/acnpres[[3]][,1]),df=4),family=multinomial(refLevel = 1))
vgamsmps4b = vgam(datapo[[1]]$y ~ sm.ps(log(olrres[[3]][,2]/olrres[[3]][,1]),df=4) + sm.ps(log(olrres[[3]][,3]/olrres[[3]][,1]),df=4),family=multinomial(refLevel = 1))
vgamsmps4c = vgam(datapo[[1]]$y ~ sm.ps(log(acpres[[3]][,2]/acpres[[3]][,1]),df=4) + sm.ps(log(acpres[[3]][,3]/acpres[[3]][,1]),df=4),family=multinomial(refLevel = 1))
vgamsmps4d = vgam(datapo[[1]]$y ~ sm.ps(log(smres[[3]][,2]/smres[[3]][,1]),df=4) + sm.ps(log(smres[[3]][,3]/smres[[3]][,1]),df=4),family=multinomial(refLevel = 1))
plotsavebin(filep2="Ordinal prediction model simulations\\ordsimk3p4bin_nonequi_unequal_moderateorc.pdf")

# Smooth calibration plots for paper
par(mfrow=c(2,2))
plotscatter3sub(preds=acnpres[[3]],obs=vgamsmps4a,mtxt=NULL,subsample0=1000,cexlab=1.4,showleg=0,showtext=1,texttoshow="MLR",cextext=1.3)
plotscatter3sub(preds=olrres[[3]],obs=vgamsmps4b,mtxt=NULL,subsample0=1000,cexlab=1.4,showleg=0,showtext=1,texttoshow="CL-PO",cextext=1.3)
plotscatter3sub(preds=acpres[[3]],obs=vgamsmps4c,mtxt=NULL,subsample0=1000,cexlab=1.4,showleg=0,showtext=1,texttoshow="AC-PO",cextext=1.3)
plotscatter3sub(preds=smres[[3]],obs=vgamsmps4d,mtxt=NULL,subsample0=1000,cexlab=1.4,showleg=0,showtext=1,texttoshow="SLM",cextext=1.3)
# Scatter plot vs true risk for paper
plotvstruerisk3(subsample=1000,cexpts=0.8)

# dev.off()
# olrx <- orm(y ~ x1 + x2 + x3 + x4, data=datapo[[1]], x=T)
# par(mfrow=c(2,2))
# resid(olrx, 'partial',pl=T,label.curves=F) 

}


# SCENARIO 10

{
# Get large dataset
set.seed(34634) # the seed for the large datasets
datapo=getdata_k4p3b(px1=c(0.2,0.2,0.65,0.65),
                     px2=c(0.2,0.4,0.4,0.6),
                     px3=c(0.25,0.2,0.6,0.7),
                     yprev=c(40,25,20,15),
                     mult=2000)

# Run models for large dataset
mlrfit=vglm(y ~ x1 + x2 + x3, family=multinomial(refLevel = "1"), data=datapo[[1]])
acnpres=mlrfunc(datafit=datapo[[1]],dataval=datapo[[1]],ECI="nonpar",levels=4,modelformula=y ~ x1 + x2 + x3)
olrres=olrfunc(datafit=datapo[[1]],dataval=datapo[[1]],ECI="nonpar",levels=4,modelformula=y ~ x1 + x2 + x3)
acpres=acpfunc(datafit=datapo[[1]],dataval=datapo[[1]],ECI="nonpar",levels=4,modelformula=y ~ x1 + x2 + x3)
smres=smfunc(datafit=datapo[[1]],dataval=datapo[[1]],ECI="nonpar",levels=4,modelformula=y ~ x1 + x2 + x3)
comsave(filep="Ordinal prediction model simulations\\ordsimk4p3bin_xnonequi_unequal_moderateorc.RData")
vgamsmps4a = vgam(datapo[[1]]$y ~ sm.ps(log(acnpres[[3]][,2]/acnpres[[3]][,1]),df=4) + sm.ps(log(acnpres[[3]][,3]/acnpres[[3]][,1]),df=4) + sm.ps(log(acnpres[[3]][,4]/acnpres[[3]][,1]),df=4),family=multinomial(refLevel = 1))
vgamsmps4b = vgam(datapo[[1]]$y ~ sm.ps(log(olrres[[3]][,2]/olrres[[3]][,1]),df=4) + sm.ps(log(olrres[[3]][,3]/olrres[[3]][,1]),df=4) + sm.ps(log(olrres[[3]][,4]/olrres[[3]][,1]),df=4),family=multinomial(refLevel = 1))
vgamsmps4c = vgam(datapo[[1]]$y ~ sm.ps(log(acpres[[3]][,2]/acpres[[3]][,1]),df=4) + sm.ps(log(acpres[[3]][,3]/acpres[[3]][,1]),df=4) + sm.ps(log(acpres[[3]][,4]/acpres[[3]][,1]),df=4),family=multinomial(refLevel = 1))
vgamsmps4d = vgam(datapo[[1]]$y ~ sm.ps(log(smres[[3]][,2]/smres[[3]][,1]),df=4) + sm.ps(log(smres[[3]][,3]/smres[[3]][,1]),df=4) + sm.ps(log(smres[[3]][,4]/smres[[3]][,1]),df=4),family=multinomial(refLevel = 1))
plotsavebin(filep2="Ordinal prediction model simulations\\ordsimk4p3bin_xnonequi_unequal_moderateorc.pdf")

# Smooth calibration plots for paper
par(mfrow=c(2,2))
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4a,mtxt=NULL,subsample0=1000,cexlab=1.4,showleg=0,showtext=1,texttoshow="MLR",cextext=1.3)
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4b,mtxt=NULL,subsample0=1000,cexlab=1.4,showleg=0,showtext=1,texttoshow="CL-PO",cextext=1.3)
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4c,mtxt=NULL,subsample0=1000,cexlab=1.4,showleg=0,showtext=1,texttoshow="AC-PO",cextext=1.3)
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4d,mtxt=NULL,subsample0=1000,cexlab=1.4,showleg=0,showtext=1,texttoshow="SLM",cextext=1.3)
# Scatter plot vs true risk for paper
plotvstruerisk4(subsample=1000,cexpts=0.8)

# dev.off()
# olrx <- orm(y ~ x1 + x2 + x3, data=datapo[[1]], x=T)
# par(mfrow=c(2,2))
# resid(olrx, 'partial',pl=T,label.curves=F) 

}


# LOAD RESULTS OF THE 10 SCENARIOS AND SUMMARIZE THEM

load("Ordinal prediction model simulations\\MLR truth Scenario 1\\ordsimk3p4_equi_equal_moderateorc.rdata")
load("Ordinal prediction model simulations\\MLR truth Scenario 2\\ordsimk3p4_equi_unequal_moderateorc.rdata")
load("Ordinal prediction model simulations\\MLR truth Scenario 3\\ordsimk3p4_nonequi_equal_moderateorc.rdata")
load("Ordinal prediction model simulations\\MLR truth Scenario 4\\ordsimk3p4_nonequi_unequal_moderateorc.rdata")
load("Ordinal prediction model simulations\\MLR truth Scenario 5\\ordsimk3p4_xnonequi_unequal_moderateorc.RData")
load("Ordinal prediction model simulations\\MLR truth Scenario 6\\ordsimk4p3_xnonequi_unequal_moderateorc.RData")
load("Ordinal prediction model simulations\\MLR truth Scenario 7\\ordsimk4p3_xnonequi_unequal_loworc.RData")
load("Ordinal prediction model simulations\\MLR truth Scenario 8\\ordsimk4p3_xnonequi_xunequal_loworc.RData")
load("Ordinal prediction model simulations\\MLR truth Scenario 9\\ordsimk3p4bin_nonequi_unequal_moderateorc.rdata")
load("Ordinal prediction model simulations\\MLR truth Scenario 10\\ordsimk4p3bin_xnonequi_unequal_moderateorc.rdata")

load("Ordinal prediction model simulations\\MLR truth Scenario 1\\ordsimk3p4_equi_equal_moderateorc_n100.rdata")
load("Ordinal prediction model simulations\\MLR truth Scenario 2\\ordsimk3p4_equi_unequal_moderateorc_n100.rdata")
load("Ordinal prediction model simulations\\MLR truth Scenario 3\\ordsimk3p4_nonequi_equal_moderateorc_n100.rdata")
load("Ordinal prediction model simulations\\MLR truth Scenario 4\\ordsimk3p4_nonequi_unequal_moderateorc_n100.rdata")

load("Ordinal prediction model simulations\\MLR truth Scenario 1\\ordsimk3p4_equi_equal_moderateorc_n500.rdata")
load("Ordinal prediction model simulations\\MLR truth Scenario 2\\ordsimk3p4_equi_unequal_moderateorc_n500.rdata")
load("Ordinal prediction model simulations\\MLR truth Scenario 3\\ordsimk3p4_nonequi_equal_moderateorc_n500.rdata")
load("Ordinal prediction model simulations\\MLR truth Scenario 4\\ordsimk3p4_nonequi_unequal_moderateorc_n500.rdata")



#################################
# SIMULATIONS UNDER CL-PO TRUTH #
#################################


# SCENARIO 1

{
# Get large dataset PO
set.seed(34634) # the seed for the large datasets
getdata_k3p4po <- function(mult){
  x1=rnorm(100*mult,0.4,1.21)
  x2=rnorm(100*mult,0.3,1.11)
  x3=rnorm(100*mult,0.4,1.21)
  x4=rnorm(100*mult,0.3,1.11)
  lp1 = -0.18 - 0.55*x1 - 0.41*x2 - 0.55*x3 - 0.41*x4
  lp2 = 1.55 - 0.55*x1 - 0.41*x2 - 0.55*x3 - 0.41*x4
  ptrue1 = exp(lp1)/(1+exp(lp1))
  ptrue2 = exp(lp2)/(1+exp(lp2)) - ptrue1
  ptrue3 = 1 - ptrue1 - ptrue2
  y=matrix(data=NA,nrow=length(lp1),ncol=1)
  for (i in 1:length(lp1)){
    y[i]=which(grepl(1, rmultinom(1,1,c(ptrue1[i],ptrue2[i],ptrue3[i]))))
  }  
  datax = as.data.frame(cbind(x1,x2,x3,x4,y,ptrue1,ptrue2,ptrue3))
  names(datax)=c("x1","x2","x3","x4","y","ptrue1","ptrue2","ptrue3")
  datax$y = ordered(datax$y)
  return(list(datax,c(mult)))
}
datapo=getdata_k3p4po(mult=2000)

# Run models for large dataset
mlrfit=vglm(y ~ x1 + x2 + x3 + x4, family=multinomial(refLevel = "1"), data=datapo[[1]])
acnpres=mlrfunc(datafit=datapo[[1]],dataval=datapo[[1]],ECI="nonpar",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
olrres=olrfunc(datafit=datapo[[1]],dataval=datapo[[1]],ECI="nonpar",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
acpres=acpfunc(datafit=datapo[[1]],dataval=datapo[[1]],ECI="nonpar",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
smres=smfunc(datafit=datapo[[1]],dataval=datapo[[1]],ECI="nonpar",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
comsave(filep1="Ordinal prediction model simulations\\ordsimk3p4_equi_equal_moderateorc_PO.RData")
vgamsmps4a = vgam(datapo[[1]]$y ~ sm.ps(log(acnpres[[3]][,2]/acnpres[[3]][,1]),df=4) + sm.ps(log(acnpres[[3]][,3]/acnpres[[3]][,1]),df=4),family=multinomial(refLevel = 1))
vgamsmps4b = vgam(datapo[[1]]$y ~ sm.ps(log(olrres[[3]][,2]/olrres[[3]][,1]),df=4) + sm.ps(log(olrres[[3]][,3]/olrres[[3]][,1]),df=4),family=multinomial(refLevel = 1))
vgamsmps4c = vgam(datapo[[1]]$y ~ sm.ps(log(acpres[[3]][,2]/acpres[[3]][,1]),df=4) + sm.ps(log(acpres[[3]][,3]/acpres[[3]][,1]),df=4),family=multinomial(refLevel = 1))
vgamsmps4d = vgam(datapo[[1]]$y ~ sm.ps(log(smres[[3]][,2]/smres[[3]][,1]),df=4) + sm.ps(log(smres[[3]][,3]/smres[[3]][,1]),df=4),family=multinomial(refLevel = 1))
plotsave(filep2="Ordinal prediction model simulations\\ordsimk3p4_equi_equal_moderateorc_PO.pdf")

# Smooth calibration plots for paper
dev.off()
par(mfrow=c(2,2))
plotlines3sub(preds=acnpres[[3]],obs=vgamsmps4a,mtxt=NULL,subsample0=1000,cexlab=1.4,showleg=0,showtext=1,texttoshow="MLR",cextext=1.3)
plotlines3sub(preds=olrres[[3]],obs=vgamsmps4b,mtxt=NULL,subsample0=1000,cexlab=1.4,showleg=0,showtext=1,texttoshow="CL-PO",cextext=1.3)
plotlines3sub(preds=acpres[[3]],obs=vgamsmps4c,mtxt=NULL,subsample0=1000,cexlab=1.4,showleg=0,showtext=1,texttoshow="AC-PO",cextext=1.3)
plotlines3sub(preds=smres[[3]],obs=vgamsmps4d,mtxt=NULL,subsample0=1000,cexlab=1.4,showleg=0,showtext=1,texttoshow="SLM",cextext=1.3)
# Scatter plot vs true risk for paper
plotvstruerisk3(subsample=1000)

# dev.off()
# olrx <- orm(y ~ x1 + x2 + x3 + x4, data=datapo[[1]], x=T)
# par(mfrow=c(2,2))
# resid(olrx, 'partial',pl=T,label.curves=F) 

# small sample simulations to check overfitting (N=100)
set.seed(4682)
mlrsim=matrix(data=NA,nrow=200,ncol=27)
olrsim=matrix(data=NA,nrow=200,ncol=27)
acpsim=matrix(data=NA,nrow=200,ncol=27)
smsim=matrix(data=NA,nrow=200,ncol=27)
for(i in 1:200) {
  datatr=getdata_k3p4po(mult=1)
  mlrrun=mlrfunc(datafit=datatr[[1]],dataval=datapo[[1]],ECI="rcs3",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
  mlrsim[i,1:length(mlrrun[[1]])]=mlrrun[[1]]
  olrrun=olrfunc(datafit=datatr[[1]],dataval=datapo[[1]],ECI="rcs3",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
  olrsim[i,1:length(olrrun[[1]])]=olrrun[[1]]
  acprun=acpfunc(datafit=datatr[[1]],dataval=datapo[[1]],ECI="rcs3",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
  acpsim[i,1:length(acprun[[1]])]=acprun[[1]]
  smrun=smfunc(datafit=datatr[[1]],dataval=datapo[[1]],ECI="rcs3",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
  smsim[i,1:length(smrun[[1]])]=smrun[[1]]
  print(i)
}
comsave_smalln(filep3="Ordinal prediction model simulations\\ordsimk3p4_equi_equal_moderateorc_PO_n100.RData")

# small sample simulations to check overfitting (N=500)
set.seed(4682)
mlrsim=matrix(data=NA,nrow=200,ncol=27)
olrsim=matrix(data=NA,nrow=200,ncol=27)
acpsim=matrix(data=NA,nrow=200,ncol=27)
smsim=matrix(data=NA,nrow=200,ncol=27)
for(i in 1:200) {
  datatr=getdata_k3p4po(mult=5)
  mlrrun=mlrfunc(datafit=datatr[[1]],dataval=datapo[[1]],ECI="rcs3",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
  mlrsim[i,1:length(mlrrun[[1]])]=mlrrun[[1]]
  olrrun=olrfunc(datafit=datatr[[1]],dataval=datapo[[1]],ECI="rcs3",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
  olrsim[i,1:length(olrrun[[1]])]=olrrun[[1]]
  acprun=acpfunc(datafit=datatr[[1]],dataval=datapo[[1]],ECI="rcs3",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
  acpsim[i,1:length(acprun[[1]])]=acprun[[1]]
  smrun=smfunc(datafit=datatr[[1]],dataval=datapo[[1]],ECI="rcs3",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
  smsim[i,1:length(smrun[[1]])]=smrun[[1]]
  print(i)
}
comsave_smalln(filep3="Ordinal prediction model simulations\\ordsimk3p4_equi_equal_moderateorc_PO_n500.RData")

}


# SCENARIO 2

{
# Get large dataset PO
set.seed(34634) # the seed for the large datasets
getdata_k3p4po <- function(mult){
  x1=rnorm(100*mult,0.4,1.28)
  x2=rnorm(100*mult,0.3,1.16)
  x3=rnorm(100*mult,0.4,1.28)
  x4=rnorm(100*mult,0.3,1.16)
  lp1 = 0.92 - 0.53*x1 - 0.39*x2 - 0.53*x3 - 0.39*x4
  lp2 = 2.80 - 0.53*x1 - 0.39*x2 - 0.53*x3 - 0.39*x4
  ptrue1 = exp(lp1)/(1+exp(lp1))
  ptrue2 = exp(lp2)/(1+exp(lp2)) - ptrue1
  ptrue3 = 1 - ptrue1 - ptrue2
  y=matrix(data=NA,nrow=length(lp1),ncol=1)
  for (i in 1:length(lp1)){
    y[i]=which(grepl(1, rmultinom(1,1,c(ptrue1[i],ptrue2[i],ptrue3[i]))))
  }  
  datax = as.data.frame(cbind(x1,x2,x3,x4,y,ptrue1,ptrue2,ptrue3))
  names(datax)=c("x1","x2","x3","x4","y","ptrue1","ptrue2","ptrue3")
  datax$y = ordered(datax$y)
  return(list(datax,c(mult)))
}
datapo=getdata_k3p4po(mult=2000)

# Run models for large dataset
mlrfit=vglm(y ~ x1 + x2 + x3 + x4, family=multinomial(refLevel = "1"), data=datapo[[1]])
acnpres=mlrfunc(datafit=datapo[[1]],dataval=datapo[[1]],ECI="nonpar",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
olrres=olrfunc(datafit=datapo[[1]],dataval=datapo[[1]],ECI="nonpar",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
acpres=acpfunc(datafit=datapo[[1]],dataval=datapo[[1]],ECI="nonpar",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
smres=smfunc(datafit=datapo[[1]],dataval=datapo[[1]],ECI="nonpar",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
comsave(filep1="Ordinal prediction model simulations\\ordsimk3p4_equi_unequal_moderateorc_PO.RData")
vgamsmps4a = vgam(datapo[[1]]$y ~ sm.ps(log(acnpres[[3]][,2]/acnpres[[3]][,1]),df=4) + sm.ps(log(acnpres[[3]][,3]/acnpres[[3]][,1]),df=4),family=multinomial(refLevel = 1))
vgamsmps4b = vgam(datapo[[1]]$y ~ sm.ps(log(olrres[[3]][,2]/olrres[[3]][,1]),df=4) + sm.ps(log(olrres[[3]][,3]/olrres[[3]][,1]),df=4),family=multinomial(refLevel = 1))
vgamsmps4c = vgam(datapo[[1]]$y ~ sm.ps(log(acpres[[3]][,2]/acpres[[3]][,1]),df=4) + sm.ps(log(acpres[[3]][,3]/acpres[[3]][,1]),df=4),family=multinomial(refLevel = 1))
vgamsmps4d = vgam(datapo[[1]]$y ~ sm.ps(log(smres[[3]][,2]/smres[[3]][,1]),df=4) + sm.ps(log(smres[[3]][,3]/smres[[3]][,1]),df=4),family=multinomial(refLevel = 1))
plotsave(filep2="Ordinal prediction model simulations\\ordsimk3p4_equi_unequal_moderateorc_PO.pdf")

# Smooth calibration plots for paper
dev.off()
par(mfrow=c(2,2))
plotlines3sub(preds=acnpres[[3]],obs=vgamsmps4a,mtxt=NULL,subsample0=1000,cexlab=1.4,showleg=0,showtext=1,texttoshow="MLR",cextext=1.3)
plotlines3sub(preds=olrres[[3]],obs=vgamsmps4b,mtxt=NULL,subsample0=1000,cexlab=1.4,showleg=0,showtext=1,texttoshow="CL-PO",cextext=1.3)
plotlines3sub(preds=acpres[[3]],obs=vgamsmps4c,mtxt=NULL,subsample0=1000,cexlab=1.4,showleg=0,showtext=1,texttoshow="AC-PO",cextext=1.3)
plotlines3sub(preds=smres[[3]],obs=vgamsmps4d,mtxt=NULL,subsample0=1000,cexlab=1.4,showleg=0,showtext=1,texttoshow="SLM",cextext=1.3)
# Scatter plot vs true risk for paper
plotvstruerisk3(subsample=1000)

# dev.off()
# olrx <- orm(y ~ x1 + x2 + x3 + x4, data=datapo[[1]], x=T)
# par(mfrow=c(2,2))
# resid(olrx, 'partial',pl=T,label.curves=F) 

# small sample simulations to check overfitting (N=100)
set.seed(4682)
mlrsim=matrix(data=NA,nrow=200,ncol=27)
olrsim=matrix(data=NA,nrow=200,ncol=27)
acpsim=matrix(data=NA,nrow=200,ncol=27)
smsim=matrix(data=NA,nrow=200,ncol=27)
for(i in 1:200) {
  datatr=getdata_k3p4po(mult=1)
  mlrrun=mlrfunc(datafit=datatr[[1]],dataval=datapo[[1]],ECI="rcs3",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
  mlrsim[i,1:length(mlrrun[[1]])]=mlrrun[[1]]
  olrrun=olrfunc(datafit=datatr[[1]],dataval=datapo[[1]],ECI="rcs3",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
  olrsim[i,1:length(olrrun[[1]])]=olrrun[[1]]
  acprun=acpfunc(datafit=datatr[[1]],dataval=datapo[[1]],ECI="rcs3",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
  acpsim[i,1:length(acprun[[1]])]=acprun[[1]]
  smrun=smfunc(datafit=datatr[[1]],dataval=datapo[[1]],ECI="rcs3",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
  smsim[i,1:length(smrun[[1]])]=smrun[[1]]
  print(i)
}
comsave_smalln(filep3="Ordinal prediction model simulations\\ordsimk3p4_equi_unequal_moderateorc_PO_n100.RData")

# small sample simulations to check overfitting (N=500)
set.seed(4682)
mlrsim=matrix(data=NA,nrow=200,ncol=27)
olrsim=matrix(data=NA,nrow=200,ncol=27)
acpsim=matrix(data=NA,nrow=200,ncol=27)
smsim=matrix(data=NA,nrow=200,ncol=27)
for(i in 1:200) {
  datatr=getdata_k3p4po(mult=5)
  mlrrun=mlrfunc(datafit=datatr[[1]],dataval=datapo[[1]],ECI="rcs3",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
  mlrsim[i,1:length(mlrrun[[1]])]=mlrrun[[1]]
  olrrun=olrfunc(datafit=datatr[[1]],dataval=datapo[[1]],ECI="rcs3",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
  olrsim[i,1:length(olrrun[[1]])]=olrrun[[1]]
  acprun=acpfunc(datafit=datatr[[1]],dataval=datapo[[1]],ECI="rcs3",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
  acpsim[i,1:length(acprun[[1]])]=acprun[[1]]
  smrun=smfunc(datafit=datatr[[1]],dataval=datapo[[1]],ECI="rcs3",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
  smsim[i,1:length(smrun[[1]])]=smrun[[1]]
  print(i)
}
comsave_smalln(filep3="Ordinal prediction model simulations\\ordsimk3p4_equi_unequal_moderateorc_PO_n500.RData")

}


# SCENARIO 3

{
# Get large dataset PO
set.seed(34634) # the seed for the large datasets
getdata_k3p4po <- function(mult){
  x1=rnorm(100*mult,0.4,1.29)
  x2=rnorm(100*mult,0.3,1.16)
  x3=rnorm(100*mult,0.4,1.29)
  x4=rnorm(100*mult,0.3,1.16)
  lp1 = 1.73 - 0.53*x1 - 0.39*x2 - 0.53*x3 - 0.39*x4
  lp2 = 4.15 - 0.53*x1 - 0.39*x2 - 0.53*x3 - 0.39*x4
  ptrue1 = exp(lp1)/(1+exp(lp1))
  ptrue2 = exp(lp2)/(1+exp(lp2)) - ptrue1
  ptrue3 = 1 - ptrue1 - ptrue2
  y=matrix(data=NA,nrow=length(lp1),ncol=1)
  for (i in 1:length(lp1)){
    y[i]=which(grepl(1, rmultinom(1,1,c(ptrue1[i],ptrue2[i],ptrue3[i]))))
  }  
  datax = as.data.frame(cbind(x1,x2,x3,x4,y,ptrue1,ptrue2,ptrue3))
  names(datax)=c("x1","x2","x3","x4","y","ptrue1","ptrue2","ptrue3")
  datax$y = ordered(datax$y)
  return(list(datax,c(mult)))
}
datapo=getdata_k3p4po(mult=2000)

# Run models for large dataset
mlrfit=vglm(y ~ x1 + x2 + x3 + x4, family=multinomial(refLevel = "1"), data=datapo[[1]])
acnpres=mlrfunc(datafit=datapo[[1]],dataval=datapo[[1]],ECI="nonpar",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
olrres=olrfunc(datafit=datapo[[1]],dataval=datapo[[1]],ECI="nonpar",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
acpres=acpfunc(datafit=datapo[[1]],dataval=datapo[[1]],ECI="nonpar",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
smres=smfunc(datafit=datapo[[1]],dataval=datapo[[1]],ECI="nonpar",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
comsave(filep1="Ordinal prediction model simulations\\ordsimk3p4_equi_xunequal_moderateorc_PO.RData")
vgamsmps4a = vgam(datapo[[1]]$y ~ sm.ps(log(acnpres[[3]][,2]/acnpres[[3]][,1]),df=4) + sm.ps(log(acnpres[[3]][,3]/acnpres[[3]][,1]),df=4),family=multinomial(refLevel = 1))
vgamsmps4b = vgam(datapo[[1]]$y ~ sm.ps(log(olrres[[3]][,2]/olrres[[3]][,1]),df=4) + sm.ps(log(olrres[[3]][,3]/olrres[[3]][,1]),df=4),family=multinomial(refLevel = 1))
vgamsmps4c = vgam(datapo[[1]]$y ~ sm.ps(log(acpres[[3]][,2]/acpres[[3]][,1]),df=4) + sm.ps(log(acpres[[3]][,3]/acpres[[3]][,1]),df=4),family=multinomial(refLevel = 1))
vgamsmps4d = vgam(datapo[[1]]$y ~ sm.ps(log(smres[[3]][,2]/smres[[3]][,1]),df=4) + sm.ps(log(smres[[3]][,3]/smres[[3]][,1]),df=4),family=multinomial(refLevel = 1))
plotsave(filep2="Ordinal prediction model simulations\\ordsimk3p4_equi_xunequal_moderateorc_PO.pdf")

# Smooth calibration plots for paper
dev.off()
par(mfrow=c(2,2))
plotlines3sub(preds=acnpres[[3]],obs=vgamsmps4a,mtxt=NULL,subsample0=1000,cexlab=1.4,showleg=0,showtext=1,texttoshow="MLR",cextext=1.3)
plotlines3sub(preds=olrres[[3]],obs=vgamsmps4b,mtxt=NULL,subsample0=1000,cexlab=1.4,showleg=0,showtext=1,texttoshow="CL-PO",cextext=1.3)
plotlines3sub(preds=acpres[[3]],obs=vgamsmps4c,mtxt=NULL,subsample0=1000,cexlab=1.4,showleg=0,showtext=1,texttoshow="AC-PO",cextext=1.3)
plotlines3sub(preds=smres[[3]],obs=vgamsmps4d,mtxt=NULL,subsample0=1000,cexlab=1.4,showleg=0,showtext=1,texttoshow="SLM",cextext=1.3)
# Scatter plot vs true risk for paper
plotvstruerisk3(subsample=1000)

# dev.off()
# olrx <- orm(y ~ x1 + x2 + x3 + x4, data=datapo[[1]], x=T)
# par(mfrow=c(2,2))
# resid(olrx, 'partial',pl=T,label.curves=F) 

# small sample simulations to check overfitting (N=100) - GIVES ERRORS SOMETIMES SO I DO 8 RUNS TO MAKE SURE I GET 200 RESULTS
set.seed(4682)
mlrsim=matrix(data=NA,nrow=200,ncol=27)
olrsim=matrix(data=NA,nrow=200,ncol=27)
acpsim=matrix(data=NA,nrow=200,ncol=27)
smsim=matrix(data=NA,nrow=200,ncol=27)
i=1
while (i<=200) { # use while function to avoid datasets with 0 or 1 case in category 3
  datatr=getdata_k3p4po(mult=1)
  if (sum(table(datatr[[1]]$y)[1:2])<=97){
    mlrrun=mlrfunc(datafit=datatr[[1]],dataval=datapo[[1]],ECI="none",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
    mlrsim[i,1:length(mlrrun[[1]])]=mlrrun[[1]]
    olrrun=olrfunc(datafit=datatr[[1]],dataval=datapo[[1]],ECI="none",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
    olrsim[i,1:length(olrrun[[1]])]=olrrun[[1]]
    acprun=acpfunc(datafit=datatr[[1]],dataval=datapo[[1]],ECI="none",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
    acpsim[i,1:length(acprun[[1]])]=acprun[[1]]
    smrun=smfunc(datafit=datatr[[1]],dataval=datapo[[1]],ECI="none",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
    smsim[i,1:length(smrun[[1]])]=smrun[[1]]
    print(i)
    i=i+1
    comsave_smalln(filep3="Ordinal prediction model simulations\\ordsimk3p4_equi_xunequal_moderateorc_PO_n100a.RData")
  }
}
set.seed(8495)
mlrsim=matrix(data=NA,nrow=200,ncol=27)
olrsim=matrix(data=NA,nrow=200,ncol=27)
acpsim=matrix(data=NA,nrow=200,ncol=27)
smsim=matrix(data=NA,nrow=200,ncol=27)
i=1
while (i<=200) { # use while function to avoid datasets with 0 or 1 case in category 3
  datatr=getdata_k3p4po(mult=1)
  if (sum(table(datatr[[1]]$y)[1:2])<=97){
    mlrrun=mlrfunc(datafit=datatr[[1]],dataval=datapo[[1]],ECI="none",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
    mlrsim[i,1:length(mlrrun[[1]])]=mlrrun[[1]]
    olrrun=olrfunc(datafit=datatr[[1]],dataval=datapo[[1]],ECI="none",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
    olrsim[i,1:length(olrrun[[1]])]=olrrun[[1]]
    acprun=acpfunc(datafit=datatr[[1]],dataval=datapo[[1]],ECI="none",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
    acpsim[i,1:length(acprun[[1]])]=acprun[[1]]
    smrun=smfunc(datafit=datatr[[1]],dataval=datapo[[1]],ECI="none",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
    smsim[i,1:length(smrun[[1]])]=smrun[[1]]
    print(i)
    i=i+1
    comsave_smalln(filep3="Ordinal prediction model simulations\\ordsimk3p4_equi_xunequal_moderateorc_PO_n100b.RData")
  }
}
set.seed(6652)
mlrsim=matrix(data=NA,nrow=200,ncol=27)
olrsim=matrix(data=NA,nrow=200,ncol=27)
acpsim=matrix(data=NA,nrow=200,ncol=27)
smsim=matrix(data=NA,nrow=200,ncol=27)
i=1
while (i<=200) { # use while function to avoid datasets with 0 or 1 case in category 3
  datatr=getdata_k3p4po(mult=1)
  if (sum(table(datatr[[1]]$y)[1:2])<=97){
    mlrrun=mlrfunc(datafit=datatr[[1]],dataval=datapo[[1]],ECI="none",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
    mlrsim[i,1:length(mlrrun[[1]])]=mlrrun[[1]]
    olrrun=olrfunc(datafit=datatr[[1]],dataval=datapo[[1]],ECI="none",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
    olrsim[i,1:length(olrrun[[1]])]=olrrun[[1]]
    acprun=acpfunc(datafit=datatr[[1]],dataval=datapo[[1]],ECI="none",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
    acpsim[i,1:length(acprun[[1]])]=acprun[[1]]
    smrun=smfunc(datafit=datatr[[1]],dataval=datapo[[1]],ECI="none",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
    smsim[i,1:length(smrun[[1]])]=smrun[[1]]
    print(i)
    i=i+1
    comsave_smalln(filep3="Ordinal prediction model simulations\\ordsimk3p4_equi_xunequal_moderateorc_PO_n100c.RData")
  }
}
set.seed(4989)
mlrsim=matrix(data=NA,nrow=200,ncol=27)
olrsim=matrix(data=NA,nrow=200,ncol=27)
acpsim=matrix(data=NA,nrow=200,ncol=27)
smsim=matrix(data=NA,nrow=200,ncol=27)
i=1
while (i<=200) { # use while function to avoid datasets with 0 or 1 case in category 3
  datatr=getdata_k3p4po(mult=1)
  if (sum(table(datatr[[1]]$y)[1:2])<=97){
    mlrrun=mlrfunc(datafit=datatr[[1]],dataval=datapo[[1]],ECI="none",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
    mlrsim[i,1:length(mlrrun[[1]])]=mlrrun[[1]]
    olrrun=olrfunc(datafit=datatr[[1]],dataval=datapo[[1]],ECI="none",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
    olrsim[i,1:length(olrrun[[1]])]=olrrun[[1]]
    acprun=acpfunc(datafit=datatr[[1]],dataval=datapo[[1]],ECI="none",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
    acpsim[i,1:length(acprun[[1]])]=acprun[[1]]
    smrun=smfunc(datafit=datatr[[1]],dataval=datapo[[1]],ECI="none",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
    smsim[i,1:length(smrun[[1]])]=smrun[[1]]
    print(i)
    i=i+1
    comsave_smalln(filep3="Ordinal prediction model simulations\\ordsimk3p4_equi_xunequal_moderateorc_PO_n100d.RData")
  }
}
set.seed(9000)
mlrsim=matrix(data=NA,nrow=200,ncol=27)
olrsim=matrix(data=NA,nrow=200,ncol=27)
acpsim=matrix(data=NA,nrow=200,ncol=27)
smsim=matrix(data=NA,nrow=200,ncol=27)
i=1
while (i<=200) { # use while function to avoid datasets with 0 or 1 case in category 3
  datatr=getdata_k3p4po(mult=1)
  if (sum(table(datatr[[1]]$y)[1:2])<=97){
    mlrrun=mlrfunc(datafit=datatr[[1]],dataval=datapo[[1]],ECI="none",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
    mlrsim[i,1:length(mlrrun[[1]])]=mlrrun[[1]]
    olrrun=olrfunc(datafit=datatr[[1]],dataval=datapo[[1]],ECI="none",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
    olrsim[i,1:length(olrrun[[1]])]=olrrun[[1]]
    acprun=acpfunc(datafit=datatr[[1]],dataval=datapo[[1]],ECI="none",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
    acpsim[i,1:length(acprun[[1]])]=acprun[[1]]
    smrun=smfunc(datafit=datatr[[1]],dataval=datapo[[1]],ECI="none",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
    smsim[i,1:length(smrun[[1]])]=smrun[[1]]
    print(i)
    i=i+1
    comsave_smalln(filep3="Ordinal prediction model simulations\\ordsimk3p4_equi_xunequal_moderateorc_PO_n100e.RData")
  }
}
set.seed(5520)
mlrsim=matrix(data=NA,nrow=200,ncol=27)
olrsim=matrix(data=NA,nrow=200,ncol=27)
acpsim=matrix(data=NA,nrow=200,ncol=27)
smsim=matrix(data=NA,nrow=200,ncol=27)
i=1
while (i<=200) { # use while function to avoid datasets with 0 or 1 case in category 3
  datatr=getdata_k3p4po(mult=1)
  if (sum(table(datatr[[1]]$y)[1:2])<=97){
    mlrrun=mlrfunc(datafit=datatr[[1]],dataval=datapo[[1]],ECI="none",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
    mlrsim[i,1:length(mlrrun[[1]])]=mlrrun[[1]]
    olrrun=olrfunc(datafit=datatr[[1]],dataval=datapo[[1]],ECI="none",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
    olrsim[i,1:length(olrrun[[1]])]=olrrun[[1]]
    acprun=acpfunc(datafit=datatr[[1]],dataval=datapo[[1]],ECI="none",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
    acpsim[i,1:length(acprun[[1]])]=acprun[[1]]
    smrun=smfunc(datafit=datatr[[1]],dataval=datapo[[1]],ECI="none",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
    smsim[i,1:length(smrun[[1]])]=smrun[[1]]
    print(i)
    i=i+1
    comsave_smalln(filep3="Ordinal prediction model simulations\\ordsimk3p4_equi_xunequal_moderateorc_PO_n100f.RData")
  }
}
set.seed(1118)
mlrsim=matrix(data=NA,nrow=200,ncol=27)
olrsim=matrix(data=NA,nrow=200,ncol=27)
acpsim=matrix(data=NA,nrow=200,ncol=27)
smsim=matrix(data=NA,nrow=200,ncol=27)
i=1
while (i<=200) { # use while function to avoid datasets with 0 or 1 case in category 3
  datatr=getdata_k3p4po(mult=1)
  if (sum(table(datatr[[1]]$y)[1:2])<=97){
    mlrrun=mlrfunc(datafit=datatr[[1]],dataval=datapo[[1]],ECI="none",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
    mlrsim[i,1:length(mlrrun[[1]])]=mlrrun[[1]]
    olrrun=olrfunc(datafit=datatr[[1]],dataval=datapo[[1]],ECI="none",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
    olrsim[i,1:length(olrrun[[1]])]=olrrun[[1]]
    acprun=acpfunc(datafit=datatr[[1]],dataval=datapo[[1]],ECI="none",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
    acpsim[i,1:length(acprun[[1]])]=acprun[[1]]
    smrun=smfunc(datafit=datatr[[1]],dataval=datapo[[1]],ECI="none",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
    smsim[i,1:length(smrun[[1]])]=smrun[[1]]
    print(i)
    i=i+1
    comsave_smalln(filep3="Ordinal prediction model simulations\\ordsimk3p4_equi_xunequal_moderateorc_PO_n100g.RData")
  }
}
set.seed(2987)
mlrsim=matrix(data=NA,nrow=200,ncol=27)
olrsim=matrix(data=NA,nrow=200,ncol=27)
acpsim=matrix(data=NA,nrow=200,ncol=27)
smsim=matrix(data=NA,nrow=200,ncol=27)
i=1
while (i<=200) { # use while function to avoid datasets with 0 or 1 case in category 3
  datatr=getdata_k3p4po(mult=1)
  if (sum(table(datatr[[1]]$y)[1:2])<=97){
    mlrrun=mlrfunc(datafit=datatr[[1]],dataval=datapo[[1]],ECI="none",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
    mlrsim[i,1:length(mlrrun[[1]])]=mlrrun[[1]]
    olrrun=olrfunc(datafit=datatr[[1]],dataval=datapo[[1]],ECI="none",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
    olrsim[i,1:length(olrrun[[1]])]=olrrun[[1]]
    acprun=acpfunc(datafit=datatr[[1]],dataval=datapo[[1]],ECI="none",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
    acpsim[i,1:length(acprun[[1]])]=acprun[[1]]
    smrun=smfunc(datafit=datatr[[1]],dataval=datapo[[1]],ECI="none",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
    smsim[i,1:length(smrun[[1]])]=smrun[[1]]
    print(i)
    i=i+1
    comsave_smalln(filep3="Ordinal prediction model simulations\\ordsimk3p4_equi_xunequal_moderateorc_PO_n100h.RData")
  }
}

# small sample simulations to check overfitting (N=500)
set.seed(4682)
mlrsim=matrix(data=NA,nrow=200,ncol=27)
olrsim=matrix(data=NA,nrow=200,ncol=27)
acpsim=matrix(data=NA,nrow=200,ncol=27)
smsim=matrix(data=NA,nrow=200,ncol=27)
for(i in 1:200) {
  datatr=getdata_k3p4po(mult=5)
  mlrrun=mlrfunc(datafit=datatr[[1]],dataval=datapo[[1]],ECI="rcs3",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
  mlrsim[i,1:length(mlrrun[[1]])]=mlrrun[[1]]
  olrrun=olrfunc(datafit=datatr[[1]],dataval=datapo[[1]],ECI="rcs3",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
  olrsim[i,1:length(olrrun[[1]])]=olrrun[[1]]
  acprun=acpfunc(datafit=datatr[[1]],dataval=datapo[[1]],ECI="rcs3",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
  acpsim[i,1:length(acprun[[1]])]=acprun[[1]]
  smrun=smfunc(datafit=datatr[[1]],dataval=datapo[[1]],ECI="rcs3",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
  smsim[i,1:length(smrun[[1]])]=smrun[[1]]
  print(i)
}
comsave_smalln(filep3="Ordinal prediction model simulations\\ordsimk3p4_equi_xunequal_moderateorc_PO_n500.RData")

}


# SCENARIO 4

{
# Get large dataset PO
set.seed(34634) # the seed for the large datasets
getdata_k4p3po <- function(mult){
  x1=rnorm(100*mult,0.3,1.50)
  x2=rnorm(100*mult,0.2,1.29)
  x3=rnorm(100*mult,0.3,1.50)
  lp1 = -0.12 - 0.54*x1 - 0.47*x2 - 0.51*x3
  lp2 = 1.22 - 0.54*x1 - 0.47*x2 - 0.51*x3
  lp3 = 2.62 - 0.54*x1 - 0.47*x2 - 0.51*x3
  ptrue1 = exp(lp1)/(1+exp(lp1))
  ptrue2 = exp(lp2)/(1+exp(lp2)) - ptrue1
  ptrue3 = exp(lp3)/(1+exp(lp3)) - ptrue1 - ptrue2
  ptrue4 = 1 - ptrue1 - ptrue2 - ptrue3
  y=matrix(data=NA,nrow=length(lp1),ncol=1)
  for (i in 1:length(lp1)){
    y[i]=which(grepl(1, rmultinom(1,1,c(ptrue1[i],ptrue2[i],ptrue3[i],ptrue4[i]))))
  }  
  datax = as.data.frame(cbind(x1,x2,x3,y,ptrue1,ptrue2,ptrue3,ptrue4))
  names(datax)=c("x1","x2","x3","y","ptrue1","ptrue2","ptrue3","ptrue4")
  datax$y = ordered(datax$y)
  return(list(datax,c(mult)))
}
datapo=getdata_k4p3po(mult=2000)

# Run models for large dataset
mlrfit=vglm(y ~ x1 + x2 + x3, family=multinomial(refLevel = "1"), data=datapo[[1]])
acnpres=mlrfunc(datafit=datapo[[1]],dataval=datapo[[1]],ECI="nonpar",levels=4,modelformula=y ~ x1 + x2 + x3)
olrres=olrfunc(datafit=datapo[[1]],dataval=datapo[[1]],ECI="nonpar",levels=4,modelformula=y ~ x1 + x2 + x3)
acpres=acpfunc(datafit=datapo[[1]],dataval=datapo[[1]],ECI="nonpar",levels=4,modelformula=y ~ x1 + x2 + x3)
smres=smfunc(datafit=datapo[[1]],dataval=datapo[[1]],ECI="nonpar",levels=4,modelformula=y ~ x1 + x2 + x3)
comsave(filep1="Ordinal prediction model simulations\\ordsimk4p3_equi_unequal_moderateorc_PO.RData")
vgamsmps4a = vgam(datapo[[1]]$y ~ sm.ps(log(acnpres[[3]][,2]/acnpres[[3]][,1]),df=4) + sm.ps(log(acnpres[[3]][,3]/acnpres[[3]][,1]),df=4) + sm.ps(log(acnpres[[3]][,4]/acnpres[[3]][,1]),df=4),family=multinomial(refLevel = 1))
vgamsmps4b = vgam(datapo[[1]]$y ~ sm.ps(log(olrres[[3]][,2]/olrres[[3]][,1]),df=4) + sm.ps(log(olrres[[3]][,3]/olrres[[3]][,1]),df=4) + sm.ps(log(olrres[[3]][,4]/olrres[[3]][,1]),df=4),family=multinomial(refLevel = 1))
vgamsmps4c = vgam(datapo[[1]]$y ~ sm.ps(log(acpres[[3]][,2]/acpres[[3]][,1]),df=4) + sm.ps(log(acpres[[3]][,3]/acpres[[3]][,1]),df=4) + sm.ps(log(acpres[[3]][,4]/acpres[[3]][,1]),df=4),family=multinomial(refLevel = 1))
vgamsmps4d = vgam(datapo[[1]]$y ~ sm.ps(log(smres[[3]][,2]/smres[[3]][,1]),df=4) + sm.ps(log(smres[[3]][,3]/smres[[3]][,1]),df=4) + sm.ps(log(smres[[3]][,4]/smres[[3]][,1]),df=4),family=multinomial(refLevel = 1))
plotsave(filep2="Ordinal prediction model simulations\\ordsimk4p3_equi_unequal_moderateorc_PO.pdf")

# Smooth calibration plots for paper
dev.off()
par(mfrow=c(2,2))
plotlines4sub(preds=acnpres[[3]],obs=vgamsmps4a,mtxt=NULL,subsample0=1000,cexlab=1.4,showleg=0,showtext=1,texttoshow="MLR",cextext=1.3)
plotlines4sub(preds=olrres[[3]],obs=vgamsmps4b,mtxt=NULL,subsample0=1000,cexlab=1.4,showleg=0,showtext=1,texttoshow="CL-PO",cextext=1.3)
plotlines4sub(preds=acpres[[3]],obs=vgamsmps4c,mtxt=NULL,subsample0=1000,spar2=1.2,cexlab=1.4,showleg=0,showtext=1,texttoshow="AC-PO",cextext=1.3)
plotlines4sub(preds=smres[[3]],obs=vgamsmps4d,mtxt=NULL,subsample0=1000,cexlab=1.4,showleg=0,showtext=1,texttoshow="SLM",cextext=1.3)
# Scatter plot vs true risk for paper
plotvstruerisk4(subsample=1000)

# dev.off()
# olrx <- orm(y ~ x1 + x2 + x3, data=datapo[[1]], x=T)
# par(mfrow=c(2,2))
# resid(olrx, 'partial',pl=T,label.curves=F) 

}


# SCENARIO 5

{
# Get large dataset PO
set.seed(34634) # the seed for the large datasets
getdata_k4p3po <- function(mult){
  x1=rnorm(100*mult,0.3,0.90)
  x2=rnorm(100*mult,0.2,0.80)
  x3=rnorm(100*mult,0.3,0.90)
  lp1 = -0.05 - 0.54*x1 - 0.47*x2 - 0.51*x3
  lp2 = 1.1 - 0.54*x1 - 0.47*x2 - 0.51*x3
  lp3 = 2.35 - 0.54*x1 - 0.47*x2 - 0.51*x3
  ptrue1 = exp(lp1)/(1+exp(lp1))
  ptrue2 = exp(lp2)/(1+exp(lp2)) - ptrue1
  ptrue3 = exp(lp3)/(1+exp(lp3)) - ptrue1 - ptrue2
  ptrue4 = 1 - ptrue1 - ptrue2 - ptrue3
  y=matrix(data=NA,nrow=length(lp1),ncol=1)
  for (i in 1:length(lp1)){
    y[i]=which(grepl(1, rmultinom(1,1,c(ptrue1[i],ptrue2[i],ptrue3[i],ptrue4[i]))))
  }  
  datax = as.data.frame(cbind(x1,x2,x3,y,ptrue1,ptrue2,ptrue3,ptrue4))
  names(datax)=c("x1","x2","x3","y","ptrue1","ptrue2","ptrue3","ptrue4")
  datax$y = ordered(datax$y)
  return(list(datax,c(mult)))
}
datapo=getdata_k4p3po(mult=2000)

# Run models for large dataset
mlrfit=vglm(y ~ x1 + x2 + x3, family=multinomial(refLevel = "1"), data=datapo[[1]])
acnpres=mlrfunc(datafit=datapo[[1]],dataval=datapo[[1]],ECI="nonpar",levels=4,modelformula=y ~ x1 + x2 + x3)
olrres=olrfunc(datafit=datapo[[1]],dataval=datapo[[1]],ECI="nonpar",levels=4,modelformula=y ~ x1 + x2 + x3)
acpres=acpfunc(datafit=datapo[[1]],dataval=datapo[[1]],ECI="nonpar",levels=4,modelformula=y ~ x1 + x2 + x3)
smres=smfunc(datafit=datapo[[1]],dataval=datapo[[1]],ECI="nonpar",levels=4,modelformula=y ~ x1 + x2 + x3)
comsave(filep1="Ordinal prediction model simulations\\ordsimk4p3_equi_unequal_loworc_PO.RData")
vgamsmps4a = vgam(datapo[[1]]$y ~ sm.ps(log(acnpres[[3]][,2]/acnpres[[3]][,1]),df=4) + sm.ps(log(acnpres[[3]][,3]/acnpres[[3]][,1]),df=4) + sm.ps(log(acnpres[[3]][,4]/acnpres[[3]][,1]),df=4),family=multinomial(refLevel = 1))
vgamsmps4b = vgam(datapo[[1]]$y ~ sm.ps(log(olrres[[3]][,2]/olrres[[3]][,1]),df=4) + sm.ps(log(olrres[[3]][,3]/olrres[[3]][,1]),df=4) + sm.ps(log(olrres[[3]][,4]/olrres[[3]][,1]),df=4),family=multinomial(refLevel = 1))
vgamsmps4c = vgam(datapo[[1]]$y ~ sm.ps(log(acpres[[3]][,2]/acpres[[3]][,1]),df=4) + sm.ps(log(acpres[[3]][,3]/acpres[[3]][,1]),df=4) + sm.ps(log(acpres[[3]][,4]/acpres[[3]][,1]),df=4),family=multinomial(refLevel = 1))
vgamsmps4d = vgam(datapo[[1]]$y ~ sm.ps(log(smres[[3]][,2]/smres[[3]][,1]),df=4) + sm.ps(log(smres[[3]][,3]/smres[[3]][,1]),df=4) + sm.ps(log(smres[[3]][,4]/smres[[3]][,1]),df=4),family=multinomial(refLevel = 1))
plotsave(filep2="Ordinal prediction model simulations\\ordsimk4p3_equi_unequal_loworc_PO.pdf")

# Smooth calibration plots for paper
dev.off()
par(mfrow=c(2,2))
plotlines4sub(preds=acnpres[[3]],obs=vgamsmps4a,mtxt=NULL,subsample0=1000,cexlab=1.4,showleg=0,showtext=1,texttoshow="MLR",cextext=1.3)
plotlines4sub(preds=olrres[[3]],obs=vgamsmps4b,mtxt=NULL,subsample0=1000,cexlab=1.4,showleg=0,showtext=1,texttoshow="CL-PO",cextext=1.3)
plotlines4sub(preds=acpres[[3]],obs=vgamsmps4c,mtxt=NULL,subsample0=1000,spar2=1.2,cexlab=1.4,showleg=0,showtext=1,texttoshow="AC-PO",cextext=1.3)
plotlines4sub(preds=smres[[3]],obs=vgamsmps4d,mtxt=NULL,subsample0=1000,cexlab=1.4,showleg=0,showtext=1,texttoshow="SLM",cextext=1.3)
# Scatter plot vs true risk for paper
plotvstruerisk4(subsample=1000)

# dev.off()
# olrx <- orm(y ~ x1 + x2 + x3, data=datapo[[1]], x=T)
# par(mfrow=c(2,2))
# resid(olrx, 'partial',pl=T,label.curves=F) 

}


# SCENARIO 6

{
# Get large dataset PO
set.seed(34634) # the seed for the large datasets
getdata_k4p3po <- function(mult){
  x1=rnorm(100*mult,0.3,0.85)
  x2=rnorm(100*mult,0.2,0.75)
  x3=rnorm(100*mult,0.3,0.85)
  lp1 = 0.18 - 0.54*x1 - 0.47*x2 - 0.51*x3
  lp2 = 1.63 - 0.54*x1 - 0.47*x2 - 0.51*x3
  lp3 = 3.59 - 0.54*x1 - 0.47*x2 - 0.51*x3
  ptrue1 = exp(lp1)/(1+exp(lp1))
  ptrue2 = exp(lp2)/(1+exp(lp2)) - ptrue1
  ptrue3 = exp(lp3)/(1+exp(lp3)) - ptrue1 - ptrue2
  ptrue4 = 1 - ptrue1 - ptrue2 - ptrue3
  y=matrix(data=NA,nrow=length(lp1),ncol=1)
  for (i in 1:length(lp1)){
    y[i]=which(grepl(1, rmultinom(1,1,c(ptrue1[i],ptrue2[i],ptrue3[i],ptrue4[i]))))
  }  
  datax = as.data.frame(cbind(x1,x2,x3,y,ptrue1,ptrue2,ptrue3,ptrue4))
  names(datax)=c("x1","x2","x3","y","ptrue1","ptrue2","ptrue3","ptrue4")
  datax$y = ordered(datax$y)
  return(list(datax,c(mult)))
}
datapo=getdata_k4p3po(mult=2000)

# Run models for large dataset
mlrfit=vglm(y ~ x1 + x2 + x3, family=multinomial(refLevel = "1"), data=datapo[[1]])
acnpres=mlrfunc(datafit=datapo[[1]],dataval=datapo[[1]],ECI="nonpar",levels=4,modelformula=y ~ x1 + x2 + x3)
olrres=olrfunc(datafit=datapo[[1]],dataval=datapo[[1]],ECI="nonpar",levels=4,modelformula=y ~ x1 + x2 + x3)
acpres=acpfunc(datafit=datapo[[1]],dataval=datapo[[1]],ECI="nonpar",levels=4,modelformula=y ~ x1 + x2 + x3)
smres=smfunc(datafit=datapo[[1]],dataval=datapo[[1]],ECI="nonpar",levels=4,modelformula=y ~ x1 + x2 + x3)
comsave(filep1="Ordinal prediction model simulations\\ordsimk4p3_equi_xunequal_loworc_PO.RData")
vgamsmps4a = vgam(datapo[[1]]$y ~ sm.ps(log(acnpres[[3]][,2]/acnpres[[3]][,1]),df=4) + sm.ps(log(acnpres[[3]][,3]/acnpres[[3]][,1]),df=4) + sm.ps(log(acnpres[[3]][,4]/acnpres[[3]][,1]),df=4),family=multinomial(refLevel = 1))
vgamsmps4b = vgam(datapo[[1]]$y ~ sm.ps(log(olrres[[3]][,2]/olrres[[3]][,1]),df=4) + sm.ps(log(olrres[[3]][,3]/olrres[[3]][,1]),df=4) + sm.ps(log(olrres[[3]][,4]/olrres[[3]][,1]),df=4),family=multinomial(refLevel = 1))
vgamsmps4c = vgam(datapo[[1]]$y ~ sm.ps(log(acpres[[3]][,2]/acpres[[3]][,1]),df=4) + sm.ps(log(acpres[[3]][,3]/acpres[[3]][,1]),df=4) + sm.ps(log(acpres[[3]][,4]/acpres[[3]][,1]),df=4),family=multinomial(refLevel = 1))
vgamsmps4d = vgam(datapo[[1]]$y ~ sm.ps(log(smres[[3]][,2]/smres[[3]][,1]),df=4) + sm.ps(log(smres[[3]][,3]/smres[[3]][,1]),df=4) + sm.ps(log(smres[[3]][,4]/smres[[3]][,1]),df=4),family=multinomial(refLevel = 1))
plotsave(filep2="Ordinal prediction model simulations\\ordsimk4p3_equi_xunequal_loworc_PO.pdf")

# Smooth calibration plots for paper
dev.off()
par(mfrow=c(2,2))
plotlines4sub(preds=acnpres[[3]],obs=vgamsmps4a,mtxt=NULL,subsample0=1000,cexlab=1.4,showleg=0,showtext=1,texttoshow="MLR",cextext=1.3)
plotlines4sub(preds=olrres[[3]],obs=vgamsmps4b,mtxt=NULL,subsample0=1000,cexlab=1.4,showleg=0,showtext=1,texttoshow="CL-PO",cextext=1.3)
plotlines4sub(preds=acpres[[3]],obs=vgamsmps4c,mtxt=NULL,subsample0=1000,cexlab=1.4,showleg=0,showtext=1,texttoshow="AC-PO",cextext=1.3)
plotlines4sub(preds=smres[[3]],obs=vgamsmps4d,mtxt=NULL,subsample0=1000,cexlab=1.4,showleg=0,showtext=1,texttoshow="SLM",cextext=1.3)
# Scatter plot vs true risk for paper
plotvstruerisk4(subsample=1000)

# dev.off()
# olrx <- orm(y ~ x1 + x2 + x3, data=datapo[[1]], x=T)
# par(mfrow=c(2,2))
# resid(olrx, 'partial',pl=T,label.curves=F) 

}


# SCENARIO 7

{
# Get large dataset PO
set.seed(34634) # the seed for the large datasets
getdata_k3p4pob <- function(mult,px1,px2,px3,px4){
  x1=sample(0:1, 100*mult, replace=TRUE,prob=c(1-px1,px1))
  x2=sample(0:1, 100*mult, replace=TRUE,prob=c(1-px2,px2))
  x3=sample(0:1, 100*mult, replace=TRUE,prob=c(1-px3,px3))
  x4=sample(0:1, 100*mult, replace=TRUE,prob=c(1-px4,px4))
  lp1 = 2.04 - 1.38*x1 - 1.09*x2 - 1.22*x3 - 0.84*x4
  lp2 = 3.95 - 1.38*x1 - 1.09*x2 - 1.22*x3 - 0.84*x4
  ptrue1 = exp(lp1)/(1+exp(lp1))
  ptrue2 = exp(lp2)/(1+exp(lp2)) - ptrue1
  ptrue3 = 1 - ptrue1 - ptrue2
  y=matrix(data=NA,nrow=length(lp1),ncol=1)
  for (i in 1:length(lp1)){
    y[i]=which(grepl(1, rmultinom(1,1,c(ptrue1[i],ptrue2[i],ptrue3[i]))))
  }  
  datax = as.data.frame(cbind(x1,x2,x3,x4,y,ptrue1,ptrue2,ptrue3))
  names(datax)=c("x1","x2","x3","x4","y","ptrue1","ptrue2","ptrue3")
  datax$y = ordered(datax$y)
  return(list(datax,c(mult)))
}
datapo=getdata_k3p4pob(mult=2000,px1=0.45,px2=0.4,px3=0.4,px4=0.3)

# Run models for large dataset
mlrfit=vglm(y ~ x1 + x2 + x3 + x4, family=multinomial(refLevel = "1"), data=datapo[[1]])
acnpres=mlrfunc(datafit=datapo[[1]],dataval=datapo[[1]],ECI="nonpar",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
olrres=olrfunc(datafit=datapo[[1]],dataval=datapo[[1]],ECI="nonpar",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
acpres=acpfunc(datafit=datapo[[1]],dataval=datapo[[1]],ECI="nonpar",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
smres=smfunc(datafit=datapo[[1]],dataval=datapo[[1]],ECI="nonpar",levels=3,modelformula=y ~ x1 + x2 + x3 + x4)
comsave(filep1="Ordinal prediction model simulations\\ordsimk3p4bin_equi_unequal_moderateorc_PO.RData")
vgamsmps4a = vgam(datapo[[1]]$y ~ sm.ps(log(acnpres[[3]][,2]/acnpres[[3]][,1]),df=4) + sm.ps(log(acnpres[[3]][,3]/acnpres[[3]][,1]),df=4),family=multinomial(refLevel = 1))
vgamsmps4b = vgam(datapo[[1]]$y ~ sm.ps(log(olrres[[3]][,2]/olrres[[3]][,1]),df=4) + sm.ps(log(olrres[[3]][,3]/olrres[[3]][,1]),df=4),family=multinomial(refLevel = 1))
vgamsmps4c = vgam(datapo[[1]]$y ~ sm.ps(log(acpres[[3]][,2]/acpres[[3]][,1]),df=4) + sm.ps(log(acpres[[3]][,3]/acpres[[3]][,1]),df=4),family=multinomial(refLevel = 1))
vgamsmps4d = vgam(datapo[[1]]$y ~ sm.ps(log(smres[[3]][,2]/smres[[3]][,1]),df=4) + sm.ps(log(smres[[3]][,3]/smres[[3]][,1]),df=4),family=multinomial(refLevel = 1))
plotsave(filep2="Ordinal prediction model simulations\\ordsimk3p4bin_equi_unequal_moderateorc_PO.pdf")

# Smooth calibration plots for paper
dev.off()
par(mfrow=c(2,2))
plotscatter3sub(preds=acnpres[[3]],obs=vgamsmps4a,mtxt=NULL,subsample0=1000,cexlab=1.4,showleg=0,showtext=1,texttoshow="MLR",cextext=1.3)
plotscatter3sub(preds=olrres[[3]],obs=vgamsmps4b,mtxt=NULL,subsample0=1000,cexlab=1.4,showleg=0,showtext=1,texttoshow="CL-PO",cextext=1.3)
plotscatter3sub(preds=acpres[[3]],obs=vgamsmps4c,mtxt=NULL,subsample0=1000,cexlab=1.4,showleg=0,showtext=1,texttoshow="AC-PO",cextext=1.3)
plotscatter3sub(preds=smres[[3]],obs=vgamsmps4d,mtxt=NULL,subsample0=1000,cexlab=1.4,showleg=0,showtext=1,texttoshow="SLM",cextext=1.3)
# Scatter plot vs true risk for paper
plotvstruerisk3(subsample=1000,cexpts=0.8)

# dev.off()
# olrx <- orm(y ~ x1 + x2 + x3 + x4, data=datapo[[1]], x=T)
# par(mfrow=c(2,2))
# resid(olrx, 'partial',pl=T,label.curves=F) 

}


# SCENARIO 8

{
# Get large dataset PO
set.seed(34634) # the seed for the large datasets
getdata_k4p3pob <- function(mult,px1,px2,px3){
  x1=sample(0:1, 100*mult, replace=TRUE,prob=c(1-px1,px1))
  x2=sample(0:1, 100*mult, replace=TRUE,prob=c(1-px2,px2))
  x3=sample(0:1, 100*mult, replace=TRUE,prob=c(1-px3,px3))
  lp1 = 1.52 - 1.77*x1 - 1.2*x2 - 1.46*x3
  lp2 = 2.89 - 1.77*x1 - 1.2*x2 - 1.46*x3
  lp3 = 4.31 - 1.77*x1 - 1.2*x2 - 1.46*x3
  ptrue1 = exp(lp1)/(1+exp(lp1))
  ptrue2 = exp(lp2)/(1+exp(lp2)) - ptrue1
  ptrue3 = exp(lp3)/(1+exp(lp3)) - ptrue1 - ptrue2
  ptrue4 = 1 - ptrue1 - ptrue2 - ptrue3
  y=matrix(data=NA,nrow=length(lp1),ncol=1)
  for (i in 1:length(lp1)){
    y[i]=which(grepl(1, rmultinom(1,1,c(ptrue1[i],ptrue2[i],ptrue3[i],ptrue4[i]))))
  }  
  datax = as.data.frame(cbind(x1,x2,x3,y,ptrue1,ptrue2,ptrue3,ptrue4))
  names(datax)=c("x1","x2","x3","y","ptrue1","ptrue2","ptrue3","ptrue4")
  datax$y = ordered(datax$y)
  return(list(datax,c(mult)))
}
datapo=getdata_k4p3pob(mult=2000,px1=0.4,px2=0.4,px3=0.6)

# Run models for large dataset
mlrfit=vglm(y ~ x1 + x2 + x3, family=multinomial(refLevel = "1"), data=datapo[[1]])
acnpres=mlrfunc(datafit=datapo[[1]],dataval=datapo[[1]],ECI="nonpar",levels=4,modelformula=y ~ x1 + x2 + x3)
olrres=olrfunc(datafit=datapo[[1]],dataval=datapo[[1]],ECI="nonpar",levels=4,modelformula=y ~ x1 + x2 + x3)
acpres=acpfunc(datafit=datapo[[1]],dataval=datapo[[1]],ECI="nonpar",levels=4,modelformula=y ~ x1 + x2 + x3)
smres=smfunc(datafit=datapo[[1]],dataval=datapo[[1]],ECI="nonpar",levels=4,modelformula=y ~ x1 + x2 + x3)
comsave(filep1="Ordinal prediction model simulations\\ordsimk4p3bin_equi_unequal_moderateorc_PO.RData")
vgamsmps4a = vgam(datapo[[1]]$y ~ sm.ps(log(acnpres[[3]][,2]/acnpres[[3]][,1]),df=4) + sm.ps(log(acnpres[[3]][,3]/acnpres[[3]][,1]),df=4),family=multinomial(refLevel = 1))
vgamsmps4b = vgam(datapo[[1]]$y ~ sm.ps(log(olrres[[3]][,2]/olrres[[3]][,1]),df=4) + sm.ps(log(olrres[[3]][,3]/olrres[[3]][,1]),df=4),family=multinomial(refLevel = 1))
vgamsmps4c = vgam(datapo[[1]]$y ~ sm.ps(log(acpres[[3]][,2]/acpres[[3]][,1]),df=4) + sm.ps(log(acpres[[3]][,3]/acpres[[3]][,1]),df=4),family=multinomial(refLevel = 1))
vgamsmps4d = vgam(datapo[[1]]$y ~ sm.ps(log(smres[[3]][,2]/smres[[3]][,1]),df=4) + sm.ps(log(smres[[3]][,3]/smres[[3]][,1]),df=4),family=multinomial(refLevel = 1))
plotsave(filep2="Ordinal prediction model simulations\\ordsimk4p3bin_equi_unequal_moderateorc_PO.pdf")

# Smooth calibration plots for paper
dev.off()
par(mfrow=c(2,2))
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4a,mtxt=NULL,subsample0=1000,cexlab=1.4,showleg=0,showtext=1,texttoshow="MLR",cextext=1.3)
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4b,mtxt=NULL,subsample0=1000,cexlab=1.4,showleg=0,showtext=1,texttoshow="CL-PO",cextext=1.3)
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4c,mtxt=NULL,subsample0=1000,cexlab=1.4,showleg=0,showtext=1,texttoshow="AC-PO",cextext=1.3)
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4d,mtxt=NULL,subsample0=1000,cexlab=1.4,showleg=0,showtext=1,texttoshow="SLM",cextext=1.3)
# Scatter plot vs true risk for paper
plotvstruerisk4(subsample=1000,cexpts=0.8)

# dev.off()
# olrx <- orm(y ~ x1 + x2 + x3 + x4, data=datapo[[1]], x=T)
# par(mfrow=c(2,2))
# resid(olrx, 'partial',pl=T,label.curves=F) 

}


# LOAD RESULTS OF THE 8 SCENARIOS AND SUMMARIZE THEM

load("Ordinal prediction model simulations\\CLPO truth Scenario 1\\ordsimk3p4_equi_equal_moderateorc_PO.rdata")
load("Ordinal prediction model simulations\\CLPO truth Scenario 2\\ordsimk3p4_equi_unequal_moderateorc_PO.rdata")
load("Ordinal prediction model simulations\\CLPO truth Scenario 3\\ordsimk3p4_equi_xunequal_moderateorc_PO.rdata")
load("Ordinal prediction model simulations\\CLPO truth Scenario 4\\ordsimk4p3_equi_unequal_moderateorc_PO.rdata")
load("Ordinal prediction model simulations\\CLPO truth Scenario 5\\ordsimk4p3_equi_unequal_loworc_PO.rdata")
load("Ordinal prediction model simulations\\CLPO truth Scenario 6\\ordsimk4p3_equi_xunequal_loworc_PO.rdata")
load("Ordinal prediction model simulations\\CLPO truth Scenario 7\\ordsimk3p4bin_equi_unequal_moderateorc_PO.rdata")
load("Ordinal prediction model simulations\\CLPO truth Scenario 8\\ordsimk4p3bin_equi_unequal_moderateorc_PO.rdata")

load("Ordinal prediction model simulations\\CLPO truth Scenario 1\\ordsimk3p4_equi_equal_moderateorc_PO_n100.rdata")
load("Ordinal prediction model simulations\\CLPO truth Scenario 2\\ordsimk3p4_equi_unequal_moderateorc_PO_n100.rdata")

# SMALL SAMPLE SIMULATIONS (N=100) FOR SCENARIO 3 RESULTED IN SOME ERRORS
# SPECIFICALLY: FOR RUNS WITH AN ERROR WERE REPLACED TO HAVE 200 RUNS WITHOUT ERRORS
# RUNS ARE RELATED TO SEPARATION ISSUES
a=load("Ordinal prediction model simulations\\CLPO truth Scenario 3\\ordsimk3p4_equi_xunequal_moderateorc_PO_n100a.rdata")
mlrsima=mlrsim[1:28,]
olrsima=olrsim[1:28,]
acpsima=acpsim[1:28,]
smsima=smsim[1:28,]
b=load("Ordinal prediction model simulations\\CLPO truth Scenario 3\\ordsimk3p4_equi_xunequal_moderateorc_PO_n100b.rdata")
mlrsimb=mlrsim[1:37,]
olrsimb=olrsim[1:37,]
acpsimb=acpsim[1:37,]
smsimb=smsim[1:37,]
c=load("Ordinal prediction model simulations\\CLPO truth Scenario 3\\ordsimk3p4_equi_xunequal_moderateorc_PO_n100c.rdata")
mlrsimc=mlrsim[1:46,]
olrsimc=olrsim[1:46,]
acpsimc=acpsim[1:46,]
smsimc=smsim[1:46,]
d=load("Ordinal prediction model simulations\\CLPO truth Scenario 3\\ordsimk3p4_equi_xunequal_moderateorc_PO_n100d.rdata")
mlrsimd=mlrsim[1:69,]
olrsimd=olrsim[1:69,]
acpsimd=acpsim[1:69,]
smsimd=smsim[1:69,]
e=load("Ordinal prediction model simulations\\CLPO truth Scenario 3\\ordsimk3p4_equi_xunequal_moderateorc_PO_n100e.rdata")
mlrsime=mlrsim[1:20,] # rest not needed
olrsime=olrsim[1:20,]
acpsime=acpsim[1:20,]
smsime=smsim[1:20,]
mlrsim=rbind(mlrsima,mlrsimb,mlrsimc,mlrsimd,mlrsime)
olrsim=rbind(olrsima,olrsimb,olrsimc,olrsimd,olrsime)
acpsim=rbind(acpsima,acpsimb,acpsimc,acpsimd,acpsime)
smsim=rbind(smsima,smsimb,smsimc,smsimd,smsime)
res=cbind(colMeans(mlrsim[,1:19]),colMeans(olrsim[,1:19]),colMeans(acpsim[,1:19]),colMeans(smsim[,1:19]),
          colMedians(mlrsim[,1:19]),colMedians(olrsim[,1:19]),colMedians(acpsim[,1:19]),colMedians(smsim[,1:19]))
rownames(res) <- c("ci_cat1","ci_cat2","ci_cat3","cs_cat1","cs_cat2","cs_cat3","ci_gt1","ci_gt2","cs_gt1","cs_gt2","ci_lp1","ci_lp2","cs_lp1","cs_lp2","ci_E","cs_E","rMSPE","Brier","ORC")
colnames(res) <- c("ACNPmean","OLRmean","ACPmean","SMLRmean","ACNPmedian","OLRmedian","ACPmedian","SMLRmedian")

load("Ordinal prediction model simulations\\CLPO truth Scenario 1\\ordsimk3p4_equi_equal_moderateorc_PO_n500.rdata")
load("Ordinal prediction model simulations\\CLPO truth Scenario 2\\ordsimk3p4_equi_unequal_moderateorc_PO_n500.rdata")
load("Ordinal prediction model simulations\\CLPO truth Scenario 3\\ordsimk3p4_equi_xunequal_moderateorc_PO_n500.rdata")