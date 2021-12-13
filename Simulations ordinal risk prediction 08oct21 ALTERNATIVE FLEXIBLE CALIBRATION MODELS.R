####################################################
####################################################
# ORDINAL RISK PREDICTION MODELS: SIMULATION STUDY #
####################################################
####################################################

###########################################
# COMPARING FLEXIBLE RECALIBRATION MODELS #
###########################################

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
source("R code\\ordcalfunctions.R")
source("R code\\ordcalfunctions_simulationstudy.R")


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

flexrecalk3(filenm = "MLR SCENARIO 1 ALTERNATIVE FLEXIBLE CALIBRATION MODELS.RData")

load("MLR SCENARIO 1 ALTERNATIVE FLEXIBLE CALIBRATION MODELS.RData")

round(cbind(ecifrm[1:4],ecifrm[13:16],ecifrm[5:8],ecifrm[17:20],ecifrm[9:12],ecifrm[21:24]), digits = 4)

par(mfrow=c(3,2))
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4azz,mtxt="Ref - MLR")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4azzz,mtxt="Ref - CRNP")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4axx,mtxt="Dich - MLR")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4axxx,mtxt="Dich - CRNP")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4ayy,mtxt="Cat - MLR")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4ayyy,mtxt="Cat - CRNP")
dev.off()
par(mfrow=c(3,2))
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4bzz,mtxt="Ref - MLR")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4bzzz,mtxt="Ref - CRNP")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4bxx,mtxt="Dich - MLR")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4bxxx,mtxt="Dich - CRNP")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4byy,mtxt="Cat - MLR")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4byyy,mtxt="Cat - CRNP")
dev.off()
par(mfrow=c(3,2))
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4czz,mtxt="Ref - MLR")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4czzz,mtxt="Ref - CRNP")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4cxx,mtxt="Dich - MLR")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4cxxx,mtxt="Dich - CRNP")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4cyy,mtxt="Cat - MLR")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4cyyy,mtxt="Cat - CRNP")
dev.off()
par(mfrow=c(3,2))
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dzz,mtxt="Ref - MLR")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dzzz,mtxt="Ref - CRNP")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dxx,mtxt="Dich - MLR")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dxxx,mtxt="Dich - CRNP")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dyy,mtxt="Cat - MLR")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dyyy,mtxt="Cat - CRNP")
dev.off()

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

flexrecalk3(filenm = "MLR SCENARIO 2 ALTERNATIVE FLEXIBLE CALIBRATION MODELS.RData")

load("MLR SCENARIO 2 ALTERNATIVE FLEXIBLE CALIBRATION MODELS.RData")

round(cbind(ecifrm[1:4],ecifrm[13:16],ecifrm[5:8],ecifrm[17:20],ecifrm[9:12],ecifrm[21:24]), digits = 4)

par(mfrow=c(3,2))
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4azz,mtxt="Ref - MLR")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4azzz,mtxt="Ref - CRNP")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4axx,mtxt="Dich - MLR")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4axxx,mtxt="Dich - CRNP")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4ayy,mtxt="Cat - MLR")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4ayyy,mtxt="Cat - CRNP")
dev.off()
par(mfrow=c(3,2))
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4bzz,mtxt="Ref - MLR")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4bzzz,mtxt="Ref - CRNP")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4bxx,mtxt="Dich - MLR")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4bxxx,mtxt="Dich - CRNP")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4byy,mtxt="Cat - MLR")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4byyy,mtxt="Cat - CRNP")
dev.off()
par(mfrow=c(3,2))
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4czz,mtxt="Ref - MLR")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4czzz,mtxt="Ref - CRNP")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4cxx,mtxt="Dich - MLR")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4cxxx,mtxt="Dich - CRNP")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4cyy,mtxt="Cat - MLR")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4cyyy,mtxt="Cat - CRNP")
dev.off()
par(mfrow=c(3,2))
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dzz,mtxt="Ref - MLR")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dzzz,mtxt="Ref - CRNP")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dxx,mtxt="Dich - MLR")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dxxx,mtxt="Dich - CRNP")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dyy,mtxt="Cat - MLR")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dyyy,mtxt="Cat - CRNP")
dev.off()

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

flexrecalk3(filenm = "MLR SCENARIO 3 ALTERNATIVE FLEXIBLE CALIBRATION MODELS.RData")

load("MLR SCENARIO 3 ALTERNATIVE FLEXIBLE CALIBRATION MODELS.RData")

round(cbind(ecifrm[1:4],ecifrm[13:16],ecifrm[5:8],ecifrm[17:20],ecifrm[9:12],ecifrm[21:24]), digits = 4)

par(mfrow=c(3,2))
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4azz,mtxt="Ref - MLR")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4azzz,mtxt="Ref - CRNP")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4axx,mtxt="Dich - MLR")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4axxx,mtxt="Dich - CRNP")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4ayy,mtxt="Cat - MLR")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4ayyy,mtxt="Cat - CRNP")
dev.off()
par(mfrow=c(3,2))
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4bzz,mtxt="Ref - MLR")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4bzzz,mtxt="Ref - CRNP")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4bxx,mtxt="Dich - MLR")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4bxxx,mtxt="Dich - CRNP")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4byy,mtxt="Cat - MLR")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4byyy,mtxt="Cat - CRNP")
dev.off()
par(mfrow=c(3,2))
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4czz,mtxt="Ref - MLR")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4czzz,mtxt="Ref - CRNP")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4cxx,mtxt="Dich - MLR")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4cxxx,mtxt="Dich - CRNP")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4cyy,mtxt="Cat - MLR")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4cyyy,mtxt="Cat - CRNP")
dev.off()
par(mfrow=c(3,2))
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dzz,mtxt="Ref - MLR")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dzzz,mtxt="Ref - CRNP")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dxx,mtxt="Dich - MLR")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dxxx,mtxt="Dich - CRNP")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dyy,mtxt="Cat - MLR")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dyyy,mtxt="Cat - CRNP")
dev.off()

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

flexrecalk3(filenm = "MLR SCENARIO 4 ALTERNATIVE FLEXIBLE CALIBRATION MODELS.RData")

load("MLR SCENARIO 4 ALTERNATIVE FLEXIBLE CALIBRATION MODELS.RData")

round(cbind(ecifrm[1:4],ecifrm[13:16],ecifrm[5:8],ecifrm[17:20],ecifrm[9:12],ecifrm[21:24]), digits = 4)

par(mfrow=c(3,2))
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4azz,mtxt="Ref - MLR")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4azzz,mtxt="Ref - CRNP")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4axx,mtxt="Dich - MLR")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4axxx,mtxt="Dich - CRNP")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4ayy,mtxt="Cat - MLR")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4ayyy,mtxt="Cat - CRNP")
dev.off()
par(mfrow=c(3,2))
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4bzz,mtxt="Ref - MLR")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4bzzz,mtxt="Ref - CRNP")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4bxx,mtxt="Dich - MLR")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4bxxx,mtxt="Dich - CRNP")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4byy,mtxt="Cat - MLR")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4byyy,mtxt="Cat - CRNP")
dev.off()
par(mfrow=c(3,2))
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4czz,mtxt="Ref - MLR")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4czzz,mtxt="Ref - CRNP")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4cxx,mtxt="Dich - MLR")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4cxxx,mtxt="Dich - CRNP")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4cyy,mtxt="Cat - MLR")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4cyyy,mtxt="Cat - CRNP")
dev.off()
par(mfrow=c(3,2))
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dzz,mtxt="Ref - MLR")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dzzz,mtxt="Ref - CRNP")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dxx,mtxt="Dich - MLR")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dxxx,mtxt="Dich - CRNP")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dyy,mtxt="Cat - MLR")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dyyy,mtxt="Cat - CRNP")
dev.off()

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

flexrecalk3(filenm = "MLR SCENARIO 5 ALTERNATIVE FLEXIBLE CALIBRATION MODELS.RData")

load("MLR SCENARIO 5 ALTERNATIVE FLEXIBLE CALIBRATION MODELS.RData")

round(cbind(ecifrm[1:4],ecifrm[13:16],ecifrm[5:8],ecifrm[17:20],ecifrm[9:12],ecifrm[21:24]), digits = 4)

par(mfrow=c(3,2))
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4azz,mtxt="Ref - MLR")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4azzz,mtxt="Ref - CRNP")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4axx,mtxt="Dich - MLR")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4axxx,mtxt="Dich - CRNP")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4ayy,mtxt="Cat - MLR")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4ayyy,mtxt="Cat - CRNP")
dev.off()
par(mfrow=c(3,2))
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4bzz,mtxt="Ref - MLR")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4bzzz,mtxt="Ref - CRNP")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4bxx,mtxt="Dich - MLR")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4bxxx,mtxt="Dich - CRNP")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4byy,mtxt="Cat - MLR")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4byyy,mtxt="Cat - CRNP")
dev.off()
par(mfrow=c(3,2))
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4czz,mtxt="Ref - MLR")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4czzz,mtxt="Ref - CRNP")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4cxx,mtxt="Dich - MLR")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4cxxx,mtxt="Dich - CRNP")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4cyy,mtxt="Cat - MLR")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4cyyy,mtxt="Cat - CRNP")
dev.off()
par(mfrow=c(3,2))
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dzz,mtxt="Ref - MLR")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dzzz,mtxt="Ref - CRNP")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dxx,mtxt="Dich - MLR")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dxxx,mtxt="Dich - CRNP")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dyy,mtxt="Cat - MLR")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dyyy,mtxt="Cat - CRNP")
dev.off()

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

flexrecalk4(filenm = "MLR SCENARIO 6 ALTERNATIVE FLEXIBLE CALIBRATION MODELS.RData")

load("MLR SCENARIO 6 ALTERNATIVE FLEXIBLE CALIBRATION MODELS.RData")

round(cbind(ecifrm[1:4],ecifrm[13:16],ecifrm[5:8],ecifrm[17:20],ecifrm[9:12],ecifrm[21:24]), digits = 4)

par(mfrow=c(3,2))
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4azz,mtxt="Ref - MLR")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4azzz,mtxt="Ref - CRNP")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4axx,mtxt="Dich - MLR")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4axxx,mtxt="Dich - CRNP")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4ayy,mtxt="Cat - MLR")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4ayyy,mtxt="Cat - CRNP")
dev.off()
par(mfrow=c(3,2))
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4bzz,mtxt="Ref - MLR")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4bzzz,mtxt="Ref - CRNP")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4bxx,mtxt="Dich - MLR")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4bxxx,mtxt="Dich - CRNP")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4byy,mtxt="Cat - MLR")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4byyy,mtxt="Cat - CRNP")
dev.off()
par(mfrow=c(3,2))
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4czz,mtxt="Ref - MLR")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4czzz,mtxt="Ref - CRNP")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4cxx,mtxt="Dich - MLR")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4cxxx,mtxt="Dich - CRNP")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4cyy,mtxt="Cat - MLR")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4cyyy,mtxt="Cat - CRNP")
dev.off()
par(mfrow=c(3,2))
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dzz,mtxt="Ref - MLR")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dzzz,mtxt="Ref - CRNP")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dxx,mtxt="Dich - MLR")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dxxx,mtxt="Dich - CRNP")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dyy,mtxt="Cat - MLR")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dyyy,mtxt="Cat - CRNP")
dev.off()

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

flexrecalk4(filenm = "MLR SCENARIO 7 ALTERNATIVE FLEXIBLE CALIBRATION MODELS.RData")

load("MLR SCENARIO 7 ALTERNATIVE FLEXIBLE CALIBRATION MODELS.RData")

round(cbind(ecifrm[1:4],ecifrm[13:16],ecifrm[5:8],ecifrm[17:20],ecifrm[9:12],ecifrm[21:24]), digits = 4)

par(mfrow=c(3,2))
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4azz,mtxt="Ref - MLR")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4azzz,mtxt="Ref - CRNP")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4axx,mtxt="Dich - MLR")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4axxx,mtxt="Dich - CRNP")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4ayy,mtxt="Cat - MLR")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4ayyy,mtxt="Cat - CRNP")
dev.off()
par(mfrow=c(3,2))
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4bzz,mtxt="Ref - MLR")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4bzzz,mtxt="Ref - CRNP")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4bxx,mtxt="Dich - MLR")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4bxxx,mtxt="Dich - CRNP")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4byy,mtxt="Cat - MLR")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4byyy,mtxt="Cat - CRNP")
dev.off()
par(mfrow=c(3,2))
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4czz,mtxt="Ref - MLR")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4czzz,mtxt="Ref - CRNP")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4cxx,mtxt="Dich - MLR")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4cxxx,mtxt="Dich - CRNP")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4cyy,mtxt="Cat - MLR")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4cyyy,mtxt="Cat - CRNP")
dev.off()
par(mfrow=c(3,2))
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dzz,mtxt="Ref - MLR")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dzzz,mtxt="Ref - CRNP")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dxx,mtxt="Dich - MLR")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dxxx,mtxt="Dich - CRNP")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dyy,mtxt="Cat - MLR")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dyyy,mtxt="Cat - CRNP")
dev.off()

}


# SCENARIO 8

{
# Get large dataset
set.seed(34634) # the seed for the large datasets
datapo=getdata_k4p3(mu_x1=c(0,0,0.6,0.6),
                    mu_x2=c(0,0.4,0.4,0.5),
                    mu_x3=c(0.1,0,0.6,0.7),
                    yprev=c(45,30,20,5),#yprev=c(55,30,15),
                    mult=200)

# Run models for large dataset
mlrfit=vglm(y ~ x1 + x2 + x3, family=multinomial(refLevel = "1"), data=datapo[[1]])
acnpres=mlrfunc(datafit=datapo[[1]],dataval=datapo[[1]],ECI="nonpar",levels=4,modelformula=y ~ x1 + x2 + x3)
olrres=olrfunc(datafit=datapo[[1]],dataval=datapo[[1]],ECI="nonpar",levels=4,modelformula=y ~ x1 + x2 + x3)
acpres=acpfunc(datafit=datapo[[1]],dataval=datapo[[1]],ECI="nonpar",levels=4,modelformula=y ~ x1 + x2 + x3)
smres=smfunc(datafit=datapo[[1]],dataval=datapo[[1]],ECI="nonpar",levels=4,modelformula=y ~ x1 + x2 + x3)
#comsave(filep="Ordinal prediction model simulations\\ordsimk4p3_xnonequi_xunequal_loworc.RData")

flexrecalk4(filenm = "MLR SCENARIO 8 ALTERNATIVE FLEXIBLE CALIBRATION MODELS.RData")

load("MLR SCENARIO 8 ALTERNATIVE FLEXIBLE CALIBRATION MODELS.RData")

round(cbind(ecifrm[1:4],ecifrm[13:16],ecifrm[5:8],ecifrm[17:20],ecifrm[9:12],ecifrm[21:24]), digits = 4)

par(mfrow=c(3,2))
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4azz,mtxt="Ref - MLR")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4azzz,mtxt="Ref - CRNP")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4axx,mtxt="Dich - MLR")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4axxx,mtxt="Dich - CRNP")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4ayy,mtxt="Cat - MLR")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4ayyy,mtxt="Cat - CRNP")
dev.off()
par(mfrow=c(3,2))
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4bzz,mtxt="Ref - MLR")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4bzzz,mtxt="Ref - CRNP")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4bxx,mtxt="Dich - MLR")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4bxxx,mtxt="Dich - CRNP")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4byy,mtxt="Cat - MLR")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4byyy,mtxt="Cat - CRNP")
dev.off()
par(mfrow=c(3,2))
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4czz,mtxt="Ref - MLR")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4czzz,mtxt="Ref - CRNP")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4cxx,mtxt="Dich - MLR")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4cxxx,mtxt="Dich - CRNP")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4cyy,mtxt="Cat - MLR")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4cyyy,mtxt="Cat - CRNP")
dev.off()
par(mfrow=c(3,2))
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dzz,mtxt="Ref - MLR")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dzzz,mtxt="Ref - CRNP")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dxx,mtxt="Dich - MLR")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dxxx,mtxt="Dich - CRNP")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dyy,mtxt="Cat - MLR")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dyyy,mtxt="Cat - CRNP")
dev.off()

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

flexrecalk3(filenm = "MLR SCENARIO 9 ALTERNATIVE FLEXIBLE CALIBRATION MODELS.RData")

load("MLR SCENARIO 9 ALTERNATIVE FLEXIBLE CALIBRATION MODELS.RData")

round(cbind(ecifrm[1:4],ecifrm[13:16],ecifrm[5:8],ecifrm[17:20],ecifrm[9:12],ecifrm[21:24]), digits = 4)

par(mfrow=c(3,2))
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4azz,mtxt="Ref - MLR")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4azzz,mtxt="Ref - CRNP")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4axx,mtxt="Dich - MLR")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4axxx,mtxt="Dich - CRNP")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4ayy,mtxt="Cat - MLR")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4ayyy,mtxt="Cat - CRNP")
dev.off()
par(mfrow=c(3,2))
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4bzz,mtxt="Ref - MLR")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4bzzz,mtxt="Ref - CRNP")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4bxx,mtxt="Dich - MLR")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4bxxx,mtxt="Dich - CRNP")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4byy,mtxt="Cat - MLR")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4byyy,mtxt="Cat - CRNP")
dev.off()
par(mfrow=c(3,2))
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4czz,mtxt="Ref - MLR")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4czzz,mtxt="Ref - CRNP")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4cxx,mtxt="Dich - MLR")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4cxxx,mtxt="Dich - CRNP")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4cyy,mtxt="Cat - MLR")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4cyyy,mtxt="Cat - CRNP")
dev.off()
par(mfrow=c(3,2))
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dzz,mtxt="Ref - MLR")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dzzz,mtxt="Ref - CRNP")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dxx,mtxt="Dich - MLR")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dxxx,mtxt="Dich - CRNP")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dyy,mtxt="Cat - MLR")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dyyy,mtxt="Cat - CRNP")
dev.off()

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

flexrecalk4(filenm = "MLR SCENARIO 10 ALTERNATIVE FLEXIBLE CALIBRATION MODELS.RData")

load("MLR SCENARIO 10 ALTERNATIVE FLEXIBLE CALIBRATION MODELS.RData")

round(cbind(ecifrm[1:4],ecifrm[13:16],ecifrm[5:8],ecifrm[17:20],ecifrm[9:12],ecifrm[21:24]), digits = 4)

par(mfrow=c(3,2))
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4azz,mtxt="Ref - MLR")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4azzz,mtxt="Ref - CRNP")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4axx,mtxt="Dich - MLR")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4axxx,mtxt="Dich - CRNP")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4ayy,mtxt="Cat - MLR")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4ayyy,mtxt="Cat - CRNP")
dev.off()
par(mfrow=c(3,2))
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4bzz,mtxt="Ref - MLR")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4bzzz,mtxt="Ref - CRNP")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4bxx,mtxt="Dich - MLR")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4bxxx,mtxt="Dich - CRNP")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4byy,mtxt="Cat - MLR")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4byyy,mtxt="Cat - CRNP")
dev.off()
par(mfrow=c(3,2))
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4czz,mtxt="Ref - MLR")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4czzz,mtxt="Ref - CRNP")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4cxx,mtxt="Dich - MLR")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4cxxx,mtxt="Dich - CRNP")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4cyy,mtxt="Cat - MLR")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4cyyy,mtxt="Cat - CRNP")
dev.off()
par(mfrow=c(3,2))
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dzz,mtxt="Ref - MLR")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dzzz,mtxt="Ref - CRNP")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dxx,mtxt="Dich - MLR")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dxxx,mtxt="Dich - CRNP")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dyy,mtxt="Cat - MLR")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dyyy,mtxt="Cat - CRNP")
dev.off()

}


# SCENARIO 11 (i.e. scen 4 with 4 noise predictors added)

{
# Get large dataset
set.seed(34634) # the seed for the large datasets
datapo=getdata_k3p4(mu_x1=c(0,0.7,0.8),
                    mu_x2=c(0,0.6,0.6),
                    mu_x3=c(0,0.5,0.8),
                    mu_x4=c(0,0.1,0.6),
                    yprev=c(55,30,15),#yprev=c(55,30,15),
                    mult=2000)
  
noise = as.data.frame(matrix(rnorm(dim(datapo[[1]])[1] * 4, mean = 0, sd = 1), ncol = 4))
colnames(noise) = c("x5", "x6", "x7", "x8")
datanoise = cbind(datapo[[1]][1:5], noise, datapo[[1]][6:8])
  
# Run models for large dataset
mlrfit=vglm(y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8, family=multinomial(refLevel = "1"), data=datanoise)
acnpres=mlrfunc(datafit=datanoise,dataval=datanoise,ECI="nonpar",levels=3,modelformula=y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8)
olrres=olrfunc(datafit=datanoise,dataval=datanoise,ECI="nonpar",levels=3,modelformula=y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8)
acpres=acpfunc(datafit=datanoise,dataval=datanoise,ECI="nonpar",levels=3,modelformula=y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8)
smres=smfunc(datafit=datanoise,dataval=datanoise,ECI="nonpar",levels=3,modelformula=y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8)

flexrecalk3(filenm = "MLR SCENARIO 11 ALTERNATIVE FLEXIBLE CALIBRATION MODELS.RData")

load("MLR SCENARIO 11 ALTERNATIVE FLEXIBLE CALIBRATION MODELS.RData")

round(cbind(ecifrm[1:4],ecifrm[13:16],ecifrm[5:8],ecifrm[17:20],ecifrm[9:12],ecifrm[21:24]), digits = 4)

par(mfrow=c(3,2))
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4azz,mtxt="Ref - MLR")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4azzz,mtxt="Ref - CRNP")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4axx,mtxt="Dich - MLR")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4axxx,mtxt="Dich - CRNP")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4ayy,mtxt="Cat - MLR")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4ayyy,mtxt="Cat - CRNP")
dev.off()
par(mfrow=c(3,2))
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4bzz,mtxt="Ref - MLR")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4bzzz,mtxt="Ref - CRNP")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4bxx,mtxt="Dich - MLR")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4bxxx,mtxt="Dich - CRNP")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4byy,mtxt="Cat - MLR")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4byyy,mtxt="Cat - CRNP")
dev.off()
par(mfrow=c(3,2))
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4czz,mtxt="Ref - MLR")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4czzz,mtxt="Ref - CRNP")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4cxx,mtxt="Dich - MLR")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4cxxx,mtxt="Dich - CRNP")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4cyy,mtxt="Cat - MLR")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4cyyy,mtxt="Cat - CRNP")
dev.off()
par(mfrow=c(3,2))
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dzz,mtxt="Ref - MLR")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dzzz,mtxt="Ref - CRNP")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dxx,mtxt="Dich - MLR")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dxxx,mtxt="Dich - CRNP")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dyy,mtxt="Cat - MLR")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dyyy,mtxt="Cat - CRNP")
dev.off()

}




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

flexrecalk3(filenm = "CLPO SCENARIO 1 ALTERNATIVE FLEXIBLE CALIBRATION MODELS.RData")

load("CLPO SCENARIO 1 ALTERNATIVE FLEXIBLE CALIBRATION MODELS.RData")

round(cbind(ecifrm[1:4],ecifrm[13:16],ecifrm[5:8],ecifrm[17:20],ecifrm[9:12],ecifrm[21:24]), digits = 4)

par(mfrow=c(3,2))
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4azz,mtxt="Ref - MLR")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4azzz,mtxt="Ref - CRNP")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4axx,mtxt="Dich - MLR")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4axxx,mtxt="Dich - CRNP")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4ayy,mtxt="Cat - MLR")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4ayyy,mtxt="Cat - CRNP")
dev.off()
par(mfrow=c(3,2))
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4bzz,mtxt="Ref - MLR")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4bzzz,mtxt="Ref - CRNP")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4bxx,mtxt="Dich - MLR")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4bxxx,mtxt="Dich - CRNP")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4byy,mtxt="Cat - MLR")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4byyy,mtxt="Cat - CRNP")
dev.off()
par(mfrow=c(3,2))
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4czz,mtxt="Ref - MLR")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4czzz,mtxt="Ref - CRNP")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4cxx,mtxt="Dich - MLR")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4cxxx,mtxt="Dich - CRNP")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4cyy,mtxt="Cat - MLR")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4cyyy,mtxt="Cat - CRNP")
dev.off()
par(mfrow=c(3,2))
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dzz,mtxt="Ref - MLR")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dzzz,mtxt="Ref - CRNP")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dxx,mtxt="Dich - MLR")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dxxx,mtxt="Dich - CRNP")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dyy,mtxt="Cat - MLR")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dyyy,mtxt="Cat - CRNP")
dev.off()

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

flexrecalk3(filenm = "CLPO SCENARIO 2 ALTERNATIVE FLEXIBLE CALIBRATION MODELS.RData")

load("CLPO SCENARIO 2 ALTERNATIVE FLEXIBLE CALIBRATION MODELS.RData")

round(cbind(ecifrm[1:4],ecifrm[13:16],ecifrm[5:8],ecifrm[17:20],ecifrm[9:12],ecifrm[21:24]), digits = 4)

par(mfrow=c(3,2))
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4azz,mtxt="Ref - MLR")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4azzz,mtxt="Ref - CRNP")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4axx,mtxt="Dich - MLR")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4axxx,mtxt="Dich - CRNP")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4ayy,mtxt="Cat - MLR")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4ayyy,mtxt="Cat - CRNP")
dev.off()
par(mfrow=c(3,2))
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4bzz,mtxt="Ref - MLR")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4bzzz,mtxt="Ref - CRNP")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4bxx,mtxt="Dich - MLR")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4bxxx,mtxt="Dich - CRNP")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4byy,mtxt="Cat - MLR")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4byyy,mtxt="Cat - CRNP")
dev.off()
par(mfrow=c(3,2))
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4czz,mtxt="Ref - MLR")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4czzz,mtxt="Ref - CRNP")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4cxx,mtxt="Dich - MLR")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4cxxx,mtxt="Dich - CRNP")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4cyy,mtxt="Cat - MLR")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4cyyy,mtxt="Cat - CRNP")
dev.off()
par(mfrow=c(3,2))
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dzz,mtxt="Ref - MLR")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dzzz,mtxt="Ref - CRNP")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dxx,mtxt="Dich - MLR")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dxxx,mtxt="Dich - CRNP")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dyy,mtxt="Cat - MLR")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dyyy,mtxt="Cat - CRNP")
dev.off()

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

flexrecalk3(filenm = "CLPO SCENARIO 3 ALTERNATIVE FLEXIBLE CALIBRATION MODELS.RData")

load("CLPO SCENARIO 3 ALTERNATIVE FLEXIBLE CALIBRATION MODELS.RData")

round(cbind(ecifrm[1:4],ecifrm[13:16],ecifrm[5:8],ecifrm[17:20],ecifrm[9:12],ecifrm[21:24]), digits = 4)

par(mfrow=c(3,2))
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4azz,mtxt="Ref - MLR")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4azzz,mtxt="Ref - CRNP")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4axx,mtxt="Dich - MLR")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4axxx,mtxt="Dich - CRNP")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4ayy,mtxt="Cat - MLR")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4ayyy,mtxt="Cat - CRNP")
dev.off()
par(mfrow=c(3,2))
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4bzz,mtxt="Ref - MLR")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4bzzz,mtxt="Ref - CRNP")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4bxx,mtxt="Dich - MLR")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4bxxx,mtxt="Dich - CRNP")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4byy,mtxt="Cat - MLR")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4byyy,mtxt="Cat - CRNP")
dev.off()
par(mfrow=c(3,2))
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4czz,mtxt="Ref - MLR")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4czzz,mtxt="Ref - CRNP")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4cxx,mtxt="Dich - MLR")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4cxxx,mtxt="Dich - CRNP")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4cyy,mtxt="Cat - MLR")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4cyyy,mtxt="Cat - CRNP")
dev.off()
par(mfrow=c(3,2))
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dzz,mtxt="Ref - MLR")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dzzz,mtxt="Ref - CRNP")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dxx,mtxt="Dich - MLR")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dxxx,mtxt="Dich - CRNP")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dyy,mtxt="Cat - MLR")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dyyy,mtxt="Cat - CRNP")
dev.off()

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

flexrecalk4(filenm = "CLPO SCENARIO 4 ALTERNATIVE FLEXIBLE CALIBRATION MODELS.RData")

load("CLPO SCENARIO 4 ALTERNATIVE FLEXIBLE CALIBRATION MODELS.RData")

round(cbind(ecifrm[1:4],ecifrm[13:16],ecifrm[5:8],ecifrm[17:20],ecifrm[9:12],ecifrm[21:24]), digits = 4)

par(mfrow=c(3,2))
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4azz,mtxt="Ref - MLR")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4azzz,mtxt="Ref - CRNP")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4axx,mtxt="Dich - MLR")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4axxx,mtxt="Dich - CRNP")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4ayy,mtxt="Cat - MLR")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4ayyy,mtxt="Cat - CRNP")
dev.off()
par(mfrow=c(3,2))
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4bzz,mtxt="Ref - MLR")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4bzzz,mtxt="Ref - CRNP")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4bxx,mtxt="Dich - MLR")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4bxxx,mtxt="Dich - CRNP")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4byy,mtxt="Cat - MLR")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4byyy,mtxt="Cat - CRNP")
dev.off()
par(mfrow=c(3,2))
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4czz,mtxt="Ref - MLR")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4czzz,mtxt="Ref - CRNP")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4cxx,mtxt="Dich - MLR")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4cxxx,mtxt="Dich - CRNP")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4cyy,mtxt="Cat - MLR")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4cyyy,mtxt="Cat - CRNP")
dev.off()
par(mfrow=c(3,2))
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dzz,mtxt="Ref - MLR")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dzzz,mtxt="Ref - CRNP")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dxx,mtxt="Dich - MLR")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dxxx,mtxt="Dich - CRNP")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dyy,mtxt="Cat - MLR")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dyyy,mtxt="Cat - CRNP")
dev.off()

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

flexrecalk4(filenm = "CLPO SCENARIO 5 ALTERNATIVE FLEXIBLE CALIBRATION MODELS.RData")

load("CLPO SCENARIO 5 ALTERNATIVE FLEXIBLE CALIBRATION MODELS.RData")

round(cbind(ecifrm[1:4],ecifrm[13:16],ecifrm[5:8],ecifrm[17:20],ecifrm[9:12],ecifrm[21:24]), digits = 4)

par(mfrow=c(3,2))
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4azz,mtxt="Ref - MLR")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4azzz,mtxt="Ref - CRNP")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4axx,mtxt="Dich - MLR")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4axxx,mtxt="Dich - CRNP")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4ayy,mtxt="Cat - MLR")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4ayyy,mtxt="Cat - CRNP")
dev.off()
par(mfrow=c(3,2))
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4bzz,mtxt="Ref - MLR")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4bzzz,mtxt="Ref - CRNP")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4bxx,mtxt="Dich - MLR")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4bxxx,mtxt="Dich - CRNP")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4byy,mtxt="Cat - MLR")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4byyy,mtxt="Cat - CRNP")
dev.off()
par(mfrow=c(3,2))
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4czz,mtxt="Ref - MLR")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4czzz,mtxt="Ref - CRNP")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4cxx,mtxt="Dich - MLR")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4cxxx,mtxt="Dich - CRNP")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4cyy,mtxt="Cat - MLR")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4cyyy,mtxt="Cat - CRNP")
dev.off()
par(mfrow=c(3,2))
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dzz,mtxt="Ref - MLR")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dzzz,mtxt="Ref - CRNP")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dxx,mtxt="Dich - MLR")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dxxx,mtxt="Dich - CRNP")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dyy,mtxt="Cat - MLR")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dyyy,mtxt="Cat - CRNP")
dev.off()

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

flexrecalk4(filenm = "CLPO SCENARIO 6 ALTERNATIVE FLEXIBLE CALIBRATION MODELS.RData")

load("CLPO SCENARIO 6 ALTERNATIVE FLEXIBLE CALIBRATION MODELS.RData")

round(cbind(ecifrm[1:4],ecifrm[13:16],ecifrm[5:8],ecifrm[17:20],ecifrm[9:12],ecifrm[21:24]), digits = 4)

par(mfrow=c(3,2))
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4azz,mtxt="Ref - MLR")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4azzz,mtxt="Ref - CRNP")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4axx,mtxt="Dich - MLR")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4axxx,mtxt="Dich - CRNP")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4ayy,mtxt="Cat - MLR")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4ayyy,mtxt="Cat - CRNP")
dev.off()
par(mfrow=c(3,2))
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4bzz,mtxt="Ref - MLR")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4bzzz,mtxt="Ref - CRNP")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4bxx,mtxt="Dich - MLR")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4bxxx,mtxt="Dich - CRNP")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4byy,mtxt="Cat - MLR")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4byyy,mtxt="Cat - CRNP")
dev.off()
par(mfrow=c(3,2))
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4czz,mtxt="Ref - MLR")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4czzz,mtxt="Ref - CRNP")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4cxx,mtxt="Dich - MLR")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4cxxx,mtxt="Dich - CRNP")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4cyy,mtxt="Cat - MLR")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4cyyy,mtxt="Cat - CRNP")
dev.off()
par(mfrow=c(3,2))
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dzz,mtxt="Ref - MLR")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dzzz,mtxt="Ref - CRNP")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dxx,mtxt="Dich - MLR")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dxxx,mtxt="Dich - CRNP")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dyy,mtxt="Cat - MLR")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dyyy,mtxt="Cat - CRNP")
dev.off()

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

flexrecalk3(filenm = "CLPO SCENARIO 7 ALTERNATIVE FLEXIBLE CALIBRATION MODELS.RData")

load("CLPO SCENARIO 7 ALTERNATIVE FLEXIBLE CALIBRATION MODELS.RData")

round(cbind(ecifrm[1:4],ecifrm[13:16],ecifrm[5:8],ecifrm[17:20],ecifrm[9:12],ecifrm[21:24]), digits = 4)

par(mfrow=c(3,2))
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4azz,mtxt="Ref - MLR")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4azzz,mtxt="Ref - CRNP")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4axx,mtxt="Dich - MLR")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4axxx,mtxt="Dich - CRNP")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4ayy,mtxt="Cat - MLR")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4ayyy,mtxt="Cat - CRNP")
dev.off()
par(mfrow=c(3,2))
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4bzz,mtxt="Ref - MLR")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4bzzz,mtxt="Ref - CRNP")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4bxx,mtxt="Dich - MLR")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4bxxx,mtxt="Dich - CRNP")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4byy,mtxt="Cat - MLR")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4byyy,mtxt="Cat - CRNP")
dev.off()
par(mfrow=c(3,2))
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4czz,mtxt="Ref - MLR")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4czzz,mtxt="Ref - CRNP")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4cxx,mtxt="Dich - MLR")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4cxxx,mtxt="Dich - CRNP")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4cyy,mtxt="Cat - MLR")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4cyyy,mtxt="Cat - CRNP")
dev.off()
par(mfrow=c(3,2))
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dzz,mtxt="Ref - MLR")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dzzz,mtxt="Ref - CRNP")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dxx,mtxt="Dich - MLR")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dxxx,mtxt="Dich - CRNP")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dyy,mtxt="Cat - MLR")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dyyy,mtxt="Cat - CRNP")
dev.off()

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

flexrecalk4(filenm = "CLPO SCENARIO 8 ALTERNATIVE FLEXIBLE CALIBRATION MODELS.RData")

load("CLPO SCENARIO 8 ALTERNATIVE FLEXIBLE CALIBRATION MODELS.RData")

round(cbind(ecifrm[1:4],ecifrm[13:16],ecifrm[5:8],ecifrm[17:20],ecifrm[9:12],ecifrm[21:24]), digits = 4)

par(mfrow=c(3,2))
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4azz,mtxt="Ref - MLR")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4azzz,mtxt="Ref - CRNP")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4axx,mtxt="Dich - MLR")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4axxx,mtxt="Dich - CRNP")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4ayy,mtxt="Cat - MLR")
plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4ayyy,mtxt="Cat - CRNP")
dev.off()
par(mfrow=c(3,2))
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4bzz,mtxt="Ref - MLR")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4bzzz,mtxt="Ref - CRNP")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4bxx,mtxt="Dich - MLR")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4bxxx,mtxt="Dich - CRNP")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4byy,mtxt="Cat - MLR")
plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4byyy,mtxt="Cat - CRNP")
dev.off()
par(mfrow=c(3,2))
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4czz,mtxt="Ref - MLR")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4czzz,mtxt="Ref - CRNP")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4cxx,mtxt="Dich - MLR")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4cxxx,mtxt="Dich - CRNP")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4cyy,mtxt="Cat - MLR")
plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4cyyy,mtxt="Cat - CRNP")
dev.off()
par(mfrow=c(3,2))
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dzz,mtxt="Ref - MLR")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dzzz,mtxt="Ref - CRNP")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dxx,mtxt="Dich - MLR")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dxxx,mtxt="Dich - CRNP")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dyy,mtxt="Cat - MLR")
plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dyyy,mtxt="Cat - CRNP")
dev.off()

}


# SCENARIO 9 (scenario 1 with 4 noise predictors)

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
  
  noise = as.data.frame(matrix(rnorm(dim(datapo[[1]])[1] * 4, mean = 0, sd = 1), ncol = 4))
  colnames(noise) = c("x5", "x6", "x7", "x8")
  datanoise = cbind(datapo[[1]][1:5], noise, datapo[[1]][6:8])
  
  # Run models for large dataset
  mlrfit=vglm(y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8, family=multinomial(refLevel = "1"), data=datanoise)
  acnpres=mlrfunc(datafit=datanoise,dataval=datanoise,ECI="nonpar",levels=3,modelformula=y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8)
  olrres=olrfunc(datafit=datanoise,dataval=datanoise,ECI="nonpar",levels=3,modelformula=y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8)
  acpres=acpfunc(datafit=datanoise,dataval=datanoise,ECI="nonpar",levels=3,modelformula=y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8)
  smres=smfunc(datafit=datanoise,dataval=datanoise,ECI="nonpar",levels=3,modelformula=y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8)
  
  flexrecalk3(filenm = "CLPO SCENARIO 9 ALTERNATIVE FLEXIBLE CALIBRATION MODELS.RData")
  
  load("CLPO SCENARIO 9 ALTERNATIVE FLEXIBLE CALIBRATION MODELS.RData")
  
  round(cbind(ecifrm[1:4],ecifrm[13:16],ecifrm[5:8],ecifrm[17:20],ecifrm[9:12],ecifrm[21:24]), digits = 4)
  
  par(mfrow=c(3,2))
  plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4azz,mtxt="Ref - MLR")
  plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4azzz,mtxt="Ref - CRNP")
  plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4axx,mtxt="Dich - MLR")
  plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4axxx,mtxt="Dich - CRNP")
  plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4ayy,mtxt="Cat - MLR")
  plotscatter4sub(preds=acnpres[[3]],obs=vgamsmps4ayyy,mtxt="Cat - CRNP")
  dev.off()
  par(mfrow=c(3,2))
  plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4bzz,mtxt="Ref - MLR")
  plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4bzzz,mtxt="Ref - CRNP")
  plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4bxx,mtxt="Dich - MLR")
  plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4bxxx,mtxt="Dich - CRNP")
  plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4byy,mtxt="Cat - MLR")
  plotscatter4sub(preds=olrres[[3]],obs=vgamsmps4byyy,mtxt="Cat - CRNP")
  dev.off()
  par(mfrow=c(3,2))
  plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4czz,mtxt="Ref - MLR")
  plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4czzz,mtxt="Ref - CRNP")
  plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4cxx,mtxt="Dich - MLR")
  plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4cxxx,mtxt="Dich - CRNP")
  plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4cyy,mtxt="Cat - MLR")
  plotscatter4sub(preds=acpres[[3]],obs=vgamsmps4cyyy,mtxt="Cat - CRNP")
  dev.off()
  par(mfrow=c(3,2))
  plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dzz,mtxt="Ref - MLR")
  plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dzzz,mtxt="Ref - CRNP")
  plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dxx,mtxt="Dich - MLR")
  plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dxxx,mtxt="Dich - CRNP")
  plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dyy,mtxt="Cat - MLR")
  plotscatter4sub(preds=smres[[3]],obs=vgamsmps4dyyy,mtxt="Cat - CRNP")
  dev.off()
  
}  