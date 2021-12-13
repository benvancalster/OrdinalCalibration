library(nnet) # multinom
library(gmodels) # CrossTable
library(psych) # describeBy
library(MASS) # polr
library(VGAM) # ordinal code
library(Hmisc)
library(pROC)
library(clusterPower) # expit
#install.packages("multiCalibration", repos = "http://repos.openanalytics.eu", type = "source")
library(multiCalibration)
library(ordinal)
library(DescTools) # Cstat
library(boot) # boot
library(foreign) # read.dta
library(rms) # plot.xmean.ordinaly
library(brant) # brant test for PO in CL-PO models fitted using polr

# Set working directory to save files
setwd("C:\\Ben\\Methodological research\\Ordinal prediction models")

# Load functions
source("R code\\ordcalfunctions.R")
source("R code\\ordcalfunctions_bootstrap.R")



# CARDIIGAN DATA

CAD5=read.dta("C:\\Ben\\Methodological research\\Ordinal prediction CAD\\CADimputed04.dta", convert.factors = FALSE)
cad = CAD5[CAD5$v_mi_m==1,]
cad$crp_2d=1*(cad$crp_2>=0.5)
cad$o5 = ordered(cad$CADextensive+1)


###################################
# MULTINOMIAL LOGISTIC REGRESSION #
###################################

# Fit the model
mlr <- vglm(o5 ~ agexact + male + ap_3 + diab + hypert_m + dyslip + smokstat_m + hdlchol + ldlchol + fibr_ln + crp_2d, family=multinomial(refLevel = "1"), data=cad)

# Get predicted probabilities, the linear predictors
mlrpred <- predictvglm(mlr,newdata=cad,type="response")
mlrlpred <- predictvglm(mlr,newdata=cad,type="link")

# Flexible recalibration models
flexrecalk5_cardiigan(filenm = "CARDIIGAN MLR ALTERNATIVE FLEXIBLE CALIBRATION MODELS.Rdata", res = mlrpred)

load("CARDIIGAN MLR ALTERNATIVE FLEXIBLE CALIBRATION MODELS.Rdata")

# Calibration scatter plots and curves per outcome category and per outcome dichotomy
#pdf("CARDIIGAN_MLR_Calibrationplots_comparison_scatter.pdf")
par(mfrow=c(3,2))
plotscatter5(preds=mlrpred,obs=vgamsmps4azz,cexaxisv=1,cexlabv=1.2,cexlegendv=0,maintxt="MLR - Ref LP")
plotscatter5(preds=mlrpred,obs=vgamsmps4azzz,cexaxisv=1,cexlabv=1.2,cexlegendv=0,maintxt="CR-NP - Ref LP")
plotscatter5(preds=mlrpred,obs=vgamsmps4axx,cexaxisv=1,cexlabv=1.2,cexlegendv=0,maintxt="MLR - Dich LP")
plotscatter5(preds=mlrpred,obs=vgamsmps4axxx,cexaxisv=1,cexlabv=1.2,cexlegendv=0,maintxt="CR-NP - Dich LP")
plotscatter5(preds=mlrpred,obs=vgamsmps4ayy,cexaxisv=1,cexlabv=1.2,cexlegendv=0,maintxt="MLR - Cat LP")
plotscatter5(preds=mlrpred,obs=vgamsmps4ayyy,cexaxisv=1,cexlabv=1.2,cexlegendv=0,maintxt="CR-NP - Cat LP")
dev.off()



###########################################
# CUMULATIVE LOGIT WITH PROPORTIONAL ODDS # 
###########################################

# Fit the model
olr <- vglm(o5 ~ agexact + male + ap_3 + diab + hypert_m + dyslip + smokstat_m + hdlchol + ldlchol + fibr_ln + crp_2d, family=cumulative(parallel=T, reverse=T), data=cad)

# Get predicted probabilities, the linear predictors
olrpred <- predictvglm(olr,newdata=cad,type="response")
olrlpred <- predictvglm(olr,newdata=cad,type="link")

# Flexible recalibration models
flexrecalk5_cardiigan(filenm = "CARDIIGAN CLPO ALTERNATIVE FLEXIBLE CALIBRATION MODELS.Rdata", res = olrpred)

load("CARDIIGAN CLPO ALTERNATIVE FLEXIBLE CALIBRATION MODELS.Rdata")

# Calibration scatter plots and curves per outcome category and per outcome dichotomy
#pdf("CARDIIGAN_CLPO_Calibrationplots_comparison_scatter.pdf")
par(mfrow=c(3,2))
plotscatter5(preds=olrpred,obs=vgamsmps4azz,cexaxisv=1,cexlabv=1.2,cexlegendv=0,maintxt="MLR - Ref LP")
plotscatter5(preds=olrpred,obs=vgamsmps4azzz,cexaxisv=1,cexlabv=1.2,cexlegendv=0,maintxt="CR-NP - Ref LP")
plotscatter5(preds=olrpred,obs=vgamsmps4axx,cexaxisv=1,cexlabv=1.2,cexlegendv=0,maintxt="MLR - Dich LP")
plotscatter5(preds=olrpred,obs=vgamsmps4axxx,cexaxisv=1,cexlabv=1.2,cexlegendv=0,maintxt="CR-NP - Dich LP")
plotscatter5(preds=olrpred,obs=vgamsmps4ayy,cexaxisv=1,cexlabv=1.2,cexlegendv=0,maintxt="MLR - Cat LP")
plotscatter5(preds=olrpred,obs=vgamsmps4ayyy,cexaxisv=1,cexlabv=1.2,cexlegendv=0,maintxt="CR-NP - Cat LP")
dev.off()


##################################################
# ADJACENT CATEGORY LOGIT WITH PROPORTIONAL ODDS #
##################################################

# Fit the model
acp <- vglm(o5 ~ agexact + male + ap_3 + diab + hypert_m + dyslip + smokstat_m + hdlchol + ldlchol + fibr_ln + crp_2d, family=acat(parallel=T), data=cad)

# Get predicted probabilities, the linear predictors
acppred <- predictvglm(acp,newdata=cad,type="response")
acplpred <- predictvglm(acp,newdata=cad,type="link")

# Flexible recalibration models
flexrecalk5_cardiigan(filenm = "CARDIIGAN ACPO ALTERNATIVE FLEXIBLE CALIBRATION MODELS.Rdata", res = acppred)

load("CARDIIGAN ACPO ALTERNATIVE FLEXIBLE CALIBRATION MODELS.Rdata")

# Calibration scatter plots and curves per outcome category and per outcome dichotomy
#pdf("CARDIIGAN_ACPO_Calibrationplots_comparison_scatter.pdf")
par(mfrow=c(3,2))
plotscatter5(preds=acppred,obs=vgamsmps4azz,cexaxisv=1,cexlabv=1.2,cexlegendv=0,maintxt="MLR - Ref LP")
plotscatter5(preds=acppred,obs=vgamsmps4azzz,cexaxisv=1,cexlabv=1.2,cexlegendv=0,maintxt="CR-NP - Ref LP")
plotscatter5(preds=acppred,obs=vgamsmps4axx,cexaxisv=1,cexlabv=1.2,cexlegendv=0,maintxt="MLR - Dich LP")
plotscatter5(preds=acppred,obs=vgamsmps4axxx,cexaxisv=1,cexlabv=1.2,cexlegendv=0,maintxt="CR-NP - Dich LP")
plotscatter5(preds=acppred,obs=vgamsmps4ayy,cexaxisv=1,cexlabv=1.2,cexlegendv=0,maintxt="MLR - Cat LP")
plotscatter5(preds=acppred,obs=vgamsmps4ayyy,cexaxisv=1,cexlabv=1.2,cexlegendv=0,maintxt="CR-NP - Cat LP")
dev.off()



###################################################
# CONTINUATION RATIO LOGIT WITH PROPORTIONAL ODDS #
###################################################

# Fit the model
crp <- vglm(o5 ~ agexact + male + ap_3 + diab + hypert_m + dyslip + smokstat_m + hdlchol + ldlchol + fibr_ln + crp_2d, family=cratio(parallel=T), data=cad)

# Get predicted probabilities, the linear predictors
crppred <- predictvglm(crp,newdata=cad,type="response")
crplpred <- predictvglm(crp,newdata=cad,type="link")

# Flexible recalibration models
flexrecalk5_cardiigan(filenm = "CARDIIGAN CRPO ALTERNATIVE FLEXIBLE CALIBRATION MODELS.Rdata", res = crppred)

load("CARDIIGAN CRPO ALTERNATIVE FLEXIBLE CALIBRATION MODELS.Rdata")

# Calibration scatter plots and curves per outcome category and per outcome dichotomy
#pdf("CARDIIGAN_CRPO_Calibrationplots_comparison_scatter.pdf")
par(mfrow=c(3,2))
plotscatter5(preds=crppred,obs=vgamsmps4azz,cexaxisv=1,cexlabv=1.2,cexlegendv=0,maintxt="MLR - Ref LP")
plotscatter5(preds=crppred,obs=vgamsmps4azzz,cexaxisv=1,cexlabv=1.2,cexlegendv=0,maintxt="CR-NP - Ref LP")
plotscatter5(preds=crppred,obs=vgamsmps4axx,cexaxisv=1,cexlabv=1.2,cexlegendv=0,maintxt="MLR - Dich LP")
plotscatter5(preds=crppred,obs=vgamsmps4axxx,cexaxisv=1,cexlabv=1.2,cexlegendv=0,maintxt="CR-NP - Dich LP")
plotscatter5(preds=crppred,obs=vgamsmps4ayy,cexaxisv=1,cexlabv=1.2,cexlegendv=0,maintxt="MLR - Cat LP")
plotscatter5(preds=crppred,obs=vgamsmps4ayyy,cexaxisv=1,cexlabv=1.2,cexlegendv=0,maintxt="CR-NP - Cat LP")
dev.off()



######################################################
# CONTINUATION RATIO LOGIT WITHOUT PROPORTIONAL ODDS #
######################################################

# Fit the model
crnp <- vglm(o5 ~ agexact + male + ap_3 + diab + hypert_m + dyslip + smokstat_m + hdlchol + ldlchol + fibr_ln + crp_2d, family=cratio(parallel=F), data=cad)

# Get predicted probabilities, the linear predictors
crnppred <- predictvglm(crnp,newdata=cad,type="response")
crnplpred <- predictvglm(crnp,newdata=cad,type="link")

# Flexible recalibration models
flexrecalk5_cardiigan(filenm = "CARDIIGAN CRNP ALTERNATIVE FLEXIBLE CALIBRATION MODELS.Rdata", res = crnppred)

load("CARDIIGAN CRNP ALTERNATIVE FLEXIBLE CALIBRATION MODELS.Rdata")

# Calibration scatter plots and curves per outcome category and per outcome dichotomy
#pdf("CARDIIGAN_CRNP_Calibrationplots_comparison_scatter.pdf")
par(mfrow=c(3,2))
plotscatter5(preds=crnppred,obs=vgamsmps4azz,cexaxisv=1,cexlabv=1.2,cexlegendv=0,maintxt="MLR - Ref LP")
plotscatter5(preds=crnppred,obs=vgamsmps4azzz,cexaxisv=1,cexlabv=1.2,cexlegendv=0,maintxt="CR-NP - Ref LP")
plotscatter5(preds=crnppred,obs=vgamsmps4axx,cexaxisv=1,cexlabv=1.2,cexlegendv=0,maintxt="MLR - Dich LP")
plotscatter5(preds=crnppred,obs=vgamsmps4axxx,cexaxisv=1,cexlabv=1.2,cexlegendv=0,maintxt="CR-NP - Dich LP")
plotscatter5(preds=crnppred,obs=vgamsmps4ayy,cexaxisv=1,cexlabv=1.2,cexlegendv=0,maintxt="MLR - Cat LP")
plotscatter5(preds=crnppred,obs=vgamsmps4ayyy,cexaxisv=1,cexlabv=1.2,cexlegendv=0,maintxt="CR-NP - Cat LP")
dev.off()



###############################################
# STEREOTYPE MULTINOMIAL LOGISTIC REGRESSSION #
###############################################

# Fit the model
sm=rrvglm(o5 ~ agexact + male + ap_3 + diab + hypert_m + dyslip + smokstat_m + hdlchol + ldlchol + fibr_ln + crp_2d, multinomial(refLevel = "1"), data = cad)

# Get predicted probabilities, the linear predictors, and the expected outcome E
smpred <- predictvglm(sm,newdata=cad,type="response")
smlpred <- predictvglm(sm,newdata=cad,type="link")

# Flexible recalibration models
flexrecalk5_cardiigan(filenm = "CARDIIGAN SLM ALTERNATIVE FLEXIBLE CALIBRATION MODELS.Rdata", res = smpred)

load("CARDIIGAN SLM ALTERNATIVE FLEXIBLE CALIBRATION MODELS.Rdata")

# Calibration scatter plots and curves per outcome category and per outcome dichotomy
#pdf("CARDIIGAN_SLM_Calibrationplots_comparison_scatter.pdf")
par(mfrow=c(3,2))
plotscatter5(preds=smpred,obs=vgamsmps4azz,cexaxisv=1,cexlabv=1.2,cexlegendv=0,maintxt="MLR - Ref LP")
plotscatter5(preds=smpred,obs=vgamsmps4azzz,cexaxisv=1,cexlabv=1.2,cexlegendv=0,maintxt="CR-NP - Ref LP")
plotscatter5(preds=smpred,obs=vgamsmps4axx,cexaxisv=1,cexlabv=1.2,cexlegendv=0,maintxt="MLR - Dich LP")
plotscatter5(preds=smpred,obs=vgamsmps4axxx,cexaxisv=1,cexlabv=1.2,cexlegendv=0,maintxt="CR-NP - Dich LP")
plotscatter5(preds=smpred,obs=vgamsmps4ayy,cexaxisv=1,cexlabv=1.2,cexlegendv=0,maintxt="MLR - Cat LP")
plotscatter5(preds=smpred,obs=vgamsmps4ayyy,cexaxisv=1,cexlabv=1.2,cexlegendv=0,maintxt="CR-NP - Cat LP")
dev.off()
