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
source("ordcalfunctions.R")
source("ordcalfunctions_bootstrap.R")



# CARDIIGAN DATA

CAD5=read.dta("C:\\Ben\\Methodological research\\Ordinal prediction CAD\\CADimputed04.dta", convert.factors = FALSE)
cad = CAD5[CAD5$v_mi_m==1,]
cad$crp_2d=1*(cad$crp_2>=0.5)
cad$o5 = ordered(cad$CADextensive+1)

# Descriptive statistics CARDIIGAN
CrossTable(cad$o5)
summary(cad$agexact)
describeBy(cad$agexact)
describeBy(cad$agexact,cad$o5)
summary(cad$hdlchol)
describeBy(cad$hdlchol)
describeBy(cad$hdlchol,cad$o5)
summary(cad$ldlchol)
describeBy(cad$ldlchol)
describeBy(cad$ldlchol,cad$o5)
summary(cad$fibr_ln)
describeBy(exp(cad$fibr_ln))
describeBy(cad$fibr_ln,cad$o5)
describeBy(exp(cad$fibr_ln),cad$o5)
CrossTable(cad$male,cad$o5)
CrossTable(cad$ap_3,cad$o5)
CrossTable(cad$diab,cad$o5)
CrossTable(cad$hypert_m,cad$o5)
CrossTable(cad$dyslip,cad$o5)
CrossTable(cad$smokstat_m,cad$o5)
CrossTable(cad$crp_2d,cad$o5)


######################################################################################
# LIKELIHOOD RATIO TESTS FOR THE PROPORTIONAL ODDS ASSUMPTION IN THE CL-PO FRAMEWORK #
######################################################################################

clpp0 <- vglm(o5 ~ agexact + male + ap_3 + diab + hypert_m + dyslip + smokstat_m + hdlchol + ldlchol + fibr_ln + crp_2d, 
              family=cumulative(parallel=T), data=cad)
clpp <- vglm(o5 ~ agexact + male + ap_3 + diab + hypert_m + dyslip + smokstat_m + hdlchol + ldlchol + fibr_ln + crp_2d, 
             family=cumulative(parallel=F ~ agexact), data=cad)
anova.vglm(clpp,clpp0,type=1)
clpp <- vglm(o5 ~ agexact + male + ap_3 + diab + hypert_m + dyslip + smokstat_m + hdlchol + ldlchol + fibr_ln + crp_2d, 
             family=cumulative(parallel=F ~ male), data=cad)
anova.vglm(clpp,clpp0,type=1)
clpp <- vglm(o5 ~ agexact + male + ap_3 + diab + hypert_m + dyslip + smokstat_m + hdlchol + ldlchol + fibr_ln + crp_2d, 
             family=cumulative(parallel=F ~ ap_3), data=cad)
anova.vglm(clpp,clpp0,type=1)
clpp <- vglm(o5 ~ agexact + male + ap_3 + diab + hypert_m + dyslip + smokstat_m + hdlchol + ldlchol + fibr_ln + crp_2d, 
             family=cumulative(parallel=F ~ diab), data=cad)
anova.vglm(clpp,clpp0,type=1)
clpp <- vglm(o5 ~ agexact + male + ap_3 + diab + hypert_m + dyslip + smokstat_m + hdlchol + ldlchol + fibr_ln + crp_2d, 
             family=cumulative(parallel=F ~ hypert_m), data=cad)
anova.vglm(clpp,clpp0,type=1)
clpp <- vglm(o5 ~ agexact + male + ap_3 + diab + hypert_m + dyslip + smokstat_m + hdlchol + ldlchol + fibr_ln + crp_2d, 
             family=cumulative(parallel=F ~ dyslip), data=cad)
anova.vglm(clpp,clpp0,type=1)
clpp <- vglm(o5 ~ agexact + male + ap_3 + diab + hypert_m + dyslip + smokstat_m + hdlchol + ldlchol + fibr_ln + crp_2d, 
             family=cumulative(parallel=F ~ smokstat_m), data=cad)
anova.vglm(clpp,clpp0,type=1)
clpp <- vglm(o5 ~ agexact + male + ap_3 + diab + hypert_m + dyslip + smokstat_m + hdlchol + ldlchol + fibr_ln + crp_2d, 
             family=cumulative(parallel=F ~ hdlchol), data=cad)
anova.vglm(clpp,clpp0,type=1)
clpp <- vglm(o5 ~ agexact + male + ap_3 + diab + hypert_m + dyslip + smokstat_m + hdlchol + ldlchol + fibr_ln + crp_2d, 
             family=cumulative(parallel=F ~ ldlchol), data=cad)
anova.vglm(clpp,clpp0,type=1)
clpp <- vglm(o5 ~ agexact + male + ap_3 + diab + hypert_m + dyslip + smokstat_m + hdlchol + ldlchol + fibr_ln + crp_2d, 
             family=cumulative(parallel=F ~ fibr_ln), data=cad)
anova.vglm(clpp,clpp0,type=1)
clpp <- vglm(o5 ~ agexact + male + ap_3 + diab + hypert_m + dyslip + smokstat_m + hdlchol + ldlchol + fibr_ln + crp_2d, 
             family=cumulative(parallel=F ~ crp_2d), data=cad)
anova.vglm(clpp,clpp0,type=1)


###################################
# MULTINOMIAL LOGISTIC REGRESSION #
###################################

# Fit the model
mlr <- vglm(o5 ~ agexact + male + ap_3 + diab + hypert_m + dyslip + smokstat_m + hdlchol + ldlchol + fibr_ln + crp_2d, family=multinomial(refLevel = "1"), data=cad)

# Get predicted probabilities, the linear predictors
mlrpred <- predictvglm(mlr,newdata=cad,type="response")
mlrlpred <- predictvglm(mlr,newdata=cad,type="link")

# Calibration intercepts and slopes per outcome category
mlrcalout = calout(out=cad$o5,preds=mlrpred,5)

# Calibration intercepts and slopes per outcome dichotomy
mlrcaldout = caldout(out=cad$o5,preds=mlrpred,k=5)

# Model-specific calibration intercepts and slopes
mlrrecali <- coefficients(vglm(cad$o5 ~ 1, offset = mlrlpred[,1:4], family=multinomial(refLevel = "1")))[c(1:4)]
mlrrecals <- coefficients(vglm(cad$o5 ~ mlrlpred[,1] + mlrlpred[,2] + mlrlpred[,3] + mlrlpred[,4], constraints=list("(Intercept)"=diag(4),"mlrlpred[, 1]"=rbind(1,0,0,0),"mlrlpred[, 2]"=rbind(0,1,0,0),"mlrlpred[, 3]"=rbind(0,0,1,0),"mlrlpred[, 4]"=rbind(0,0,0,1)),family=multinomial(refLevel = "1")))[c(5:8)]

# Flexible recalibration model
mlrlp1=log(mlrpred[,2]/mlrpred[,1])
mlrlp2=log(mlrpred[,3]/mlrpred[,1])
mlrlp3=log(mlrpred[,4]/mlrpred[,1])
mlrlp4=log(mlrpred[,5]/mlrpred[,1])
mlrvgamsmps4 = vgam(cad$o5 ~ sm.ps(mlrlp1,df=4) + sm.ps(c(mlrlp2),df=4) + sm.ps(c(mlrlp3),df=4) + sm.ps(c(mlrlp4),df=4),family=multinomial(refLevel = "1"))

# ECI
mlrECIr = eci_rel(calout=mlrvgamsmps4,preds=mlrpred,k=5,outc=cad$o5)

# Calibration scatter plots and curves per outcome category and per outcome dichotomy
pdf("CARDIIGAN_MLR_Calibrationplots.pdf")
par(mfrow=c(2,2))
plotscatter5(preds=mlrpred,obs=mlrvgamsmps4,cexaxisv=1,cexlabv=1.2,cexlegendv=0.8)
plotlines5(preds=mlrpred,obs=mlrvgamsmps4,cexaxisv=1,cexlabv=1.2,cexlegendv=0.8)
oplotscatter5(preds=mlrpred,obs=mlrvgamsmps4,cexaxisv=1,cexlabv=1.2,cexlegendv=0.8)
oplotlines5(preds=mlrpred,obs=mlrvgamsmps4,cexaxisv=1,cexlabv=1.2,cexlegendv=0.8)
dev.off()

# ordinal c statistic
mlrc = orc(cad$o5,mlrpred,5)

# Combine 'apparent' performance
results_app = c(mlrcalout[,1],mlrcalout[,2],mlrcaldout[,1],mlrcaldout[,2],mlrrecali,mlrrecals,mlrcale,mlrECIr,mlrc)

# Internal validation using bootstrapping
results_corr = bootmlr(seedval=3436,ncat=5,dataset=cad,outc="o5",bootnum=200,appres=results_app,modelformula=o5 ~ agexact + male + ap_3 + diab + hypert_m + dyslip + smokstat_m + hdlchol + ldlchol + fibr_ln + crp_2d)
colnames(results_corr)=c("Apparent","Corrected","Optimism")
rownames(results_corr)=c("ci_cat1","ci_cat2","ci_cat3","ci_cat4","ci_cat5","cs_cat1","cs_cat2","cs_cat3","cs_cat4","cs_cat5",
                         "ci_gt1","ci_gt2","ci_gt3","ci_gt4","cs_gt1","cs_gt2","cs_gt3","cs_gt4",
                         "ci_lp1","ci_lp2","ci_lp3","ci_lp4","cs_lp1","cs_lp2","cs_lp3","cs_lp4","ci_E","cs_E","ECIr","ORC")

save(results_corr,file="CARDIIGAN_MLR.RData")


###########################################
# CUMULATIVE LOGIT WITH PROPORTIONAL ODDS # 
###########################################

# Fit the model
olr <- vglm(o5 ~ agexact + male + ap_3 + diab + hypert_m + dyslip + smokstat_m + hdlchol + ldlchol + fibr_ln + crp_2d, family=cumulative(parallel=T), data=cad)

# Get predicted probabilities, the linear predictors
olrpred <- predictvglm(olr,newdata=cad,type="response")
olrlpred <- predictvglm(olr,newdata=cad,type="link")

# Calibration intercepts and slopes per outcome category
olrcalout = calout(out=cad$o5,preds=olrpred,k=5)

# Calibration intercepts and slopes per outcome dichotomy
olrcaldout = caldout(out=cad$o5,preds=olrpred,k=5)

# Model-specific calibration intercepts and slopes
#olrrecal <- vglm(cad$ordinal3 ~ olrlpred[,1] + olrlpred[,2], family=cumulative(parallel=T)) # NOT OF FULL RANK DUE TO PO, SO DOES NOT WORK
olrrecali <- coefficients(vglm(cad$o5 ~ 1, offset = olrlpred[,1], family=cumulative(parallel=T)))[1]
olrrecals <- coefficients(vglm(cad$o5 ~ olrlpred[,1], family=cumulative(parallel=T)))[c(5)]
olrrecali=c(olrrecali, coefficients(vglm(cad$o5 ~ 1, offset = olrlpred[,2], family=cumulative(parallel=T)))[2])
olrrecals=c(olrrecals, coefficients(vglm(cad$o5 ~ olrlpred[,2], family=cumulative(parallel=T)))[c(5)])
olrrecali=c(olrrecali, coefficients(vglm(cad$o5 ~ 1, offset = olrlpred[,3], family=cumulative(parallel=T)))[3])
olrrecals=c(olrrecals, coefficients(vglm(cad$o5 ~ olrlpred[,3], family=cumulative(parallel=T)))[c(5)])
olrrecali=c(olrrecali, coefficients(vglm(cad$o5 ~ 1, offset = olrlpred[,4], family=cumulative(parallel=T)))[4])
olrrecals=c(olrrecals, coefficients(vglm(cad$o5 ~ olrlpred[,4], family=cumulative(parallel=T)))[c(5)])

# Flexible recalibration model
olrlp1=log(olrpred[,2]/olrpred[,1])
olrlp2=log(olrpred[,3]/olrpred[,1])
olrlp3=log(olrpred[,4]/olrpred[,1])
olrlp4=log(olrpred[,5]/olrpred[,1])
olrvgamsmps4 = vgam(cad$o5 ~ sm.ps(olrlp1,df=4) + sm.ps(olrlp2,df=4) + sm.ps(olrlp3,df=4) + sm.ps(olrlp4,df=4),family=multinomial(refLevel = 1))

# ECI
olrECIr = eci_rel(calout=olrvgamsmps4,preds=olrpred,k=5,outc=cad$o5)

# Calibration scatter plots and curves per outcome category and per outcome dichotomy
pdf("CARDIIGAN_CLPO_Calibrationplots.pdf")
par(mfrow=c(2,2))
plotscatter5(preds=olrpred,obs=olrvgamsmps4,cexaxisv=1,cexlabv=1.2,cexlegendv=0.8)
plotlines5(preds=olrpred,obs=olrvgamsmps4,cexaxisv=1,cexlabv=1.2,cexlegendv=0.8)
oplotscatter5(preds=olrpred,obs=olrvgamsmps4,cexaxisv=1,cexlabv=1.2,cexlegendv=0.8)
oplotlines5(preds=olrpred,obs=olrvgamsmps4,cexaxisv=1,cexlabv=1.2,cexlegendv=0.8)
dev.off()

# ordinal c statistic
olrc = orc(out=cad$o5,preds=olrpred,k=5)

# Combine 'apparent' performance
results_app = c(olrcalout[,1],olrcalout[,2],olrcaldout[,1],olrcaldout[,2],olrrecali,olrrecals,olrcale,olrECIr,olrc)

# Internal validation using bootstrapping
results_corr = bootolr(seedval=3436,ncat=5,dataset=cad,outc="o5",bootnum=200,appres=results_app,modelformula=o5 ~ agexact + male + ap_3 + diab + hypert_m + dyslip + smokstat_m + hdlchol + ldlchol + fibr_ln + crp_2d)
colnames(results_corr)=c("Apparent","Corrected","Optimism")
rownames(results_corr)=c("ci_cat1","ci_cat2","ci_cat3","ci_cat4","ci_cat5","cs_cat1","cs_cat2","cs_cat3","cs_cat4","cs_cat5",
                         "ci_gt1","ci_gt2","ci_gt3","ci_gt4","cs_gt1","cs_gt2","cs_gt3","cs_gt4",
                         "ci_lp1","ci_lp2","ci_lp3","ci_lp4","cs_lp1","cs_lp2","cs_lp3","cs_lp4","ci_E","cs_E","ECIr","ORC")

save(results_corr,file="CARDIIGAN_CLPO.RData")


##################################################
# ADJACENT CATEGORY LOGIT WITH PROPORTIONAL ODDS #
##################################################

# Fit the model
acp <- vglm(o5 ~ agexact + male + ap_3 + diab + hypert_m + dyslip + smokstat_m + hdlchol + ldlchol + fibr_ln + crp_2d, family=acat(parallel=T), data=cad)

# Get predicted probabilities, the linear predictors
acppred <- predictvglm(acp,newdata=cad,type="response")
acplpred <- predictvglm(acp,newdata=cad,type="link")

# Calibration intercepts and slopes per outcome category
acpcalout = calout(out=cad$o5,preds=acppred,k=5)

# Calibration intercepts and slopes per outcome dichotomy
acpcaldout = caldout(out=cad$o5,preds=acppred,k=5)

# Model-specific calibration intercepts and slopes
#acprecal <- vglm(cad$ordinal3 ~ acplpred[,1] + acplpred[,2], family=acat(parallel=T)) # NOT OF FULL RANK DUE TO PO, SO DOES NOT WORK
acprecali <- coefficients(vglm(cad$o5 ~ 1, offset = acplpred[,1], family=acat(parallel=T)))[1]
acprecals <- coefficients(vglm(cad$o5 ~ acplpred[,1], family=acat(parallel=T)))[c(5)]
acprecali=c(acprecali, coefficients(vglm(cad$o5 ~ 1, offset = acplpred[,2], family=acat(parallel=T)))[2])
acprecals=c(acprecals, coefficients(vglm(cad$o5 ~ acplpred[,2], family=acat(parallel=T)))[c(5)])
acprecali=c(acprecali, coefficients(vglm(cad$o5 ~ 1, offset = acplpred[,3], family=acat(parallel=T)))[3])
acprecals=c(acprecals, coefficients(vglm(cad$o5 ~ acplpred[,3], family=acat(parallel=T)))[c(5)])
acprecali=c(acprecali, coefficients(vglm(cad$o5 ~ 1, offset = acplpred[,4], family=acat(parallel=T)))[4])
acprecals=c(acprecals, coefficients(vglm(cad$o5 ~ acplpred[,4], family=acat(parallel=T)))[c(5)])

# Flexible recalibration model
acplp1=log(acppred[,2]/acppred[,1])
acplp2=log(acppred[,3]/acppred[,1])
acplp3=log(acppred[,4]/acppred[,1])
acplp4=log(acppred[,5]/acppred[,1])
acpvgamsmps4 = vgam(cad$o5 ~ sm.ps(acplp1,df=4) + sm.ps(acplp2,df=4) + sm.ps(acplp3,df=4) + sm.ps(acplp4,df=4),family=multinomial(refLevel = 1))

# ECI
acpECIr = eci_rel(calout=acpvgamsmps4,preds=acppred,k=5,outc=cad$o5)

# Calibration plots
pdf("CARDIIGAN_ACPO_Calibrationplots.pdf")
par(mfrow=c(2,2))
plotscatter5(preds=acppred,obs=acpvgamsmps4,cexaxisv=1,cexlabv=1.2,cexlegendv=0.8)
plotlines5(preds=acppred,obs=acpvgamsmps4,cexaxisv=1,cexlabv=1.2,cexlegendv=0.8)
oplotscatter5(preds=acppred,obs=acpvgamsmps4,cexaxisv=1,cexlabv=1.2,cexlegendv=0.8)
oplotlines5(preds=acppred,obs=acpvgamsmps4,cexaxisv=1,cexlabv=1.2,cexlegendv=0.8)
dev.off()

# ordinal c statistic
acpc = orc(out=cad$o5,preds=acppred,k=5)

# Combine 'apparent' performance
results_app = c(acpcalout[,1],acpcalout[,2],acpcaldout[,1],acpcaldout[,2],acprecali,acprecals,acpcale,acpECIr,acpc)

# Internal validation using bootstrapping
results_corr = bootacp(seedval=3436,ncat=5,dataset=cad,outc="o5",bootnum=200,appres=results_app,modelformula=o5 ~ agexact + male + ap_3 + diab + hypert_m + dyslip + smokstat_m + hdlchol + ldlchol + fibr_ln + crp_2d)
colnames(results_corr)=c("Apparent","Corrected","Optimism")
rownames(results_corr)=c("ci_cat1","ci_cat2","ci_cat3","ci_cat4","ci_cat5","cs_cat1","cs_cat2","cs_cat3","cs_cat4","cs_cat5",
                         "ci_gt1","ci_gt2","ci_gt3","ci_gt4","cs_gt1","cs_gt2","cs_gt3","cs_gt4",
                         "ci_lp1","ci_lp2","ci_lp3","ci_lp4","cs_lp1","cs_lp2","cs_lp3","cs_lp4","ci_E","cs_E","ECIr","ORC")

save(results_corr,file="CARDIIGAN_ACPO.RData")


###################################################
# CONTINUATION RATIO LOGIT WITH PROPORTIONAL ODDS #
###################################################

# Fit the model
crp <- vglm(o5 ~ agexact + male + ap_3 + diab + hypert_m + dyslip + smokstat_m + hdlchol + ldlchol + fibr_ln + crp_2d, family=cratio(parallel=T), data=cad)

# Get predicted probabilities, the linear predictors
crppred <- predictvglm(crp,newdata=cad,type="response")
crplpred <- predictvglm(crp,newdata=cad,type="link")

# Calibration intercepts and slopes per outcome category
crpcalout = calout(out=cad$o5,preds=crppred,k=5)

# Calibration intercepts and slopes per outcome dichotomy
crpcaldout = caldout(out=cad$o5,preds=crppred,k=5)

# Model-specific calibration intercepts and slopes
#crprecal <- vglm(cad$o5 ~ crplpred[,1] + crplpred[,2], family=cratio(parallel=T)) # NOT OF FULL RANK  DUE TO PO, SO DOES NOT WORK
crprecali <- coefficients(vglm(cad$o5 ~ 1, offset = crplpred[,1], family=cratio(parallel=T)))[1]
crprecals <- coefficients(vglm(cad$o5 ~ crplpred[,1], family=cratio(parallel=T)))[c(5)]
crprecali=c(crprecali, coefficients(vglm(cad$o5 ~ 1, offset = crplpred[,2], family=cratio(parallel=T)))[2])
crprecals=c(crprecals, coefficients(vglm(cad$o5 ~ crplpred[,2], family=cratio(parallel=T)))[c(5)])
crprecali=c(crprecali, coefficients(vglm(cad$o5 ~ 1, offset = crplpred[,3], family=cratio(parallel=T)))[3])
crprecals=c(crprecals, coefficients(vglm(cad$o5 ~ crplpred[,3], family=cratio(parallel=T)))[c(5)])
crprecali=c(crprecali, coefficients(vglm(cad$o5 ~ 1, offset = crplpred[,4], family=cratio(parallel=T)))[4])
crprecals=c(crprecals, coefficients(vglm(cad$o5 ~ crplpred[,4], family=cratio(parallel=T)))[c(5)])

# Flexible recalibration model
crplp1=log(crppred[,2]/crppred[,1])
crplp2=log(crppred[,3]/crppred[,1])
crplp3=log(crppred[,4]/crppred[,1])
crplp4=log(crppred[,5]/crppred[,1])
crpvgamsmps4 = vgam(cad$o5 ~ sm.ps(crplp1,df=4) + sm.ps(crplp2,df=4) + sm.ps(crplp3,df=4) + sm.ps(crplp4,df=4),family=multinomial(refLevel = 1))

# ECI
crpECIr = eci_rel(calout=crpvgamsmps4,preds=crppred,k=5,outc=cad$o5)

# Calibration plots
pdf("CARDIIGAN_CRPO_Calibrationplots.pdf")
par(mfrow=c(2,2))
plotscatter5(preds=crppred,obs=crpvgamsmps4,cexaxisv=1,cexlabv=1.2,cexlegendv=0.8)
plotlines5(preds=crppred,obs=crpvgamsmps4,cexaxisv=1,cexlabv=1.2,cexlegendv=0.8)
oplotscatter5(preds=crppred,obs=crpvgamsmps4,cexaxisv=1,cexlabv=1.2,cexlegendv=0.8)
oplotlines5(preds=crppred,obs=crpvgamsmps4,cexaxisv=1,cexlabv=1.2,cexlegendv=0.8)
dev.off()

# ordinal c statistic
crpc = orc(out=cad$o5,preds=crppred,k=5)

# Combine 'apparent' performance
results_app = c(crpcalout[,1],crpcalout[,2],crpcaldout[,1],crpcaldout[,2],crprecali,crprecals,crpcale,crpECIr,crpc)

# Internal validation using bootstrapping
results_corr = bootcrp(seedval=3436,ncat=5,dataset=cad,outc="o5",bootnum=200,appres=results_app,modelformula=o5 ~ agexact + male + ap_3 + diab + hypert_m + dyslip + smokstat_m + hdlchol + ldlchol + fibr_ln + crp_2d)
colnames(results_corr)=c("Apparent","Corrected","Optimism")
rownames(results_corr)=c("ci_cat1","ci_cat2","ci_cat3","ci_cat4","ci_cat5","cs_cat1","cs_cat2","cs_cat3","cs_cat4","cs_cat5",
                         "ci_gt1","ci_gt2","ci_gt3","ci_gt4","cs_gt1","cs_gt2","cs_gt3","cs_gt4",
                         "ci_lp1","ci_lp2","ci_lp3","ci_lp4","cs_lp1","cs_lp2","cs_lp3","cs_lp4","ci_E","cs_E","ECIr","ORC")

save(results_corr,file="CARDIIGAN_CRPO.RData")


######################################################
# CONTINUATION RATIO LOGIT WITHOUT PROPORTIONAL ODDS #
######################################################

# Fit the model
crnp <- vglm(o5 ~ agexact + male + ap_3 + diab + hypert_m + dyslip + smokstat_m + hdlchol + ldlchol + fibr_ln + crp_2d, family=cratio(parallel=F), data=cad)

# Get predicted probabilities, the linear predictors
crnppred <- predictvglm(crnp,newdata=cad,type="response")
crnplpred <- predictvglm(crnp,newdata=cad,type="link")

# Calibration intercepts and slopes per outcome category
crnpcalout = calout(out=cad$o5,preds=crnppred,k=5)

# Calibration intercepts and slopes per outcome dichotomy
crnpcaldout = caldout(out=cad$o5,preds=crnppred,k=5)

# Model-specific calibration intercepts and slopes
crnprecali <- coefficients(vglm(cad$o5 ~ 1, offset = crnplpred[,1:4], family=cratio(parallel=F)))[c(1:4)]
crnprecals <- coefficients(vglm(cad$o5 ~ crnplpred[,1] + crnplpred[,2] + crnplpred[,3] + crnplpred[,4], constraints=list("(Intercept)"=diag(4),"crnplpred[, 1]"=rbind(1,0,0,0),"crnplpred[, 2]"=rbind(0,1,0,0),"crnplpred[, 3]"=rbind(0,0,1,0),"crnplpred[, 4]"=rbind(0,0,0,1)),family=cratio(parallel=F)))[c(5:8)]

# Flexible recalibration model
crnplp1=log(crnppred[,2]/crnppred[,1])
crnplp2=log(crnppred[,3]/crnppred[,1])
crnplp3=log(crnppred[,4]/crnppred[,1])
crnplp4=log(crnppred[,5]/crnppred[,1])
crnpvgamsmps4 = vgam(cad$o5 ~ sm.ps(crnplp1,df=4) + sm.ps(crnplp2,df=4) + sm.ps(crnplp3,df=4) + sm.ps(crnplp4,df=4),family=multinomial(refLevel = 1))

# ECI
crnpECIr = eci_rel(calout=crnpvgamsmps4,preds=crnppred,k=5,outc=cad$o5)

# Calibration plots
pdf("CARDIIGAN_CRNP_Calibrationplots.pdf")
par(mfrow=c(2,2))
plotscatter5(preds=crnppred,obs=crnpvgamsmps4,cexaxisv=1,cexlabv=1.2,cexlegendv=0.8)
plotlines5(preds=crnppred,obs=crnpvgamsmps4,cexaxisv=1,cexlabv=1.2,cexlegendv=0.8)
oplotscatter5(preds=crnppred,obs=crnpvgamsmps4,cexaxisv=1,cexlabv=1.2,cexlegendv=0.8)
oplotlines5(preds=crnppred,obs=crnpvgamsmps4,cexaxisv=1,cexlabv=1.2,cexlegendv=0.8)
dev.off()

# ordinal c statistic
crnpc = orc(out=cad$o5,preds=crnppred,k=5)

# Combine 'apparent' performance
results_app = c(crnpcalout[,1],crnpcalout[,2],crnpcaldout[,1],crnpcaldout[,2],crnprecali,crnprecals,crnpcale,crnpECIr,crnpc)

# Internal validation using bootstrapping
results_corr = bootcrnp(seedval=3436,ncat=5,dataset=cad,outc="o5",bootnum=200,appres=results_app,modelformula=o5 ~ agexact + male + ap_3 + diab + hypert_m + dyslip + smokstat_m + hdlchol + ldlchol + fibr_ln + crp_2d)
colnames(results_corr)=c("Apparent","Corrected","Optimism")
rownames(results_corr)=c("ci_cat1","ci_cat2","ci_cat3","ci_cat4","ci_cat5","cs_cat1","cs_cat2","cs_cat3","cs_cat4","cs_cat5",
                         "ci_gt1","ci_gt2","ci_gt3","ci_gt4","cs_gt1","cs_gt2","cs_gt3","cs_gt4",
                         "ci_lp1","ci_lp2","ci_lp3","ci_lp4","cs_lp1","cs_lp2","cs_lp3","cs_lp4","ci_E","cs_E","ECIr","ORC")

save(results_corr,file="CARDIIGAN_CRNP.RData")


###############################################
# STEREOTYPE MULTINOMIAL LOGISTIC REGRESSSION #
###############################################

# Fit the model
sm=rrvglm(o5 ~ agexact + male + ap_3 + diab + hypert_m + dyslip + smokstat_m + hdlchol + ldlchol + fibr_ln + crp_2d, multinomial(refLevel = "1"), data = cad)

# Get predicted probabilities, the linear predictors, and the expected outcome E
smpred <- predictvglm(sm,newdata=cad,type="response")
smlpred <- predictvglm(sm,newdata=cad,type="link")

# Calibration intercepts and slopes per outcome category
smcalout = calout(out=cad$o5,preds=smpred,k=5)

# Calibration intercepts and slopes per outcome dichotomy
smcaldout = caldout(out=cad$o5,preds=smpred,k=5)

# Model-specific calibration intercepts and slopes
#smrecal <- rrvglm(cad$o5 ~ smlpred[,1] + smlpred[,2], multinomial(refLevel = 1)) # NOT OF FULL RANK DUE TO PO, SO DOES NOT WORK 
smrecali <- coefficients(vglm(cad$o5 ~ 1, offset = smlpred[,1:4], family=multinomial(refLevel = "1")))[c(1:4)]
smrecals <- coefficients(vglm(cad$o5 ~ smlpred[,1] + smlpred[,2] + smlpred[,3] + smlpred[,4], constraints=list("(Intercept)"=diag(4),"smlpred[, 1]"=rbind(1,0,0,0),"smlpred[, 2]"=rbind(0,1,0,0),"smlpred[, 3]"=rbind(0,0,1,0),"smlpred[, 4]"=rbind(0,0,0,1)),family=multinomial(refLevel = "1")))[c(5:8)]

# Flexible recalibration model
smlp1=log(smpred[,2]/smpred[,1])
smlp2=log(smpred[,3]/smpred[,1])
smlp3=log(smpred[,4]/smpred[,1])
smlp4=log(smpred[,5]/smpred[,1])
smvgamsmps4 = vgam(cad$o5 ~ sm.ps(smlp1,df=4) + sm.ps(smlp2,df=4) + sm.ps(smlp3,df=4) + sm.ps(smlp4,df=4),family=multinomial(refLevel = 1))

# ECI
smECIr = eci_rel(calout=smvgamsmps4,preds=smpred,k=5,outc=cad$o5)

# Calibration plots
pdf("CARDIIGAN_SM_Calibrationplots.pdf")
par(mfrow=c(2,2))
plotscatter5(preds=smpred,obs=smvgamsmps4,cexaxisv=1,cexlabv=1.2,cexlegendv=0.8)
plotlines5(preds=smpred,obs=smvgamsmps4,cexaxisv=1,cexlabv=1.2,cexlegendv=0.8)
oplotscatter5(preds=smpred,obs=smvgamsmps4,cexaxisv=1,cexlabv=1.2,cexlegendv=0.8)
oplotlines5(preds=smpred,obs=smvgamsmps4,cexaxisv=1,cexlabv=1.2,cexlegendv=0.8)
dev.off()

# ordinal c statistic
smc = orc(out=cad$o5,preds=smpred,k=5)

# Combine 'apparent' performance
results_app = c(smcalout[,1],smcalout[,2],smcaldout[,1],smcaldout[,2],smrecali,smrecals,smcale,smECIr,smc)

# Internal validation using bootstrapping
results_corr = bootsm(seedval=3436,ncat=5,dataset=cad,outc="o5",bootnum=200,appres=results_app,modelformula=o5 ~ agexact + male + ap_3 + diab + hypert_m + dyslip + smokstat_m + hdlchol + ldlchol + fibr_ln + crp_2d)
colnames(results_corr)=c("Apparent","Corrected","Optimism")
rownames(results_corr)=c("ci_cat1","ci_cat2","ci_cat3","ci_cat4","ci_cat5","cs_cat1","cs_cat2","cs_cat3","cs_cat4","cs_cat5",
                         "ci_gt1","ci_gt2","ci_gt3","ci_gt4","cs_gt1","cs_gt2","cs_gt3","cs_gt4",
                         "ci_lp1","ci_lp2","ci_lp3","ci_lp4","cs_lp1","cs_lp2","cs_lp3","cs_lp4","ci_E","cs_E","ECIr","ORC")

save(results_corr,file="CARDIIGAN_SM.RData")


#############################################
# Descriptive statistics of estimated risks #
#############################################

summary(mlrpred)
summary(olrpred)
summary(acppred)
summary(crppred)
summary(crnppred)
summary(smpred)


#######################################################
# SCATTER PLOTS OF PREDICTED RISKS BETWEEN ALL MODELS #
#######################################################

ref=cbind(as.matrix(seq(0,1,by=0.01)),as.matrix(seq(0,1,by=0.01)),as.matrix(seq(0,1,by=0.01)),as.matrix(seq(0,1,by=0.01)),as.matrix(seq(0,1,by=0.01)))
mlrpred2=rbind(ref,mlrpred)
olrpred2=rbind(ref,olrpred)
acppred2=rbind(ref,acppred)
crppred2=rbind(ref,crppred)
crnppred2=rbind(ref,crnppred)
smpred2=rbind(ref,smpred)
type=c(rep("red",length(seq(0,1,by=0.01))),rep("black",length(mlrpred)))

pairs(~mlrpred2[,1]+olrpred2[,1]+acppred2[,1]+crppred2[,1]+crnppred2[,1]+smpred2[,1],col=type,upper.panel=NULL,xlim=c(0,1),ylim=c(0,1),cex=0.1,labels=c("MLR","CL-PO","AC-PO","CR-PO","CR-NP","SLM"),main=NULL)
pairs(~mlrpred2[,2]+olrpred2[,2]+acppred2[,2]+crppred2[,2]+crnppred2[,2]+smpred2[,2],col=type,upper.panel=NULL,xlim=c(0,1),ylim=c(0,1),cex=0.1,labels=c("MLR","CL-PO","AC-PO","CR-PO","CR-NP","SLM"),main=NULL)
pairs(~mlrpred2[,3]+olrpred2[,3]+acppred2[,3]+crppred2[,3]+crnppred2[,3]+smpred2[,3],col=type,upper.panel=NULL,xlim=c(0,1),ylim=c(0,1),cex=0.1,labels=c("MLR","CL-PO","AC-PO","CR-PO","CR-NP","SLM"),main=NULL)
pairs(~mlrpred2[,4]+olrpred2[,4]+acppred2[,4]+crppred2[,4]+crnppred2[,4]+smpred2[,4],col=type,upper.panel=NULL,xlim=c(0,1),ylim=c(0,1),cex=0.1,labels=c("MLR","CL-PO","AC-PO","CR-PO","CR-NP","SLM"),main=NULL)
pairs(~mlrpred2[,5]+olrpred2[,5]+acppred2[,5]+crppred2[,5]+crnppred2[,5]+smpred2[,5],col=type,upper.panel=NULL,xlim=c(0,1),ylim=c(0,1),cex=0.1,labels=c("MLR","CL-PO","AC-PO","CR-PO","CR-NP","SLM"),main=NULL)


############################################
# Calibration scatter plots for all models #
############################################

# FIGURE IN MAIN PAPER

par(mfrow=c(3,2),
    oma=c(5,5,0,0) + 0.0,
    mar=c(0,0,1,1) + 0.0)
ref <- rbind(c(0,0),c(1,1)) # needed to plot the ideal diagonal line
preds=mlrpred
obs=mlrvgamsmps4
plot(preds[,1],fitted(obs)[,1],type="p",pch=1,col="green",lwd=1,xlim=0:1,ylim=0:1,cex=0.8,axes=F)
axis(side=2,at=c(0,0.2,0.4,0.6,0.8,1),cex.axis=1.2)
box(which="plot")
points(preds[,2],fitted(obs)[,2],type="p",pch=1,col="orange",lwd=1,cex=0.8)
points(preds[,3],fitted(obs)[,3],type="p",pch=1,col="red",lwd=1,cex=0.8)
points(preds[,4],fitted(obs)[,4],type="p",pch=1,col="brown",lwd=1,cex=0.8)
points(preds[,5],fitted(obs)[,5],type="p",pch=1,col="black",lwd=1,cex=0.8)
lines(ref,ref,type="l",col="gray",lty=3,lwd=2) # plot the ideal diagonal line
text(1, 0, "MLR", adj=c(1,0),cex=1.5)

preds=olrpred
obs=olrvgamsmps4
plot(preds[,1],fitted(obs)[,1],type="p",pch=1,col="green",lwd=1,xlim=0:1,ylim=0:1,cex=0.8,axes=F)
box(which="plot")
points(preds[,2],fitted(obs)[,2],type="p",pch=1,col="orange",lwd=1,cex=0.8)
points(preds[,3],fitted(obs)[,3],type="p",pch=1,col="red",lwd=1,cex=0.8)
points(preds[,4],fitted(obs)[,4],type="p",pch=1,col="brown",lwd=1,cex=0.8)
points(preds[,5],fitted(obs)[,5],type="p",pch=1,col="black",lwd=1,cex=0.8)
lines(ref,ref,type="l",col="gray",lty=3,lwd=2) # plot the ideal diagonal line
text(1, 0, "CL-PO", adj=c(1,0),cex=1.5)

preds=acppred
obs=acpvgamsmps4
plot(preds[,1],fitted(obs)[,1],type="p",pch=1,col="green",lwd=1,xlim=0:1,ylim=0:1,cex=0.8,axes=F)
axis(side=2,at=c(0,0.2,0.4,0.6,0.8,1),cex.axis=1.2)
box(which="plot")
points(preds[,2],fitted(obs)[,2],type="p",pch=1,col="orange",lwd=1,cex=0.8)
points(preds[,3],fitted(obs)[,3],type="p",pch=1,col="red",lwd=1,cex=0.8)
points(preds[,4],fitted(obs)[,4],type="p",pch=1,col="brown",lwd=1,cex=0.8)
points(preds[,5],fitted(obs)[,5],type="p",pch=1,col="black",lwd=1,cex=0.8)
lines(ref,ref,type="l",col="gray",lty=3,lwd=2) # plot the ideal diagonal line
text(1, 0, "AC-PO", adj=c(1,0),cex=1.5)

preds=crppred
obs=crpvgamsmps4
plot(preds[,1],fitted(obs)[,1],type="p",pch=1,col="green",lwd=1,xlim=0:1,ylim=0:1,cex=0.8,axes=F)
box(which="plot")
points(preds[,2],fitted(obs)[,2],type="p",pch=1,col="orange",lwd=1,cex=0.8)
points(preds[,3],fitted(obs)[,3],type="p",pch=1,col="red",lwd=1,cex=0.8)
points(preds[,4],fitted(obs)[,4],type="p",pch=1,col="brown",lwd=1,cex=0.8)
points(preds[,5],fitted(obs)[,5],type="p",pch=1,col="black",lwd=1,cex=0.8)
lines(ref,ref,type="l",col="gray",lty=3,lwd=2) # plot the ideal diagonal line
text(1, 0, "CR-PO", adj=c(1,0),cex=1.5)

preds=crnppred
obs=crnpvgamsmps4
plot(preds[,1],fitted(obs)[,1],type="p",pch=1,col="green",lwd=1,xlim=0:1,ylim=0:1,cex=0.8,axes=F)
axis(side=1,at=c(0,0.2,0.4,0.6,0.8,1),cex.axis=1.2)
axis(side=2,at=c(0,0.2,0.4,0.6,0.8,1),cex.axis=1.2)
box(which="plot")
points(preds[,2],fitted(obs)[,2],type="p",pch=1,col="orange",lwd=1,cex=0.8)
points(preds[,3],fitted(obs)[,3],type="p",pch=1,col="red",lwd=1,cex=0.8)
points(preds[,4],fitted(obs)[,4],type="p",pch=1,col="brown",lwd=1,cex=0.8)
points(preds[,5],fitted(obs)[,5],type="p",pch=1,col="black",lwd=1,cex=0.8)
lines(ref,ref,type="l",col="gray",lty=3,lwd=2) # plot the ideal diagonal line
text(1, 0, "CR-NP", adj=c(1,0),cex=1.5)

preds=smpred
obs=smvgamsmps4
plot(preds[,1],fitted(obs)[,1],type="p",pch=1,col="green",lwd=1,xlim=0:1,ylim=0:1,cex=0.8,axes=F)
axis(side=1,at=c(0,0.2,0.4,0.6,0.8,1),cex.axis=1.2)
box(which="plot")
points(preds[,2],fitted(obs)[,2],type="p",pch=1,col="orange",lwd=1,cex=0.8)
points(preds[,3],fitted(obs)[,3],type="p",pch=1,col="red",lwd=1,cex=0.8)
points(preds[,4],fitted(obs)[,4],type="p",pch=1,col="brown",lwd=1,cex=0.8)
points(preds[,5],fitted(obs)[,5],type="p",pch=1,col="black",lwd=1,cex=0.8)
lines(ref,ref,type="l",col="gray",lty=3,lwd=2) # plot the ideal diagonal line
text(1, 0, "SLM", adj=c(1,0),cex=1.5)

title(xlab = list("Estimated probability",cex=2),
      ylab = list("Observed proportion",cex=2),
      outer = TRUE, line=3)


####################################################
# Ordinal calibration scatter plots for all models #
####################################################

# FIGURE IN MAIN PAPER

par(mfrow=c(3,2),
    oma=c(5,5,0,0) + 0.0,
    mar=c(0,0,1,1) + 0.0)
ref <- rbind(c(0,0),c(1,1)) # needed to plot the ideal diagonal line
preds=mlrpred
obs=mlrvgamsmps4
plot(preds[,2]+preds[,3]+preds[,4]+preds[,5],fitted(obs)[,2]+fitted(obs)[,3]+fitted(obs)[,4]+fitted(obs)[,5],type="p",pch=1,col="orange",lwd=1,xlim=0:1,ylim=0:1,cex=0.8,axes=F)
axis(side=2,at=c(0,0.2,0.4,0.6,0.8,1),cex.axis=1.2)
box(which="plot")
points(preds[,3]+preds[,4]+preds[,5],fitted(obs)[,3]+fitted(obs)[,4]+fitted(obs)[,5],type="p",pch=1,col="red",lwd=1,cex=0.8)
points(preds[,4]+preds[,5],fitted(obs)[,4]+fitted(obs)[,5],type="p",pch=1,col="brown",lwd=1,cex=0.8)
points(preds[,5],fitted(obs)[,5],type="p",pch=1,col="black",lwd=1,cex=0.8)
lines(ref,ref,type="l",col="gray",lty=3,lwd=2) # plot the ideal diagonal line
text(1, 0, "MLR", adj=c(1,0),cex=1.5)

preds=olrpred
obs=olrvgamsmps4
plot(preds[,2]+preds[,3]+preds[,4]+preds[,5],fitted(obs)[,2]+fitted(obs)[,3]+fitted(obs)[,4]+fitted(obs)[,5],type="p",pch=1,col="orange",lwd=1,xlim=0:1,ylim=0:1,cex=0.8,axes=F)
box(which="plot")
points(preds[,3]+preds[,4]+preds[,5],fitted(obs)[,3]+fitted(obs)[,4]+fitted(obs)[,5],type="p",pch=1,col="red",lwd=1,cex=0.8)
points(preds[,4]+preds[,5],fitted(obs)[,4]+fitted(obs)[,5],type="p",pch=1,col="brown",lwd=1,cex=0.8)
points(preds[,5],fitted(obs)[,5],type="p",pch=1,col="black",lwd=1,cex=0.8)
lines(ref,ref,type="l",col="gray",lty=3,lwd=2) # plot the ideal diagonal line
text(1, 0, "CL-PO", adj=c(1,0),cex=1.5)

preds=acppred
obs=acpvgamsmps4
plot(preds[,2]+preds[,3]+preds[,4]+preds[,5],fitted(obs)[,2]+fitted(obs)[,3]+fitted(obs)[,4]+fitted(obs)[,5],type="p",pch=1,col="orange",lwd=1,xlim=0:1,ylim=0:1,cex=0.8,axes=F)
axis(side=2,at=c(0,0.2,0.4,0.6,0.8,1),cex.axis=1.2)
box(which="plot")
points(preds[,3]+preds[,4]+preds[,5],fitted(obs)[,3]+fitted(obs)[,4]+fitted(obs)[,5],type="p",pch=1,col="red",lwd=1,cex=0.8)
points(preds[,4]+preds[,5],fitted(obs)[,4]+fitted(obs)[,5],type="p",pch=1,col="brown",lwd=1,cex=0.8)
points(preds[,5],fitted(obs)[,5],type="p",pch=1,col="black",lwd=1,cex=0.8)
lines(ref,ref,type="l",col="gray",lty=3,lwd=2) # plot the ideal diagonal line
text(1, 0, "AC-PO", adj=c(1,0),cex=1.5)

preds=crppred
obs=crpvgamsmps4
plot(preds[,2]+preds[,3]+preds[,4]+preds[,5],fitted(obs)[,2]+fitted(obs)[,3]+fitted(obs)[,4]+fitted(obs)[,5],type="p",pch=1,col="orange",lwd=1,xlim=0:1,ylim=0:1,cex=0.8,axes=F)
box(which="plot")
points(preds[,3]+preds[,4]+preds[,5],fitted(obs)[,3]+fitted(obs)[,4]+fitted(obs)[,5],type="p",pch=1,col="red",lwd=1,cex=0.8)
points(preds[,4]+preds[,5],fitted(obs)[,4]+fitted(obs)[,5],type="p",pch=1,col="brown",lwd=1,cex=0.8)
points(preds[,5],fitted(obs)[,5],type="p",pch=1,col="black",lwd=1,cex=0.8)
lines(ref,ref,type="l",col="gray",lty=3,lwd=2) # plot the ideal diagonal line
text(1, 0, "CR-PO", adj=c(1,0),cex=1.5)

preds=crnppred
obs=crnpvgamsmps4
plot(preds[,2]+preds[,3]+preds[,4]+preds[,5],fitted(obs)[,2]+fitted(obs)[,3]+fitted(obs)[,4]+fitted(obs)[,5],type="p",pch=1,col="orange",lwd=1,xlim=0:1,ylim=0:1,cex=0.8,axes=F)
axis(side=1,at=c(0,0.2,0.4,0.6,0.8,1),cex.axis=1.2)
axis(side=2,at=c(0,0.2,0.4,0.6,0.8,1),cex.axis=1.2)
box(which="plot")
points(preds[,3]+preds[,4]+preds[,5],fitted(obs)[,3]+fitted(obs)[,4]+fitted(obs)[,5],type="p",pch=1,col="red",lwd=1,cex=0.8)
points(preds[,4]+preds[,5],fitted(obs)[,4]+fitted(obs)[,5],type="p",pch=1,col="brown",lwd=1,cex=0.8)
points(preds[,5],fitted(obs)[,5],type="p",pch=1,col="black",lwd=1,cex=0.8)
lines(ref,ref,type="l",col="gray",lty=3,lwd=2) # plot the ideal diagonal line
text(1, 0, "CR-NP", adj=c(1,0),cex=1.5)

preds=smpred
obs=smvgamsmps4
plot(preds[,2]+preds[,3]+preds[,4]+preds[,5],fitted(obs)[,2]+fitted(obs)[,3]+fitted(obs)[,4]+fitted(obs)[,5],type="p",pch=1,col="orange",lwd=1,xlim=0:1,ylim=0:1,cex=0.8,axes=F)
axis(side=1,at=c(0,0.2,0.4,0.6,0.8,1),cex.axis=1.2)
box(which="plot")
points(preds[,3]+preds[,4]+preds[,5],fitted(obs)[,3]+fitted(obs)[,4]+fitted(obs)[,5],type="p",pch=1,col="red",lwd=1,cex=0.8)
points(preds[,4]+preds[,5],fitted(obs)[,4]+fitted(obs)[,5],type="p",pch=1,col="brown",lwd=1,cex=0.8)
points(preds[,5],fitted(obs)[,5],type="p",pch=1,col="black",lwd=1,cex=0.8)
lines(ref,ref,type="l",col="gray",lty=3,lwd=2) # plot the ideal diagonal line
text(1, 0, "SLM", adj=c(1,0),cex=1.5)

title(xlab = list("Estimated probability",cex=2),
      ylab = list("Observed proportion",cex=2),
      outer = TRUE, line=3)


###############################
# LOAD RESULTS FOR EACH MODEL #
###############################

load("CARDIIGAN case study\\CARDIIGAN_MLR.rdata")
load("CARDIIGAN case study\\CARDIIGAN_CLPO.rdata")
load("CARDIIGAN case study\\CARDIIGAN_ACPO.rdata")
load("CARDIIGAN case study\\CARDIIGAN_CRPO.rdata")
load("CARDIIGAN case study\\CARDIIGAN_CRNP.rdata")
load("CARDIIGAN case study\\CARDIIGAN_SM.rdata")
round(results_corr,digits=3)
