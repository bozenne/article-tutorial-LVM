### multivariate-regression.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar 28 2022 (10:26) 
## Version: 
## Last-Updated: mar 31 2022 (16:36) 
##           By: Brice Ozenne
##     Update #: 7
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

library(lava)
library(lavaSearch2)
library(lme4)
library(nlme)
library(multcomp)
library(LMMstar)

## * load data
dfW.SPECT <- read.table("data/data-SPECT.txt")
dfW.SPECT$group <- factor(dfW.SPECT$group, c("healthy","concussion"))
dfW.SPECT$genotype <- factor(dfW.SPECT$genotype, c("MAB","HAB"))

name.Y <- c("log.thalamus", 
            "log.pallidostriatum", 
            "log.neocortex", 
            "log.midbrain", 
            "log.pons", 
            "log.cingulateGyrus", 
            "log.hippocampus", 
            "log.supramarginalGyrus", 
            "log.corpusCallosum")
dfL.SPECT <- reshape(dfW.SPECT, idvar = "id", direction = "long",
                     varying = name.Y)
names(dfL.SPECT)[names(dfL.SPECT)=="time"] <- "region"
names(dfL.SPECT)[names(dfL.SPECT)=="log"] <- "SPECT"
rownames(dfL.SPECT) <- NULL

## * multiple linear regression
e.mlm <- lm(cbind(log.thalamus, 
                  log.pallidostriatum, 
                  log.neocortex, 
                  log.midbrain, 
                  log.pons, 
                  log.cingulateGyrus, 
                  log.hippocampus, 
                  log.supramarginalGyrus, 
                  log.corpusCallosum)~genotype+group, data = dfW.SPECT)

e.mlm <- list(thalamus = lm(log.thalamus~genotype+group, data = dfW.SPECT),
              pallidostriatum = lm(log.pallidostriatum~genotype+group, data = dfW.SPECT),
              neocortex = lm(log.neocortex~genotype+group, data = dfW.SPECT),
              midbrain = lm(log.midbrain~genotype+group, data = dfW.SPECT),
              pons = lm(log.pons~genotype+group, data = dfW.SPECT),
              cingulate.gyrus = lm(log.cingulateGyrus~genotype+group, data = dfW.SPECT),
              hippocampus = lm(log.hippocampus~genotype+group, data = dfW.SPECT),
              supramarginal.gyrus = lm(log.supramarginalGyrus~genotype+group, data = dfW.SPECT),
              corpus.callosum = lm(log.corpusCallosum~genotype+group, data = dfW.SPECT))
class(e.mlm) <- "mmm"

e.glhtmlm <- glht(e.mlm, linfct = mlf("groupconcussion=0"))
summary(e.glhtmlm, test = adjusted("bonferroni"))
summary(e.glhtmlm)
## Linear Hypotheses:
##                                           Estimate Std. Error z value Pr(>|z|)    
## thalamus: groupconcussion == 0             0.13061    0.04514   2.894  0.01788 *  
## pallidostriatum: groupconcussion == 0      0.13536    0.04175   3.242  0.00613 ** 
## neocortex: groupconcussion == 0            0.05692    0.03985   1.428  0.43082    
## midbrain: groupconcussion == 0             0.13831    0.03976   3.479  0.00251 ** 
## pons: groupconcussion == 0                 0.03494    0.04473   0.781  0.88475    
## cingulate.gyrus: groupconcussion == 0      0.20342    0.03500   5.813  < 0.001 ***
## hippocampus: groupconcussion == 0          0.15788    0.03765   4.193  < 0.001 ***
## supramarginal.gyrus: groupconcussion == 0  0.06198    0.03883   1.596  0.33216    
## corpus.callosum: groupconcussion == 0      0.16897    0.03569   4.734  < 0.001 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## (Adjusted p values reported -- single-step method)

## * Random intercept model
e.CS <- lmer(SPECT ~ region*group + region*genotype + (1|id), data = dfL.SPECT)
logLik(e.CS)


e.CS <- lmer(SPECT ~ region/group + region/genotype + (1|id), data = dfL.SPECT)
summary(e.CS)

e.glhtCS <- glht(e.CS, linfct = c("regionthalamus:groupconcussion=0",
                                  "regionpallidostriatum:groupconcussion=0",
                                  "regionneocortex:groupconcussion=0",
                                  "regionmidbrain:groupconcussion=0",
                                  "regionpons:groupconcussion=0",
                                  "regioncingulateGyrus:groupconcussion=0",
                                  "regionhippocampus:groupconcussion=0",
                                  "regionsupramarginalGyrus:groupconcussion=0",
                                  "regioncorpusCallosum:groupconcussion=0"
                                  ))
rownames(e.glhtCS$linfct) <- c("thalamus", "pallidostriatum", "neocortex", "midbrain", "pons",
                               "cingulate.gyrus", "hippocampus", "supramarginal.gyrus", "corpus.callosum")    

summary(e.glhtCS)
## Linear Hypotheses:
##                          Estimate Std. Error z value Pr(>|z|)    
## thalamus == 0             0.13061    0.03996   3.268  0.00584 ** 
## pallidostriatum == 0      0.13536    0.03996   3.387  0.00396 ** 
## neocortex == 0            0.05692    0.03996   1.424  0.45066    
## midbrain == 0             0.13831    0.03996   3.461  0.00337 ** 
## pons == 0                 0.03494    0.03996   0.874  0.84914    
## cingulate.gyrus == 0      0.20342    0.03996   5.090  < 0.001 ***
## hippocampus == 0          0.15788    0.03996   3.950  < 0.001 ***
## supramarginal.gyrus == 0  0.06198    0.03996   1.551  0.37250    
## corpus.callosum == 0      0.16897    0.03996   4.228  < 0.001 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## (Adjusted p values reported -- single-step method)

## * Mixed model with unstructured covariance matrix
e.UN <- gls(SPECT ~ region/group + region/genotype,
            correlation = corSymm(form = ~1|id),
            weights = varIdent(form = ~1|region),
            data = dfL.SPECT)

e.glhtUN <- glht(e.UN, linfct = c("regionthalamus:groupconcussion=0",
                                  "regionpallidostriatum:groupconcussion=0",
                                  "regionneocortex:groupconcussion=0",
                                  "regionmidbrain:groupconcussion=0",
                                  "regionpons:groupconcussion=0",
                                  "regioncingulateGyrus:groupconcussion=0",
                                  "regionhippocampus:groupconcussion=0",
                                  "regionsupramarginalGyrus:groupconcussion=0",
                                  "regioncorpusCallosum:groupconcussion=0"
                                  ))
rownames(e.glhtUN$linfct) <- c("thalamus", "pallidostriatum", "neocortex", "midbrain", "pons",
                               "cingulate.gyrus", "hippocampus", "supramarginal.gyrus", "corpus.callosum")    

summary(e.glhtUN)
## Linear Hypotheses:
##                          Estimate Std. Error z value Pr(>|z|)    
## thalamus == 0             0.13061    0.04514   2.894  0.01862 *  
## pallidostriatum == 0      0.13536    0.04175   3.242  0.00630 ** 
## neocortex == 0            0.05692    0.03985   1.428  0.44025    
## midbrain == 0             0.13831    0.03976   3.479  0.00281 ** 
## pons == 0                 0.03494    0.04473   0.781  0.89288    
## cingulate.gyrus == 0      0.20342    0.03500   5.813  < 0.001 ***
## hippocampus == 0          0.15788    0.03765   4.193  < 0.001 ***
## supramarginal.gyrus == 0  0.06198    0.03883   1.596  0.33946    
## corpus.callosum == 0      0.16897    0.03569   4.734  < 0.001 ***

## * Mixed model with unstructured covariance matrix
mSPECT0.lvm <- lvm(c(log.thalamus, 
                     log.pallidostriatum, 
                     log.neocortex, 
                     log.midbrain, 
                     log.pons, 
                     log.cingulateGyrus, 
                     log.hippocampus, 
                     log.supramarginalGyrus, 
                     log.corpusCallosum)~genotype+group+eta)
latent(mSPECT0.lvm) <- ~eta

eSPECT0.lvm <- estimate(mSPECT0.lvm, data = dfW.SPECT, control = list(constrain = TRUE))

mSPECT.lvm <- mSPECT0.lvm
covariance(mSPECT.lvm) <- log.supramarginalGyrus~log.neocortex

eSPECT.lvm <- estimate(mSPECT.lvm, data = dfW.SPECT, control = list(constrain = TRUE, start = coef(eSPECT0.lvm)))
modelsearch(eSPECT.lvm)




e.glhtLVM <- glht2(eSPECT.lvm, linfct = c("log.thalamus~groupconcussion=0",
                                          "log.pallidostriatum~groupconcussion=0",
                                          "log.neocortex~groupconcussion=0",
                                          "log.midbrain~groupconcussion=0",
                                          "log.pons~groupconcussion=0",
                                          "log.cingulateGyrus~groupconcussion=0",
                                          "log.hippocampus~groupconcussion=0",
                                          "log.supramarginalGyrus~groupconcussion=0",
                                          "log.corpusCallosum~groupconcussion=0"
                                          ))




##----------------------------------------------------------------------
### multivariate-regression.R ends here
