### anonymize.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar 28 2022 (09:25) 
## Version: 
## Last-Updated: mar 31 2022 (17:10) 
##           By: Brice Ozenne
##     Update #: 24
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

library(lava)

## * SPECT

## ** define LVM model
mSPECT.lvm <- lvm(c(log.thalamus, 
                    log.pallidostriatum, 
                    log.neocortex, 
                    log.midbrain, 
                    log.pons, 
                    log.cingulateGyrus, 
                    log.hippocampus, 
                    log.supramarginalGyrus, 
                    log.corpusCallosum)~genotype+group+eta)
latent(mSPECT.lvm) <- ~eta
covariance(mSPECT.lvm) <- log.supramarginalGyrus~log.neocortex  

## ** extract model parameters
if(dir.exists("source")){
    ## import data
    df.SPECT <- readRDS("source/data-SPECT.rds")
    df.SPECT$group <- as.numeric(df.SPECT$group=="concussion")
    df.SPECT$genotype <- as.numeric(df.SPECT$genotype=="HAB")

    ## fit
    eSPECT.lvm <- estimate(mSPECT.lvm, data = df.SPECT, control = list(constrain=FALSE))
    logLik(eSPECT.lvm)
    ## 'log Lik.' 257.3579 (df=46)

    ## get coefficients
    eSPECT.coef <- c(summary(eSPECT.lvm)$coef[,"Estimate"],
                     "group" = mean(df.SPECT$group),
                     "gender" = mean(df.SPECT$gender=="Male"),
                     genotype = mean(df.SPECT$genotype))
}

## ** model parameters
eSPECT.coef <- c("log.thalamus~genotype" = 0.60113337,
                 "log.thalamus~group" = 0.11998744,
                 "log.thalamus~eta" = 1,
                 "log.pallidostriatum~genotype" = 0.57196912,
                 "log.pallidostriatum~group" = 0.11358296,
                 "log.pallidostriatum~eta" = 0.83188873,
                 "log.neocortex~genotype" = 0.57947786,
                 "log.neocortex~group" = 0.0428338,
                 "log.neocortex~eta" = 0.83046515,
                 "log.midbrain~genotype" = 0.61590585,
                 "log.midbrain~group" = 0.09895344,
                 "log.midbrain~eta" = 0.85229559,
                 "log.pons~genotype" = 0.53958244,
                 "log.pons~group" = 0.01544713,
                 "log.pons~eta" = 0.87156712,
                 "log.cingulateGyrus~genotype" = 0.65551138,
                 "log.cingulateGyrus~group" = 0.15935961,
                 "log.cingulateGyrus~eta" = 0.76966478,
                 "log.hippocampus~genotype" = 0.57525276,
                 "log.hippocampus~group" = 0.11900547,
                 "log.hippocampus~eta" = 0.77094843,
                 "log.supramarginalGyrus~genotype" = 0.57436083,
                 "log.supramarginalGyrus~group" = 0.05088963,
                 "log.supramarginalGyrus~eta" = 0.86072425,
                 "log.corpusCallosum~genotype" = 0.57192471,
                 "log.corpusCallosum~group" = 0.17415555,
                 "log.corpusCallosum~eta" = 0.6775107,
                 "log.thalamus~~log.thalamus" = 0.0130141,
                 "log.pallidostriatum~~log.pallidostriatum" = 0.01010281,
                 "log.neocortex~~log.neocortex" = 0.00806946,
                 "log.neocortex~~log.supramarginalGyrus" = 0.00602613,
                 "log.midbrain~~log.midbrain" = 0.00344643,
                 "log.pons~~log.pons" = 0.00977953,
                 "log.cingulateGyrus~~log.cingulateGyrus" = 0.00420049,
                 "log.hippocampus~~log.hippocampus" = 0.01151499,
                 "log.supramarginalGyrus~~log.supramarginalGyrus" = 0.00819216,
                 "log.corpusCallosum~~log.corpusCallosum" = 0.0160105,
                 "eta~~eta" = 0.06325209,
                 "log.thalamus" = 0,
                 "log.pallidostriatum" = 0.07631049,
                 "log.neocortex" = -0.05896311,
                 "log.midbrain" = 0.24974019,
                 "log.pons" = 0.28074318,
                 "log.cingulateGyrus" = 0.09516947,
                 "log.hippocampus" = 0.08466852,
                 "log.supramarginalGyrus" = -0.09783103,
                 "log.corpusCallosum" = -0.00509232,
                 "eta" = 1.43044384,
                 "group" = 0.38888889,
                 "gender" = 0.44444444,
                 "genotype" = 0.52777778)


## ** simulation model (with covariate distribution)
mSPECTSim.lvm <- mSPECT.lvm
distribution(mSPECTSim.lvm, ~group) <- binomial.lvm(size = 1, p = eSPECT.coef["group"])
distribution(mSPECTSim.lvm, ~genotype) <- binomial.lvm(size = 1, p = eSPECT.coef["genotype"])
distribution(mSPECTSim.lvm, ~gender) <- binomial.lvm(size = 1, p = eSPECT.coef["gender"])
## mSPECTSim.lvm  <- categorical(mSPECTSim.lvm, ~genotype, p = mean(df.SPECT$genotype=="HAB"), labels = c("MAB","HAB"))
eSPECTSim.lvm <- do.call(sim, args = c(list(x = mSPECTSim.lvm), as.list(eSPECT.coef[coef(mSPECTSim.lvm)])))
## estimate(mSPECT.lvm, sim(eSPECTSim.lvm,1e3))

## ** generated data
set.seed(10)
dfsim.SPECT <- sim(eSPECTSim.lvm,100,latent=FALSE)
dfsim.SPECT$group <- factor(dfsim.SPECT$group, levels = 0:1, labels = c("healthy", "concussion"))
dfsim.SPECT$genotype <- factor(dfsim.SPECT$genotype, levels = 0:1, labels = c("MAB", "HAB"))
dfsim.SPECT$gender <- factor(dfsim.SPECT$gender, levels = 0:1, labels = c("female", "male"))
dfsim.SPECT$id <- paste0("ID",1:NROW(dfsim.SPECT))
dfsim.SPECT <- dfsim.SPECT[,c("id","gender","group","genotype", setdiff(names(dfsim.SPECT),c("id","gender","group","genotype")))]
write.table(dfsim.SPECT,"data/data-SPECT.txt", row.names = FALSE, sep = ";")
## read.table("data/data-SPECT.txt", header = TRUE, sep = ";")

## * PET

## ** define LVM model
mPET.lvm <- lvm(neocortex.log ~ eta + age+sex+sert+mr+sb.per.kg,
                hippocampus.log ~ eta + age+sex+sert+mr+sb.per.kg,
                caudate.log ~ eta + age+sex+sert+mr+sb.per.kg,
                putamen.log ~ eta + age+sex+sert+mr+sb.per.kg,
                eta  ~ group
                )
covariance(mPET.lvm) <- putamen.log~caudate.log
latent(mPET.lvm) <- ~eta

## ** extract model parameters
if(dir.exists("source")){

    ## import data
    df.PET <- read.table("source/data-PET.txt", dec = ".", header = TRUE)

    ## fit
    ePET.lvm <- estimate(mPET.lvm, data = df.PET, control = list(constrain=FALSE))
    logLik(ePET.lvm)
    ## 'log Lik.' 613.151 (df=34)

    ## get coefficients
    ePET.coef <- c(summary(ePET.lvm)$coef[,"Estimate"],
                   "age" = mean(df.PET$age),
                   "age~~age" = var(df.PET$age),
                   "sert" = mean(df.PET$sert=="LALA"),
                   "mr" = mean(df.PET$mr=="trio"),
                   "sb.per.kg" = mean(df.PET$sb.per.kg),
                   "sb.per.kg~~sb.per.kg" = mean(df.PET$sb.per.kg))

}

## ** model parameters
ePET.coef <- c("neocortex.log~eta" = 1,
               "neocortex.log~age" = -0.00202731,
               "neocortex.log~sb.per.kg" = -2.45924156,
               "neocortex.log~sexMale" = 0.07623864,
               "neocortex.log~sertnonLALA" = 0.02495424,
               "neocortex.log~mrtrio" = -0.03974207,
               "eta~groupHealthy Control" = 0.07796427,
               "hippocampus.log~eta" = 1.09846887,
               "hippocampus.log~age" = -0.00094373,
               "hippocampus.log~sb.per.kg" = -2.15394497,
               "hippocampus.log~sexMale" = 0.02981759,
               "hippocampus.log~sertnonLALA" = 0.0003854,
               "hippocampus.log~mrtrio" = -0.08665295,
               "caudate.log~eta" = 0.92934557,
               "caudate.log~age" = -0.00348555,
               "caudate.log~sb.per.kg" = -0.88957266,
               "caudate.log~sexMale" = -0.02011877,
               "caudate.log~sertnonLALA" = 0.01968885,
               "caudate.log~mrtrio" = -0.04102281,
               "putamen.log~eta" = 0.89647695,
               "putamen.log~age" = -0.00417257,
               "putamen.log~sb.per.kg" = -1.76418045,
               "putamen.log~sexMale" = 0.03510004,
               "putamen.log~sertnonLALA" = -0.01218037,
               "putamen.log~mrtrio" = 0.00109894,
               "neocortex.log~~neocortex.log" = 0.005535,
               "eta~~eta" = 0.0150819,
               "hippocampus.log~~hippocampus.log" = 0.00830567,
               "caudate.log~~caudate.log" = 0.0086576,
               "caudate.log~~putamen.log" = 0.00275813,
               "putamen.log~~putamen.log" = 0.00486161,
               "neocortex.log" = 0,
               "eta" = -0.427904,
               "hippocampus.log" = 0.55344104,
               "caudate.log" = 1.7410477,
               "putamen.log" = 1.70379965,
               "age" = 27.09720122,
               "age~~age" = 65.23446655,
               "sert" = 0.29120879,
               "mr" = 0.20879121,
               "sb.per.kg" = 0.01513636,
               "sb.per.kg~~sb.per.kg" = 0.01513636)

## ** simulation model (with covariate distribution)
mPETSim.lvm <- mPET.lvm
distribution(mPETSim.lvm, ~sert) <- binomial.lvm(size = 1, p = ePET.coef["sert"])
distribution(mPETSim.lvm, ~mr) <- binomial.lvm(size = 1, p = ePET.coef["mr"])
distribution(mPETSim.lvm, ~age) <- gaussian.lvm(size = 1, mean = ePET.coef["age"], sd = sqrt(ePET.coef["age~~age"]))
distribution(mPETSim.lvm, ~sb.per.kg) <- gaussian.lvm(size = 1, mean = ePET.coef["sb.per.kg"], sd = sqrt(ePET.coef["sb.per.kg~~sb.per.kg"]))

ePETSim.lvm <- do.call(sim, args = c(list(x = mPETSim.lvm), as.list(ePET.coef[coef(mPETSim.lvm)])))
## estimate(mPET.lvm, sim(ePETSim.lvm,1e3))

## ** generated data
set.seed(10)
dfsim.PET <- sim(ePETSim.lvm,100,latent=FALSE)
dfsim.PET$sert <- factor(dfsim.PET$sert, levels = 0:1, labels = c("non-LALA", "LALA"))
dfsim.PET$mr <- factor(dfsim.PET$mr, levels = 0:1, labels = c("prisma", "trio"))
dfsim.PET$id <- paste0("ID",1:NROW(dfsim.PET))
dfsim.PET <- dfsim.PET[,c("id","age","sb.per.kg","sert","mr", setdiff(names(dfsim.PET),c("id","age","sb.per.kg","sert","mr")))]

write.table(dfsim.PET,"data/data-PET.txt", row.names = FALSE, sep = ";")
## read.table("data/data-PET.txt", header = TRUE, sep = ";")
 
##----------------------------------------------------------------------
### anonymize.R ends here
