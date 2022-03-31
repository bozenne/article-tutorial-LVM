### anonymize.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar 28 2022 (09:25) 
## Version: 
## Last-Updated: mar 28 2022 (18:59) 
##           By: Brice Ozenne
##     Update #: 13
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
write.table(dfsim.SPECT,"data/data-SPECT.txt")

##----------------------------------------------------------------------
### anonymize.R ends here
