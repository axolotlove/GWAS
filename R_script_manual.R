library(PhenotypeSimulator)

genotypefile  <- system.file("/home/axolotlove/GWAS/EUR_1KG/",
                               "EUR_1KG_genome_nodup.bed",
                               package = "PhenotypeSimulator")
genotypefile <- gsub("\\.bed", "", genotypefile)

genVar <- 0.4
noiseVar <- 1 - genVar
totalSNPeffect <- 0.1
h2s <- totalSNPeffect/genVar
phi <- 0.6 
rho <- 0.2
delta <- 0.2
shared <- 1
independent <- 1 - shared


total <- shared * h2s *genVar +  independent * h2s * genVar +
  shared * (1-h2s) * genVar +   independent * (1-h2s) * genVar +
  shared * phi* noiseVar +  independent * phi* noiseVar +
  rho * noiseVar +
  shared * delta * noiseVar +  independent * delta * noiseVar

total == 1

phenotype2 <- runSimulation(N = 100, P = 3, genotypefile="/home/axolotlove/GWAS/EUR_1KG/EUR_1KG_genome_nodup", 
                            format ="plink",
                            cNrSNP=10, 
                            mBetaGenetic = 0, sdBetaGenetic = 0.2,
                            theta=1,
                            genVar = genVar, h2s = h2s,
                            phi = 0.6, delta = 0.2, rho=0.2,
                            NrFixedEffects = 2, NrConfounders = c(2, 2),
                            distConfounders = c("bin", "norm"),
                            probConfounders = 0.2,
                            verbose = TRUE )

varGenFixed <- t(phenotype2$varComponents
                 [grepl("var_genFix", names(phenotype2$varComponents))])
varGenBg <- t(phenotype2$varComponents
              [grepl("var_genBg", names(phenotype2$varComponents))])

varNoiseFixed <- t(phenotype2$varComponents
                   [grepl("var_noiseFixed", names(phenotype2$varComponents))])
varNoiseBg <- t(phenotype2$varComponents
                [grepl("var_noiseBg", names(phenotype2$varComponents))])
varNoiseCorr <- t(phenotype2$varComponents
                  [grepl("var_noiseCor", names(phenotype2$varComponents))])

totalPropVariance <- as.matrix(t(data.frame(varGenFixed, 
                                            varGenBg, 
                                            varNoiseFixed,
                                            varNoiseBg, 
                                            varNoiseCorr=c(varNoiseCorr, 0))))
totalPropVariance <- cbind(totalPropVariance, rowSums(totalPropVariance))
totalPropVariance <- rbind(totalPropVariance, 
                           sumProportion=colSums(totalPropVariance))

colnames(totalPropVariance) <- c("shared effect", "independent effect", 
                                 "total effect")

knitr::kable(totalPropVariance, caption="Proportion of variance explained
             by the different phenotype components")

par(mar = rep(3, 6))
image(t(phenotype$phenoComponentsFinal$Y), main="Phenotype: [samples x traits]", 
      axes=FALSE, cex.main=0.8)
mtext(side = 1, text = "Samples", line = 1)
mtext(side = 2, text = "Traits", line = 1)
image(cor(phenotype$phenoComponentsFinal$Y), 
      main="Correlation of traits [traits x traits]", axes=FALSE, cex.main=0.8)
mtext(side = 1, text = "Traits", line = 1)
mtext(side = 2, text = "Traits", line = 1)

out <- savePheno(phenotype2, directory="/home/axolotlove/GWAS/PhenotypeSimulator",
                 outstring="pheno_simulation",
                 format=c("csv","plink"),
                 verbose=TRUE)

