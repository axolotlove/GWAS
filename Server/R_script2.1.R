library(PhenotypeSimulator)
setwd('~/GWAS/Server/pathway/')


genotypes <- readStandardGenotypes(N=1000, filename='genotypes_hapgen.controls.caspase',
                                   format="plink")


causalSNPs <-getCausalSNPs(N=1000, genotypes = genotypes$genotypes,NrCausalSNPs = 15, verbose = FALSE)
genFixed <-geneticFixedEffects(N = 1000, P = 10, X_causal = causalSNPs)

genotypes <- readStandardGenotypes(N=1000, filename='genotypes_subset_hapgen.controls',
                                   format="plink")

genotypes_sd <-standardiseGenotypes(genotypes$genotypes)
kinship <- getKinship(N=1000, X=genotypes_sd, verbose = FALSE)

genBg <-geneticBgEffects(N=1000, kinship = kinship, P = 10)

noiseFixed <-noiseFixedEffects(N = 1000, P = 10, NrFixedEffects = 4,
                               NrConfounders =c(1, 2, 1, 2),pIndependentConfounders =c(0, 1, 1, 0.5),
                               distConfounders =c("bin", "cat_norm","cat_unif", "norm"),
                               probConfounders = 0.2,catConfounders =c(3, 4))

correlatedBg <-correlatedBgEffects(N = 1000, P = 10, pcorr = 0.8)
noiseBg <-noiseBgEffects(N = 1000, P = 10)

genVar <- 0.3
noiseVar <- 1 - genVar
totalSNPeffect
h2s <- 0.8
phi <- 0.6 
rho <- 0.1
delta <- 0.3
shared <- 0.95
independent <- 1 - shared

genFixed_shared_scaled <- rescaleVariance(genFixed$shared, shared * h2s *genVar)
genFixed_independent_scaled <- rescaleVariance(genFixed$independent, 
                                               independent * h2s *genVar)
genBg_shared_scaled <- rescaleVariance(genBg$shared, shared * (1-h2s) *genVar)
genBg_independent_scaled <- rescaleVariance(genBg$independent, 
                                            independent * (1-h2s) * genVar)

noiseBg_shared_scaled <- rescaleVariance(noiseBg$shared, shared * phi* noiseVar)
noiseBg_independent_scaled <- rescaleVariance(noiseBg$independent, 
                                              independent * phi* noiseVar)
correlatedBg_scaled <- rescaleVariance(correlatedBg$correlatedBg, 
                                       shared * rho * noiseVar)

noiseFixed_shared_scaled <- rescaleVariance(noiseFixed$shared, shared * delta * 
                                              noiseVar)
noiseFixed_independent_scaled <- rescaleVariance(noiseFixed$independent, 
                                                 independent * delta * noiseVar)

total <- shared * h2s *genVar +  independent * h2s * genVar +
  shared * (1-h2s) * genVar +   independent * (1-h2s) * genVar +
  shared * phi* noiseVar +  independent * phi* noiseVar +
  rho * noiseVar +
  shared * delta * noiseVar +  independent * delta * noiseVar

total == 1

Y <- scale(genBg_shared_scaled$component + noiseBg_shared_scaled$component +
             genBg_independent_scaled$component + noiseBg_independent_scaled$component +
             genFixed_shared_scaled$component + noiseFixed_shared_scaled$component +
             genFixed_independent_scaled$component + noiseFixed_independent_scaled$component +
             correlatedBg_scaled$component)        


row.names(Y) <- genotypes[["id_samples"]]


image(t(Y), main="Phenotype: [samples x traits]", 
      axes=FALSE, cex.main=0.8)


write.table(Y, 'Y_caspase_second', sep = '\t', col.names = FALSE)
