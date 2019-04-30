## Preliminary data preparation
library(bayou)
library(foreach)
library(doParallel)
registerDoParallel(cores=length(priors))
source("./regressionModels_2019.R")
td <- readRDS("../../output/newtd_growthform.rds")

## Analysis of vessel size vs. vessel number. 
tree <- td$phy
dat <- td$dat
tree <- multi2di(tree, random=FALSE) # Resolve polytomies
tree$edge.length[tree$edge.length==0] <- .Machine$double.eps # Bayou doesn't like 0 length branches

pred <- dat
.pred <- dplyr::select(pred, logN, growthC, tempK, precip)
.pred[,c(1,3,4)] <- as.data.frame(apply(.pred[,c(1,3,4)], 2, scale))
dat <- scale(td[['logS']])
dat <- dat[,1]

mod <- names(priors)[15]
mod
foreach(i=1:15) %dopar% {
  mcmc <- bayou.makeMCMC(tree, dat, .pred, SE=0.01, model=models[[mod]]$model, samp=100, prior = priors[[mod]], 
                         new.dir = "../../output/lnVs_bayou/", outname=paste(mod, "_r1", sep=""), startpar=models[[mod]]$startpar, plot.freq=NULL)
  
  mcmc$run(10000)
  
  chain <- mcmc$load()
  
  shiftsum <- shiftSummaries(chain, mcmc, pp.cutoff=0.3)
  pdf(paste("../../output/lnVs_bayou/shiftsum_", mod, "_r1.pdf", sep=""))
  plotShiftSummaries(shiftsum, pal=viridis::viridis)
  dev.off()
  
  saveRDS(chain, paste("../../output/lnVs_bayou/chain_", mod, "_r1.pdf", sep=""))
  saveRDS(mcmc, paste("../../output/lnVs_bayou/mcmc_", mod, "_r1.pdf", sep=""))
}
