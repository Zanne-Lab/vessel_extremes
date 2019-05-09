## Preliminary data preparation
library(bayou)
library(treeplyr)
library(foreach)
library(doParallel)
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

source("./regressionModels_2019.R")
registerDoParallel(cores=length(priors))

foreach(i=2:5) %dopar% {
  mod <- names(priors)[i]
  mcmc <- bayou.makeMCMC(tree, dat, .pred, SE=0.01, model=models[[mod]]$model, samp=100, prior = priors[[mod]], 
                         new.dir = "../../output/lnVs_bayou/", outname=paste(mod, "_r2", sep=""), startpar=models[[mod]]$startpar, plot.freq=NULL)
  
  mcmc$run(300000)
  
  chain <- mcmc$load()
  chain <- set.burnin(chain, 0.3)
  
  if(length(grep("N", mod))==1){
    shiftsum <- shiftSummaries(chain, mcmc, pp.cutoff=0.3)
    pdf(paste("../../output/lnVs_bayou/shiftsum_", mod, "_r1.pdf", sep=""))
    plotSimmap.mcmc(chain, burnin=0.3, pal=viridis::viridis)
    plotShiftSummaries(shiftsum, pal=viridis::viridis)
    dev.off()
  }
  
  
  saveRDS(chain, paste("../../output/lnVs_bayou/chain_", mod, "_r1.rds", sep=""))
  saveRDS(mcmc, paste("../../output/lnVs_bayou/mcmc_", mod, "_r1.rds", sep=""))
}


registerDoParallel(cores=50)
Bk <- qbeta(seq(0,1, length.out=50), 0.3,1)

mlnL <- foreach(i=1:length(priors)) %do% {
  mod <- names(priors)[i]
  mcmc <-  readRDS(paste("../../output/lnVs_bayou/mcmc_", mod, "_r1.rds", sep=""))
  chain <- readRDS(paste("../../output/lnVs_bayou/chain_", mod, "_r1.rds", sep=""))
  ss <- mcmc$steppingstone(300000, chain, Bk, burnin=0.3, plot=FALSE)
  saveRDS(ss, file=paste("../../output/lnVs_bayou/ss_", mod, "_r2.rds", sep=""))
  ss$lnr
}

filenames <- list.files("../../output/lnVs_bayou/")
ssfn <- filenames[grep("ss_", filenames)]
mlnL <- list()
for(i in 1:length(ssfn)){
  mlnL[[ssfn[i]]] <- readRDS(paste("../../output/lnVs_bayou/", ssfn[i], sep=""))$lnr
}
