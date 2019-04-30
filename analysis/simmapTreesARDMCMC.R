library(devtools)
library(bayou)
library(treeplyr)
#library(pez)
library(caper)
library(dplyr)
#install_github("richfitz/phyndr")
#install_github("wcornwell/taxonlookup")
library(phyndr)
library(taxonlookup)
library(diversitree)
library(OUwie)
library(foreach)
library(doParallel)
library(iterators)
library(reshape2)

setwd("~/repos/vessel_extremes/analysis/")

## Read in trees, climate data, taxonomy and growth form data. 
tree <- read.tree("../datasets/zanne1.1.tre")
AEZ3 <- read.csv("../speciesTraitDataAEZ3.csv")
#climate <- read.csv("../WoodyClimateHemi.12May.csv")
species.summaries <- read.csv("../datasets/species_summaries_all.csv")
plant_lookup <- plant_lookup()

## Prepare names for matching
AEZ3$species <- gsub(" ", "_", AEZ3$gs)
#climate$species <- gsub(" ", "_", climate$species)
#dat.old <- left_join(climate, AEZ3)
#colnames(species.summaries)[1] <- "species"
#rownames(dat.old) <- dat.old$species
colnames(species.summaries)[1] <- "species"
species.summaries$species <- gsub(" ", "_", species.summaries$species)
dat.new <- left_join(species.summaries, AEZ3)
rownames(dat.new) <- gsub(" ", "_", dat.new$species)

tax <- lookup_table(unique(c(tree$tip.label, rownames(dat.new))), by_species = TRUE)
tree.missing <- setdiff(tree$tip.label, rownames(tax))
dat.missing <- setdiff(rownames(dat.new), rownames(tax))
ptree <- drop.tip(tree, tree.missing)
#pdat <- rownames(dat.new)[-match(dat.missing, dat.new$species)]
nH <- nodeHeights(ptree)
ptree$edge.length[ptree$edge[,2] <= length(ptree$tip.label)] <- 
  ptree$edge.length[ptree$edge[,2] <= length(ptree$tip.label)] + (max(branching.times(ptree)) - nH[ptree$edge[,2] <= length(ptree$tip.label), 2])
is.ultrametric(ptree)

tdOUwie <- make.treedata(ptree, dat.new) 
tdOUwie
tdOUwie <- make.treedata(ptree, dat.new) 
tdOUwie <- filter(tdOUwie, !is.na(vesselSize), !is.na(phenology),!is.na(tmin.01))
tdOUwie <- dplyr::select(tdOUwie, vesselSize, phenology, starts_with("tmin"))
tdOUwie <- mutate(tdOUwie, OU2freeze.01 = c("NF", "F")[as.numeric(tmin.01 < 0)+1],
                  OU2freeze.025 = c("NF", "F")[as.numeric(tmin.025 < 0)+1],
                  OU2freeze.05 = c("NF", "F")[as.numeric(tmin.05 < 0)+1],
                  OU2freeze.50 = c("NF", "F")[as.numeric(tmin.50 < 0)+1],
                  OU2freeze.95 = c("NF", "F")[as.numeric(tmin.95 < 0)+1], 
                  OU2freeze.975 = c("NF", "F")[as.numeric(tmin.975 < 0)+1], 
                  OU2freeze.99 = c("NF", "F")[as.numeric(tmin.99 < 0)+1]
)
tdOUwie <- filter(tdOUwie, phenology != "D_EV")
#,
#OU2freeze.mean = c("NF", "F")[as.numeric(tmin.mean < 0)+1], 
#OU2freeze.lower = c("NF", "F")[as.numeric(tmin.lower < 0)+1],
#OU2freeze.upper = c("NF", "F")[as.numeric(tmin.upper < 0)+1])

pasteNA <- function(x, y, sep=" ", collapse=NULL){
  ifelse(is.na(y), NA, paste(x,y, sep=sep, collapse=collapse))
}

tdOUwie <- mutate(tdOUwie, OU1="global", OU2phenology=phenology, lnVs=log(vesselSize))
tdOUwie <- mutate(tdOUwie, OU4freeze.01= paste(phenology, OU2freeze.01, sep="_"), 
                  OU4freeze.025= paste(phenology, OU2freeze.025, sep="_"),
                  OU4freeze.05= paste(phenology, OU2freeze.05, sep="_"),
                  OU4freeze.50= paste(phenology, OU2freeze.50, sep="_"),
                  OU4freeze.95= paste(phenology, OU2freeze.95, sep="_"),
                  OU4freeze.975= paste(phenology, OU2freeze.975, sep="_"),
                  OU4freeze.99= paste(phenology, OU2freeze.99, sep="_")
)#,
#OU4freeze.mean= paste(phenology, OU2freeze.mean, sep="_"), 
#OU4freeze.lower = pasteNA(phenology, OU2freeze.lower, sep="_"),
#OU4freeze.upper = pasteNA(phenology, OU2freeze.upper, sep="_")

tdOUwie <- treeply(tdOUwie, reorder, "postorder")
tdOUwieFreeze <- dplyr::select(tdOUwie, lnVs, OU1, starts_with("OU"))
par(mfrow=c(5,4), mar=c(2, 2, 0.1, 0.1))
dum <- lapply(2:16, function(x) boxplot(tdOUwieFreeze[['lnVs']] ~ tdOUwieFreeze[[x]]))

#

tdOUwieFreeze$phy$edge.length <- tdOUwieFreeze$phy$edge.length/max(branching.times(tdOUwieFreeze$phy))

require(foreach)
require(doParallel)
registerDoParallel(cores=80)
itd <- iter(tdOUwieFreeze$dat[,3:ncol(tdOUwieFreeze$dat)], by="column")

smtds <- foreach(i=itd) %do% {
  idat <- data.frame(tdOUwieFreeze$phy$tip.label, "lnVs"=tdOUwieFreeze$dat$lnVs,i)
  itdx <- make.treedata(tdOUwieFreeze$phy, idat) 
  itdx <- filter(itdx, !is.na(itdx$dat[[2]]))
}

aic <- list()
for(i in 1:length(smtds)){
  ard.fit <- fitDiscrete(smtds[[i]]$phy, smtds[[i]][[2]], model="ARD")
  er.fit <- fitDiscrete(smtds[[i]]$phy, smtds[[i]][[2]], model="ER")
  aic[[i]] <- cbind(ard.fit$opt$aicc, er.fit$opt$aicc)
  #ace(smtds[[i]][[2]], smtds[[i]]$phy, type="discrete", method="ML", model="ARD", marginal=TRUE)
}
aic <- do.call(rbind, aic)
rownames(aic) <- colnames(tdOUwieFreeze$dat)[3:17]

models <- ifelse(aic[,1] < aic[,2] - 4, "ARD", "ER")

itd <- iter(tdOUwieFreeze$dat[,6:ncol(tdOUwieFreeze$dat)], by="column")

smtrees.ARDmcmc <- foreach(i=itd) %dopar% {
  idat <- cbind(tdOUwieFreeze$phy$tip.label, i)
  itdx <- make.treedata(tdOUwieFreeze$phy, idat) 
  itdx <- filter(itdx, !is.na(itdx$dat[[1]]))
  make.simmap(itdx$phy, itdx[[1]], model="ARD", nsim=20, Q="mcmc", prior=list(beta=0.1, use.empirical=TRUE))
}


itd <- iter(tdOUwieFreeze$dat[,3:5], by="column")
smtrees.ERmcmc <- foreach(i=itd) %dopar% {
  idat <- cbind(tdOUwieFreeze$phy$tip.label, i)
  itdx <- make.treedata(tdOUwieFreeze$phy, idat) 
  itdx <- filter(itdx, !is.na(itdx$dat[[1]]))
  make.simmap(itdx$phy, itdx[[1]], model="ER", nsim=20, Q="mcmc", prior=list(beta=0.1, use.empirical=TRUE))
}

saveRDS(smtrees.ARDmcmc, file="../output/OUwie/smtrees.ARDmcmc2.rds")
saveRDS(smtrees.ERmcmc, file="../output/OUwie/smtrees.ERmcmc2.rds")
saveRDS(aic, file="../output/OUwie/ARDERquantileAICs2.rds")
#smtrees.ER <- readRDS("../output/smtrees_ER.rds")
#plotSimmap(smtrees.ARDmcmc[[9]][[6]])
#plotSimmap(smtrees.ER[[9]][[6]])

