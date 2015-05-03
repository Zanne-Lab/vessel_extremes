#First pass at conduit analysis
#Will Pearse - 2014-05-15
#NOTE: I don't really know what climate varable(s) I'm meant to be using! I've had a quick look and think, but these could be silly
#NOTE: Waaaaay too many NAs in growthForm, but I'm pushing through as this is all preliminary (and smaller --> more tractable on a plane!)
###############################
#HEADERS#######################
###############################
require(caper)
require(phytools)
require(OUwie)
require(geiger)
#Below can be found online at github.com/willpearse/willeerd/R/clean.taxonomy.R
source("~/Code/willeerd/R/clean.taxonomy.R")

###############################
#DATA##########################
###############################
#Load
tree <- read.tree("~/Dropbox/ConduitAnatomyGFClimate/Vascular_Plants_rooted.dated.tre")
climate <- read.csv("~/Dropbox/ConduitAnatomyGFClimate/WoodyClimateHemi.12May.csv")
traits <- read.csv("~/Dropbox/ConduitAnatomyGFClimate/CombinedGFConduit.csv")
#Match traits and bind missing onto phylogeny
traits$gs <- gsub(" ", "_", traits$gs)
missing.spp <- setdiff(traits$gs, tree$tip.label)
missing.gen <- gsub("_[a-z]*", "", missing.spp)
tree <- congeneric.merge(missing.spp, tree, split="_", cite=FALSE)
tree <- drop.tip(tree, setdiff(tree$tip.label, traits$gs))
tree$node.label <- NULL
#Match climate data onto phylogeny and traits
climate$species <- gsub(" ", "_", as.character(climate$species))
climate <- climate[climate$species %in% tree$tip.label,]
#Match all onto one-another and exclude unneeded data
traits <- traits[traits$gs %in% climate$species, 1:4]
traits <- traits[order(traits$gs), ]
climate <- climate[order(climate$species), c("species", "tcoldq.me", "pdryq.me", "alt.me")]
identical(climate$species, traits$gs)
data <- cbind(traits[,-1], climate)
c.data <- comparative.data(tree, data, names.col=species)
#Log response variables
c.data$data$log.size <- log10(c.data$data$vesselSize)
c.data$data$log.no <- log10(c.data$data$vesselNumber)
save.image("wip.RData")

###############################
#ANALYSIS######################
###############################

###############################
#Conduit ~ climate + growth####
###############################
basic <- lm(log.size ~ tcoldq.me + pdryq.me + alt.me + growthForm, data=c.data$data)
lam.size <- pgls(log.size ~ tcoldq.me + pdryq.me + alt.me + growthForm, data=c.data, lambda="ML")
del.size <- pgls(log.size ~ tcoldq.me + pdryq.me + alt.me + growthForm, data=c.data, delta="ML")
kap.size <- pgls(log.size ~ tcoldq.me + pdryq.me + alt.me + growthForm, data=c.data, kappa="ML")
lam.no <- pgls(log.no ~ tcoldq.me + pdryq.me + alt.me + growthForm, data=c.data, lambda="ML")
del.no <- pgls(log.no ~ tcoldq.me + pdryq.me + alt.me + growthForm, data=c.data, delta="ML")
kap.no <- pgls(log.no ~ tcoldq.me + pdryq.me + alt.me + growthForm, data=c.data, kappa="ML")
#...Will would be pleased; they all pick things up (different models?).
#I think I'm hitting the edge of parameter space with Delta/kappa so that might need some work.

###############################
#Model fitting and peaks#######
###############################
o.size <- data.frame(species=rownames(c.data$data), regimes=as.numeric(c.data$data$growthForm), size=c.data$data$log.size)
size.ou <- OUwie(c.data$phy, o.size, "BM1")
o.size
