## # 
require(bayou)
require(aRbor)
require(pez)
require(caper)
setwd("~/Dropbox/ConduitAnatomyGFClimate/")

tree <- read.tree("./Vascular_Plants_rooted.dated.tre")
climate <- read.csv("./WoodyClimateHemi.12May.csv")
traits <- read.csv("./CombinedGFConduit.csv")
#Match traits and bind missing onto phylogeny
traits$gs <- gsub(" ", "_", traits$gs)
missing.spp <- setdiff(traits$gs, tree$tip.label)
missing.gen <- gsub("_[a-z]*", "", missing.spp)
tree <- congeneric.merge(tree, missing.spp, split="_")
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
td <- make.treedata(c.data$phy, c.data$data)
saveRDS(c.data, file="./datasets/matcheddata.rds")




