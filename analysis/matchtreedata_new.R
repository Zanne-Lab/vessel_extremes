## # 
require(devtools)
require(bayou)
require(aRbor)
require(pez)
require(caper)
require(dplyr)
install_github("richfitz/phyndr")
install_github("wcornwell/TaxonLookup")
require(phyndr)
library(TaxonLookup)
setwd("~/repos/growthforms/analysis/")
data(plant_lookup)

tree <- read.tree("../Vascular_Plants_rooted.dated.tre")
climate <- read.csv("../WoodyClimateHemi.12May.csv")
#traits <- read.csv("../CombinedGFConduit.csv")
growthforms <- read.csv("../datasets/growthFormData.csv")
climate2 <- read.table("../datasets/spp_summaries.txt")
vesselanatomy <- read.csv("../datasets/T&MGFConduit.csv")
additional <- read.csv("../datasets/growthFormDataTropicosAdditions.csv")
vesselanatomy <- vesselanatomy[!duplicated(vesselanatomy$gs.tpl1.1),]
colnames(additional)[1] <- colnames(growthforms)[1] <- colnames(climate[1]) <- colnames(vesselanatomy)[1] <- "species"
growthforms <- rbind(growthforms, additional)
growthforms <- group_by(growthforms, species)
ugrowthforms <- summarize(growthforms, support=unique(support))
vac <- left_join(vesselanatomy, climate2)
vacf <- left_join(vac, growthforms)
vacf <- vacf[!duplicated(vacf$species),]
climate$species <- gsub(" ", "_", as.character(climate$species))
climate <- climate[climate$species %in% tree$tip.label,]
vacf$species <- gsub(" ", "_", as.character(vacf$species))
dat <- left_join(vacf,dplyr::select(climate, species, tcoldq.me, pdryq.me, alt.me))
rownames(dat) <- dat$species

tax <- lookup_table(unique(c(tree$tip.label, rownames(dat))), by_species = TRUE)
tree.missing <- setdiff(tree$tip.label, rownames(tax))
dat.missing <- setdiff(rownames(dat), rownames(tax))

phylook <- phyndr_taxonomy(ptree, data_species=pdat, taxonomy=tax)

td <- make.treedata(tree, dat)
apply(td$dat, 2, function(x) (length(x)-sum(is.na(x))))
tot =0
x=0
for(i in 1:100){
   phy <- phyndr::phyndr_sample(phylook)
   td <- make.treedata(phy, dat)
   suppN <- sum(td$dat$support=="C", na.rm=TRUE)
   suppNtot <- apply(td$dat, 2, function(x) (length(x)-sum(is.na(x))))['support']
   if(suppN == x){
     if(suppNtot>tot){
       TD <- td
       tot = suppNtot
       print(suppNtot)
     }
   }
   if(suppN > x){
     TD <- td
     x = suppN
     print(suppN)
   }

}
td <- TD
#
#

head(td$dat)
setdiff(td$phy$tip.label[!is.na(td$dat$pdryq.me)], td$phy$tip.label[!is.na(td$dat$precip)])

saveRDS(td, file="../output/cleandata/matchednewdata.rds")

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
#c.data$data$log.size <- log10(c.data$data$vesselSize)
#c.data$data$log.no <- log10(c.data$data$vesselNumber)
#td <- make.treedata(c.data$phy, c.data$data)





o