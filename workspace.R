#Headers
require(pez)
require(caper)

#Load traits and tree, match
tree <- read.tree("data/Vascular_Plants_rooted.dated.tre")
data <- read.csv("data/speciesTraitDataAEZ3.csv")
data <- data[,c("gs", "phenology", "woodiness", "vesselSize")]
data <- na.omit(data)
data$gs <- gsub(" ", "_", data$gs)
tree <- congeneric.merge(tree, data$gs)

#Load environmental data
env <- read.table("data/spp_summaries.txt")
env$species <- gsub(" ", "_", env$species)
data <- merge(data, env, by.x="gs", by.y="species")
tree$node.label <- NULL

#Load geographical data
lats <- read.table("data/species_lats.txt")
lats$species <- gsub(" ", "_", lats$species)
data <- merge(data, lats, by.x="gs", by.y="species")
data$pole.lim <- apply(data, 1, function(x) max(abs(as.numeric(x[c("lat.25","lat.975")]))))

#Create comparative.data
c.data <- comparative.data(tree, data, gs)
