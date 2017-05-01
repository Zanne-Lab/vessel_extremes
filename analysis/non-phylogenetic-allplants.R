# Non-phylogenetic plots and models
# Will Pearse - 2017-04-30

# Load data
growth.form <- read.csv("../speciesTraitDataAEZ3.csv", as.is=TRUE)
env <- read.csv("../datasets/species_summaries.csv", as.is=TRUE)
data <- merge(growth.form, env, by.x="gs", by.y="X")
data <- data[!duplicated(data$gs),]
data$woodiness <- 

# Subset data to only the needed variables; rescale
data <- data[,c("gs","tmin.025","pmin.025","tseas.975","pseas.975","gbif.records","woodiness")]
orig.data <- data
data$tmin.025 <- as.numeric(scale(data$tmin.025))
data$pmin.025 <- as.numeric(scale(data$pmin.025))
data$tseas.975 <- as.numeric(scale(data$tseas.975))
data$pseas.975 <- as.numeric(scale(data$pseas.975))

# Plot data across bands
.prop.woody <- function(x, var, data){
    x <- data[,var] <= x
    if(sum(x)==0)
        return(NA)
    return(sum(data$woodiness[x])/sum(x))
}

range <- with(data, range(c(tmin.025,pmin.025,tseas.975,pseas.975)))
plot(0, xlim=range, ylim=c(0,1), type="n")
for(prop in seq(range[1],range[2],by=.1)){
    points(.prop.woody(prop, "tmin.025", data) ~ prop, col="black", cex=.5)
    points(.prop.woody(prop, "pmin.025", data) ~ prop, col="red", cex=.5)
    points(.prop.woody(prop, "tseas.975", data) ~ prop, col="blue", cex=.5)
    points(.prop.woody(prop, "pseas.975", data) ~ prop, col="green", cex=.5)
}

