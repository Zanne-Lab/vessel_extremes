# Non-phylogenetic plots and models
# Will Pearse - 2017-04-30

# Load data
growth.form <- read.csv("../speciesTraitDataAEZ3.csv", as.is=TRUE)
env <- read.csv("../datasets/species_summaries_all.csv", as.is=TRUE)
env$species <- read.csv("../datasets/species_list.csv", as.is=TRUE)[,2]
data <- merge(growth.form, env, by.x="gs", by.y="species")
data <- data[!duplicated(data$gs),]

# Subset data to only the needed variables; recode
data <- data[,c("gs","tmin.025","pmin.025","tseas.975","pseas.975","decimallatitude.025","decimallatitude.975","gbif.records","woodiness")]
data$tmin.025 <- data$tmin.025/10
data$herb <- data$woodiness!="W"
data <- na.omit(data)

# Do each separate plot
.prop.herb <- function(x, var, data){
    x <- data[,var] <= x
    if(sum(x)==0)
        return(NA)
    return(1-sum(data$herb[x])/nrow(data))
}
.plot.herb <- function(var, data, ...){
    range <- range(data[,var])
    plot(0, xlim=range, type="n", ylim=c(0,1),axes=FALSE, ...)
    props <- seq(range[1], range[2], length.out=200)
    for(i in seq(2, length(props)))
        lines(props[(i-1):i], c(.prop.herb(props[i-1],var,data),.prop.herb(props[i],var,data)), col="black")
    axis(1)
    axis(2)
}
pdf("../output/figures/figure_1_temp.pdf")
.plot.herb("tmin.025", data, xlab=expression(paste("Minimum temperature (",degree,")")), ylab="Proportion herbaceous")
dev.off()
pdf("../output/figures/figure_1_precip.pdf")
.plot.herb("pmin.025", data, xlab="Minimum precipitation (mm)", ylab="Proportion herbaceous")
dev.off()
pdf("../output/figures/figure_1_temp_seas.pdf")
.plot.herb("tseas.975", data, xlab="Temperature seasonality", ylab="Proportion herbaceous")
dev.off()
pdf("../output/figures/figure_1_precip_seas.pdf")
.plot.herb("pseas.975", data, xlab="Temperature seasonality", ylab="Proportion herbaceous")
dev.off()
pdf("../output/figures/figure_1_lat_lower.pdf")
.plot.herb("decimallatitude.025", data, xlab=expression(paste("Minimum Southern latitude (",degree,")")), ylab="Proportion herbaceous")
dev.off()
pdf("../output/figures/figure_1_lat_upper.pdf")
.plot.herb("decimallatitude.975", data, xlab=expression(paste("Maximum Northern latitude (",degree,")")), ylab="Proportion herbaceous")
dev.off()
