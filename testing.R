# Load original data
data <- read.csv("WoodyClimateHemi.12May.csv")
data$herb <- data$woodiness=="H"
original <- data

# Load new data
growth.form <- read.csv("speciesTraitDataAEZ3.csv", as.is=TRUE)
env <- read.csv("datasets/species_summaries_all.csv", as.is=TRUE)
env$species <- read.csv("datasets/species_list.csv", as.is=TRUE)[,2]
data <- merge(growth.form, env, by.x="gs", by.y="species")
data <- data[!duplicated(data$gs),]


missing <- setdiff(data$gs, original$species)
common <- intersect(data$gs, original$species)

original <- setNames(original$southLat.lo, original$species)
new <- setNames(data$decimallatitude.025, data$gs)
new <- new[names(new) %in% names(original)]
original <- original[names(original) %in% names(new)]
new <- new[order(names(new))]
original <- original[order(names(original))]
identical(names(new), names(original))
#new <- new/10

plot(new, original)
cor.test(new, original)

plot.herb <- function(var, data, by, ...){
    bins <- with(data, seq(min(data[,var]),max(data[,var]),by=by))
    x <- with(data, table(cut(data[,var],bins), herb))[-1,]
    weight <- as.numeric(scale(rowSums(x)[-1]))
    weight <- (weight - min(weight) + 1) * 2
    x <- apply(x, 1, function(y) y[2]/sum(y))
    plot(x ~ bins[-1:-2], type="n", axes=FALSE, ylab="Proportion herbaceous", ylim=c(0,1), lwd=weight, ...)
    for(i in seq_along(x))
        lines(bins[(i+2):(i+3)], x[i:(i+1)], lwd=weight[i])
    axis(1)
    axis(2)
}

data <- data[,c("species","herb","tmin.lo")]
data <- na.omit(data)

pdf("~/Desktop/demo.pdf")
plot.herb("tmin.lo", data, by=5)
dev.off()
