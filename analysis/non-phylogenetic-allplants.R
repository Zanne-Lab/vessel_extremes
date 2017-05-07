# Non-phylogenetic plots and models
# Will Pearse - 2017-04-30

library(MuMIn)

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

pdf("../output/figures/figure_1_temp.pdf")
plot.herb("tmin.025", data, by=5, xlab=expression(paste("Minimum temperature (",degree,")")))
dev.off()
pdf("../output/figures/figure_1_precip.pdf")
plot.herb("pmin.025", data, by=20, xlim=c(0,250), xlab="Minimum precipitation (mm)")
dev.off()
pdf("../output/figures/figure_1_temp_seas.pdf")
plot.herb("tseas.975", data, by=500, xlim=c(0,15000), xlab="Temperature seasonality")
dev.off()
pdf("../output/figures/figure_1_precip_seas.pdf")
plot.herb("pseas.975", data, by=20, xlim=c(0,200), xlab="Temperature seasonality")
dev.off()
pdf("../output/figures/figure_1_lat_lower.pdf")
plot.herb("decimallatitude.025", data, by=5, xlim=c(-60,0), xlab=expression(paste("Minimum Southern latitude (",degree,")")))
dev.off()
pdf("../output/figures/figure_1_lat_upper.pdf")
plot.herb("decimallatitude.975", data, by=5, xlim=c(0,90), xlab=expression(paste("Maximum Northern latitude (",degree,")")))
dev.off()

# Let's do some modelling (?)

model <- glm(herb ~ tmin.025 + I(log10(pmin.025+1)) + pseas.975 + tseas.975, data=data, na.action="na.pass")
model.set <- dredge(model)
model.avg <- model.avg(model.set, subset=delta<4)
summary(model.avg)

with(data, boxplot(log10(pmin.025+1) ~ herb))
with(data, boxplot(tmin.025 ~ herb))
