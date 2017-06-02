# Non-phylogenetic plots and models
# Will Pearse - 2017-04-30

library(MuMIn)

# Load data
growth.form <- read.csv("../speciesTraitDataAEZ3.csv", as.is=TRUE)
env <- read.csv("../datasets/species_summaries_all.csv", as.is=TRUE)
data <- merge(growth.form, env, by.x="gs", by.y="X")
data <- data[!duplicated(data$gs),]

# Subset data to only the needed variables; recode
data <- data[,c("gs","tmin.025","pmin.025","tseas.975","pseas.975","decimallatitude.025","decimallatitude.975","gbif.records","woodiness")]
data$tmin.025 <- data$tmin.025/10
data$herb <- data$woodiness!="W"
data <- na.omit(data)

# Plot them all out
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
plot.herb("decimallatitude.025", data, by=5(-60,0), xlab=expression(paste("Minimum Southern latitude (",degree,")")))
dev.off()
pdf("../output/figures/figure_1_lat_upper.pdf")
plot.herb("decimallatitude.975", data, by=5, xlim=c(0,90), xlab=expression(paste("Maximum Northern latitude (",degree,")")))
dev.off()

# Let's see if we can get them all on the same plot...
s.data <- as.data.frame(scale(data[,2:7]))
s.data$herb <- data$herb
s.data$pmin.025 <- as.numeric(scale(log10(data$pmin.025 + 1)))

line.herb <- function(var, data, by, ...){
    bins <- with(data, seq(min(data[,var]),max(data[,var]),by=by))
    x <- with(data, table(cut(data[,var],bins), herb))[-1,]
    weight <- as.numeric(scale(rowSums(x)[-1]))
    weight <- (weight - min(weight) + 1) * 2
    x <- apply(x, 1, function(y) y[2]/sum(y))
    for(i in seq_along(x))
        lines(bins[(i+2):(i+3)], x[i:(i+1)], lwd=weight[i])
}

plot(1, type="n", xlim=c(-5,5), ylim=c(0,1), xlab="", ylab="Proportion herbaceous")
for(var in names(s.data)[1:6])
    line.herb(var, s.data, .25)
#...OK, it works, but it looks dreadful...



# Let's do some modelling

model <- glm(herb ~ tmin.025 + I(log10(pmin.025+1)) + pseas.975 + tseas.975, data=data, family=binomial, na.action="na.pass")
model.set <- dredge(model)
# Hooray! Only one model! The model wih everything in it :D
summary(model)

env.model <- glm(herb ~ decimallatitude.025 + I(decimallatitude.025^2) + decimallatitude.975 + I(decimallatitude.975^2), data=data, family=binomial)
anova(model, env.model, test="Chi")
#...uh-oh. No difference in the models...


data$decimallatitude.975.sq <- data$decimallatitude.975^2
data$decimallatitude.025.sq <- data$decimallatitude.025^2
model <- glm(herb ~ tmin.025 + I(log10(pmin.025+1)) + pseas.975 + tseas.975 + decimallatitude.025 + decimallatitude.025.sq + decimallatitude.975 + decimallatitude.975.sq, data=data, family=binomial, na.action="na.pass")
model.set <- dredge(model, subset=dc(decimallatitude.975,decimallatitude.975.sq,decimallatitude.975) & dc(decimallatitude.025,decimallatitude.025.sq,decimallatitude.025))
model.avg <- model.avg(model.set, subset=delta<4)
sink("../output/nonphylo_models/extremes_unscaled.txt")
summary(model)
sink(NULL)
#...everything matters...

# Let's look at the scaled effects to get relative importance
s.data <- as.data.frame(scale(data[,2:7]))
s.data$decimallatitude.975.sq <- s.data$decimallatitude.975^2
s.data$decimallatitude.025.sq <- s.data$decimallatitude.025^2
s.data$herb <- data$herb
s.data$pmin.025 <- as.numeric(scale(log10(data$pmin.025 + 1)))
model <- glm(herb ~ tmin.025 + pmin.025 + pseas.975 + tseas.975 + decimallatitude.025 + I(decimallatitude.025^2) + decimallatitude.975 + I(decimallatitude.975^2), data=s.data, family=binomial, na.action="na.pass")
#...OK, they're all important, but temperature is the most important, and is 3x more important than anything else
sink("../output/nonphylo_models/extremes_scaled.txt")
summary(model)
sink(NULL)


pdf("../output/figures/herb_temp.pdf")
with(data, boxplot(tmin.025 ~ herb, horizontal=TRUE, names=c("Woody",'Herbaceous'), cex=.5, pch=20, xlab=expression(paste("Minimum temperature (",degree,")"))))
model <- glm(herb ~ tmin.025 + I(log10(pmin.025+1)) + pseas.975 + tseas.975 + decimallatitude.025 + decimallatitude.025.sq + decimallatitude.975 + decimallatitude.975.sq, data=data, family=binomial, na.action="na.pass")
pred <- with(data, predict(model, data.frame(
                              tmin.025=seq(min(tmin.025),max(tmin.025),by=.01), pmin.025=mean(pmin.025), pseas.975=mean(pseas.975), tseas.975=mean(tseas.975), decimallatitude.025=mean(decimallatitude.025), decimallatitude.025.sq=mean(decimallatitude.025.sq), decimallatitude.975=mean(decimallatitude.975), decimallatitude.975.sq=mean(decimallatitude.975.sq))
                   , type="response"))
lines(pred+1 ~ seq(min(data$tmin.025), max(data$tmin.025), by=.01), lwd=3, col="red")
dev.off()
