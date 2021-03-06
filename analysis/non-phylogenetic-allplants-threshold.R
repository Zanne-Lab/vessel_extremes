# Non-phylogenetic plots and models
# - now with only the top 50% species
# Will Pearse - 2017-04-30

library(MuMIn)

# Load data
growth.form <- read.csv("../speciesTraitDataAEZ3.csv", as.is=TRUE)
env <- read.csv("../datasets/species_summaries_all.csv", as.is=TRUE)
data <- merge(growth.form, env, by.x="gs", by.y="X")
data <- data[!duplicated(data$gs),]

# Remove non-angiosperms (using list from TPL)
non.angio <- read.table("../datasets/non-angio-genera.txt", header=FALSE, as.is=TRUE)[,1]
data <- data[!data$g %in% non.angio,]

# Subset data to only have those in the top half of samples on GBIF
data <- data[data$gbif.record >= median(data$gbif.records),]

# Subset data to only the needed variables; recode
data <- data[,c("gs","tmin.025","pmin.025","tseas.975","pseas.975","decimallatitude.025","decimallatitude.975","gbif.records","woodiness")]
data$tmin.025 <- data$tmin.025/10
data$herb <- data$woodiness!="W"
data <- na.omit(data)

# Plot them all out
plot.herb <- function(var, data, by, ...){
    bins <- with(data, seq(min(data[,var]),max(data[,var]),by=by))
    x <- with(data, table(cut(data[,var],bins), herb))
    #weight <- as.numeric(scale(rowSums(x)[-1]))
    #weight <- (weight - min(weight) + 1) * 2
    weight <- (log(rowSums(x)[-1]) - 1) * 3
    x <- apply(x, 1, function(y) y[2]/sum(y))
    plot(x ~ I(bins[-length(bins)]+by/2), axes=FALSE, ylab="Proportion Herbaceous", ylim=c(0,1), lwd=weight, type="n", ...)
    for(i in seq(1,length(x)-1))
        lines((bins[-length(bins)]+by/2)[i:(i+1)], x[c(i,i+1)], lwd=weight[i])
    axis(1)
    axis(2)
}
setEPS()
postscript("../output/figures/figure_1_temp_top50.eps")
par(mar=c(5.1,4.6,4.1,2.1), cex.lab=1.25, cex.axis=1.25, cex=1.25)
plot.herb("tmin.025", data, by=5, xlab=expression(paste("Minimum Temperature (",degree,")")))
text(-35, 1, "(c)", font=2, cex=2)
dev.off()
setEPS()
postscript("../output/figures/figure_1_precip_top50.eps")
par(mar=c(5.1,4.6,4.1,2.1), cex.lab=1.25, cex.axis=1.25, cex=1.25)
data$log.pmin.025 <- log(data$pmin.025+exp(1))
plot.herb("log.pmin.025", data, by=.5, xlab="Log Minimum Precipitation (log[mm])", xlim=c(1,6))
text(1.5, 1, "(d)", font=2, cex=2)
dev.off()
setEPS()
postscript("../output/figures/figure_1_temp_seas_top50.eps")
par(mar=c(5.1,4.6,4.1,2.1), cex.lab=1.25, cex.axis=1.25, cex=1.25)
plot.herb("tseas.975", data, by=500, xlim=c(0,15000), xlab="Temperature Seasonality")
text(1000, 1, "(e)", font=2, cex=2)
dev.off()
setEPS()
postscript("../output/figures/figure_1_precip_seas_top50.eps")
par(mar=c(5.1,4.6,4.1,2.1), cex.lab=1.25, cex.axis=1.25, cex=1.25)
plot.herb("pseas.975", data, by=20, xlim=c(0,200), xlab="Precipitation Seasonality")
text(10, 1, "(f)", font=2, cex=2)
dev.off()
setEPS()
postscript("../output/figures/figure_1_lat_lower_top50.eps")
par(mar=c(5.1,4.6,4.1,2.1), cex.lab=1.25, cex.axis=1.25, cex=1.25)
plot.herb("decimallatitude.025", data, by=5, xlim=c(-60,0), xlab=expression(paste("Minimum Latitude (",degree,")")))
legend("topright", legend=c(100,1000,10000), lwd=log(c(100,1000,10000)-1)*3, bty="n", title="Number of species", horiz=TRUE)
text(-55, 1, "(a)", font=2, cex=2)
dev.off()
setEPS()
postscript("../output/figures/figure_1_lat_upper_top50.eps")
par(mar=c(5.1,4.6,4.1,2.1), cex.lab=1.25, cex.axis=1.25, cex=1.25)
plot.herb("decimallatitude.975", data, by=5, xlim=c(0,90), xlab=expression(paste("Maximum Latitude (",degree,")")))
text(5, 1, "(b)", font=2, cex=2)
dev.off()

# Let's see if we can get them all on the same plot...
s.data <- as.data.frame(scale(data[,2:7]))
s.data$herb <- data$herb
s.data$pmin.025 <- as.numeric(scale(log(data$pmin.025 + 1)))

line.herb <- function(var, data, by, ...){
    bins <- with(data, seq(min(data[,var]),max(data[,var]),by=by))
    x <- with(data, table(cut(data[,var],bins), herb))[-1,]
    weight <- as.numeric(scale(rowSums(x)[-1]))
    weight <- (weight - min(weight) + 1) * 2
    x <- apply(x, 1, function(y) y[2]/sum(y))
    for(i in seq_along(x))
        lines(bins[(i+2):(i+3)], x[i:(i+1)], lwd=weight[i])
}

plot(1, type="n", xlim=c(-5,5), ylim=c(0,1), xlab="", ylab="Proportion Hherbaceous")
for(var in names(s.data)[1:6])
    line.herb(var, s.data, .25)
#...OK, it works, but it looks dreadful...



# Let's do some modelling
# Variation partitioning ---------
# Using McFadden's adjusted psuedo R2
# https://stats.idre.ucla.edu/other/mult-pkg/faq/general/faq-what-are-pseudo-r-squareds/
adjR2 = function(m)  {
    1 - ((m$deviance - length(m$coef) + 1) / m$null.deviance)
}

full <- glm(herb ~ tmin.025 + I(log(pmin.025+1)) + pseas.975 + tseas.975 +
            decimallatitude.025 + I(decimallatitude.025^2) +
            decimallatitude.975 + I(decimallatitude.975^2),
            data=data, family=binomial, na.action="na.pass")
clim <- glm(herb ~ tmin.025 + I(log(pmin.025+1)) + pseas.975 + tseas.975,
            data=data, family=binomial, na.action="na.pass")
geo <- glm(herb ~ decimallatitude.025 + I(decimallatitude.025^2) +
           decimallatitude.975 + I(decimallatitude.975^2),
           data=data, family=binomial, na.action="na.pass")


clim_only = adjR2(full) - adjR2(geo)
geo_only = adjR2(full) - adjR2(clim)
shared = adjR2(full) - clim_only - geo_only
unexp = 1 - adjR2(full)

vegan::showvarparts(2)
clim_only   #[a] = 0.03 
shared      #[b] = 0.12 
geo_only    #[c] = 0.008 
unexp       #[d] = 0.83

# AIC dredging --------------
model <- glm(herb ~ tmin.025 + I(log(pmin.025+1)) + pseas.975 + tseas.975, data=data, family=binomial, na.action="na.pass")
model.set <- dredge(model)
# Hooray! Only one model! The model wih everything in it :D
summary(model)

env.model <- glm(herb ~ decimallatitude.025 + I(decimallatitude.025^2) + decimallatitude.975 + I(decimallatitude.975^2), data=data, family=binomial)
anova(model, env.model, test="Chi")
#...uh-oh. No difference in the models...


data$decimallatitude.975.sq <- data$decimallatitude.975^2
data$decimallatitude.025.sq <- data$decimallatitude.025^2
model <- glm(herb ~ tmin.025 + I(log(pmin.025+1)) + pseas.975 + tseas.975 + decimallatitude.025 + decimallatitude.025.sq + decimallatitude.975 + decimallatitude.975.sq, data=data, family=binomial, na.action="na.pass")
model.set <- dredge(model, subset=dc(decimallatitude.975,decimallatitude.975.sq,decimallatitude.975) & dc(decimallatitude.025,decimallatitude.025.sq,decimallatitude.025))
model.avg <- model.avg(model.set, subset=delta<4)
sink("../output/nonphylo_models/extremes_unscaled_top50.txt")
summary(model.avg)
sink(NULL)
#...roughly the same as what we had before, but the seasonality terms
#   have dropped. Probably some information in that, and (if anything)
#   drives home that we're looking at the right things in the
#   macro-evolutionary stuff

# Let's look at the scaled effects to get relative importance
s.data <- as.data.frame(scale(data[,2:7]))
s.data$decimallatitude.975.sq <- s.data$decimallatitude.975^2
s.data$decimallatitude.025.sq <- s.data$decimallatitude.025^2
s.data$herb <- data$herb
s.data$pmin.025 <- as.numeric(scale(log(data$pmin.025 + 1)))
model <- glm(herb ~ tmin.025 + pmin.025 + pseas.975 + tseas.975 + decimallatitude.025 + I(decimallatitude.025^2) + decimallatitude.975 + I(decimallatitude.975^2), data=s.data, family=binomial, na.action="na.pass")

model.set <- dredge(model, subset=dc(decimallatitude.975,decimallatitude.975.sq,decimallatitude.975) & dc(decimallatitude.025,decimallatitude.025.sq,decimallatitude.025))
model.avg <- model.avg(model.set, subset=delta<4)

#...exactly the same take-home as above, which is (if nothing else) re-assuring
sink("../output/nonphylo_models/extremes_scaled_top50.txt")
summary(model.avg)
sink(NULL)

# ...variance partitioning (ish)...
library(relaimpo)
model <- glm(herb ~ tmin.025 + pmin.025 + pseas.975 + tseas.975 + decimallatitude.025 + I(decimallatitude.025^2) + decimallatitude.975 + I(decimallatitude.975^2), data=s.data, na.action="na.pass")
calc.relimp(model)
