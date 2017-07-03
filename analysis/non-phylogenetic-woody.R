# Non-phylogenetic plots and models - vessel stuff
# Will Pearse - 2017-04-30

# Headers
library(MuMIn)
library(relaimpo)

# Load data
growth.form <- read.csv("../speciesTraitDataAEZ3.csv", as.is=TRUE)
growth.form <- growth.form[growth.form$woodiness=="W",]
growth.form <- growth.form[,c("gs","phenology","vesselSize")]
env <- read.csv("../datasets/species_summaries_all.csv", as.is=TRUE)
data <- merge(growth.form, env, by.x="gs", by.y="X")
data <- data[!duplicated(data$gs),]

# Subset data to only the needed variables; rescale/re-code
data <- data[,c("gs","tmin.025","pmin.025","tseas.975","pseas.975","gbif.records","decimallatitude.025","decimallatitude.975","vesselSize","phenology")]
data$tmin.025 <- data$tmin.025/10
data$phenology <- as.numeric(data$phenology=="D")
data <- na.omit(data)

# Conduit diameter modelling
data$decimallatitude.975.sq <- data$decimallatitude.975^2
data$decimallatitude.025.sq <- data$decimallatitude.025^2
conduit.model <- lm(log(vesselSize) ~ tmin.025 + pmin.025 + pseas.975 + tseas.975 + phenology + decimallatitude.025 + decimallatitude.025.sq + decimallatitude.975 + decimallatitude.975.sq, data=data, na.action="na.pass")
conduit.model.set <- dredge(conduit.model, subset=dc(decimallatitude.975,decimallatitude.975.sq,decimallatitude.975) & dc(decimallatitude.025,decimallatitude.025.sq,decimallatitude.025))
conduit.model.avg <- model.avg(conduit.model.set, subset=delta<4)
summary(conduit.model.avg)
sink("../output/nonphylo_models/conduit_unscaled.txt")
summary(conduit.model.avg)
sink(NULL)
#...looks fine. Let's get relative importance...

# Scaled conduit diameter models
s.data <- as.data.frame(scale(data[,c("tmin.025", "pmin.025", "tseas.975", "pseas.975", "decimallatitude.025", "decimallatitude.975")]))
s.data$phenology <- data$phenology
s.data$vesselSize <- data$vesselSize
s.data$decimallatitude.975.sq <- s.data$decimallatitude.975^2
s.data$decimallatitude.025.sq <- s.data$decimallatitude.025^2
s.conduit.model <- lm(log(vesselSize) ~ tmin.025 + pmin.025 + pseas.975 + tseas.975 + phenology + decimallatitude.025 + decimallatitude.025.sq + decimallatitude.975 + decimallatitude.975.sq, data=s.data, na.action="na.pass")
s.conduit.model.set <- dredge(s.conduit.model, subset=dc(decimallatitude.975,decimallatitude.975.sq,decimallatitude.975) & dc(decimallatitude.025,decimallatitude.025.sq,decimallatitude.025))
s.conduit.model.avg <- model.avg(s.conduit.model.set, subset=delta<4)
summary(s.conduit.model.avg)
sink("../output/nonphylo_models/conduit_scaled.txt")
summary(s.conduit.model.avg)
sink(NULL)
#...all good.

# Plots

pdf("../output/figures/figure_2_temp.pdf")
par(mar=c(5.1,4.6,4.1,2.1), cex.lab=1.25, cex.axis=1.25, cex=1.25)
with(data, plot(log(vesselSize) ~ tmin.025, data=data, bg=ifelse(phenology, "black", "white"), xlab=expression(paste("Minimum temperature (",degree,"C)")), ylab="Log Conduit Area", type="n"))
abline(h=log(0.001520531), col="grey50", lwd=3, lty=2)
with(data, points(log(vesselSize) ~ tmin.025, data=data, pch=21, bg=ifelse(phenology, "black", "white")))
dev.off()

pdf("../output/figures/figure_2_precip.pdf")
par(mar=c(5.1,4.6,4.1,2.1), cex.lab=1.25, cex.axis=1.25, cex=1.25)
with(data, plot(log(vesselSize) ~ I(log(pmin.025+exp(1))), data=data, xlab="Log minimum precipitation (log[mm])", ylab="Log Conduit Area", type="n"))
fg <- bg <- ifelse(data$tmin.025>0, "red", "blue")
bg <- ifelse(data$phenology, bg, "white")
abline(h=log(0.001520531), col="grey50", lwd=3, lty=2)
with(data, points(log(vesselSize) ~ I(log(pmin.025+exp(1))), data=data, pch=21, col=fg, bg=bg))
dev.off()

pdf("../output/figures/figure_2_temp_seas.pdf")
par(mar=c(5.1,4.6,4.1,2.1), cex.lab=1.25, cex.axis=1.25, cex=1.25)
with(data, plot(log(vesselSize) ~ tseas.975, data=data, pch=21, bg=ifelse(phenology, "black", "white"), xlab="Temperature seasonality", ylab="Log Conduit Area", type="n"))
abline(h=log(0.001520531), col="grey50", lwd=3, lty=2)
with(data, points(log(vesselSize) ~ tseas.975, data=data, pch=21, bg=ifelse(phenology, "black", "white")))
dev.off()

pdf("../output/figures/figure_2_precip_seas.pdf")
par(mar=c(5.1,4.6,4.1,2.1), cex.lab=1.25, cex.axis=1.25, cex=1.25)
with(data, plot(log(vesselSize) ~ pseas.975, data=data, pch=21, bg=ifelse(phenology, "black", "white"), xlab="Precipitation seasonality", ylab="Log Conduit Area", type="n"))
abline(h=log(0.001520531), col="grey50", lwd=3, lty=2)
with(data, points(log(vesselSize) ~ pseas.975, data=data, pch=21, bg=ifelse(phenology, "black", "white")))
dev.off()

pdf("../output/figures/figure_2_lat_lower.pdf")
par(mar=c(5.1,4.6,4.1,2.1), cex.lab=1.25, cex.axis=1.25, cex=1.25)
with(data, plot(log(vesselSize) ~ decimallatitude.025, data=data, pch=21, bg=ifelse(phenology, "black", "white"), xlab=expression(paste("Minimum latitude (",degree,")")), ylab="Log Conduit Area", type="n", xlim=c(-60,0)))
abline(h=log(0.001520531), col="grey50", lwd=3, lty=2)
with(data, points(log(vesselSize) ~ decimallatitude.025, data=data, pch=21, bg=ifelse(phenology, "black", "white")))
dev.off()

pdf("../output/figures/figure_2_lat_upper.pdf")
par(mar=c(5.1,4.6,4.1,2.1), cex.lab=1.25, cex.axis=1.25, cex=1.25)
with(data, plot(log(vesselSize) ~ decimallatitude.975, data=data, pch=21, bg=ifelse(phenology, "black", "white"), xlab=expression(paste("Maximum latitude (",degree,")")), ylab="Log Conduit Area", type="n", xlim=c(0,70)))
abline(h=log(0.001520531), col="grey50", lwd=3, lty=2)
with(data, points(log(vesselSize) ~ decimallatitude.975, data=data, pch=21, bg=ifelse(phenology, "black", "white")))
dev.off()




#...just to check...
# Phenology modelling
phenol.model <- glm(phenology ~ tmin.025 + pmin.025 + pseas.975 + tseas.975, data=data, family=binomial, na.action="na.pass")
phenol.model.set <- dredge(phenol.model)
phenol.model.avg <- model.avg(phenol.model.set, subset=delta<4)
summary(phenol.model.avg)

test.model <- glm(phenology ~ pmin.025 + tseas.975, data=data)
plt.mod(test.model)
with(data, boxplot(tseas.975 ~ phenology))
#...I think the seasonality effect passes the sense check - it looks like a real thing to me...
