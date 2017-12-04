# Non-phylogenetic plots and models - vessel stuff
# Will Pearse - 2017-04-30

# Headers
library(MuMIn)
library(relaimpo)

# Load data (Remove non-angiosperms (using list from TPL))
growth.form <- read.csv("../speciesTraitDataAEZ3.csv", as.is=TRUE)
growth.form <- growth.form[growth.form$woodiness=="W",]
#non.angio <- read.table("../datasets/non-angio-genera.txt", header=FALSE, as.is=TRUE)[,1]
#growth.form <- growth.form[!growth.form$g %in% non.angio,]
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
s.conduit.model <- lm(log(vesselSize) ~ phenology + tmin.025 + pmin.025 + pseas.975 + tseas.975 + decimallatitude.025 + decimallatitude.025.sq + decimallatitude.975 + decimallatitude.975.sq, data=s.data, na.action="na.pass")
s.conduit.model.set <- dredge(s.conduit.model, subset=dc(decimallatitude.975,decimallatitude.975.sq,decimallatitude.975) & dc(decimallatitude.025,decimallatitude.025.sq,decimallatitude.025))
s.conduit.model.avg <- model.avg(s.conduit.model.set, subset=delta<4)
summary(s.conduit.model.avg)
sink("../output/nonphylo_models/conduit_scaled.txt")
summary(s.conduit.model.avg)
sink(NULL)
#...all good.

# Plots
colours <- rep("forestgreen", nrow(data))
colours[data$phenology==FALSE & data$tmin.025<=0] <- "skyblue2"
colours[data$phenology==TRUE & data$tmin.025<=0] <- "darkorange"
colours[data$phenology==FALSE & data$tmin.025>0] <- "forestgreen"
colours[data$phenology==TRUE & data$tmin.025>0] <- "orange"

setEPS()
postscript("../output/figures/figure_2_temp.eps")
par(mar=c(5.1,4.6,4.1,2.1), cex.lab=1.25, cex.axis=1.25, cex=1.25)
with(data, plot(log(vesselSize) ~ tmin.025, data=data, xlab=expression(paste("Minimum Temperature (",degree,"C)")), ylab=expression(paste("Log Conduit Area (", mm^2, ")")), type="n"))
abline(h=log(0.001520531), col="grey50", lwd=3, lty=2)
with(data, points(log(vesselSize) ~ tmin.025, data=data, pch=20, col=colours))
text(-35, -10, "(c)", font=2)
dev.off()
setEPS()
postscript("../output/figures/figure_2_precip.eps")
par(mar=c(5.1,4.6,4.1,2.1), cex.lab=1.25, cex.axis=1.25, cex=1.25)
with(data, plot(log(vesselSize) ~ I(log(pmin.025+exp(1))), data=data, xlab="Log Minimum Precipitation (log[mm])", ylab=expression(paste("Log Conduit Area (", mm^2, ")")), type="n"))
abline(h=log(0.001520531), col="grey50", lwd=3, lty=2)
with(data, points(log(vesselSize) ~ I(log(pmin.025+exp(1))), data=data, pch=20, col=colours))
text(1.25, -10, "(d)", font=2)
dev.off()
setEPS()
postscript("../output/figures/figure_2_temp_seas.eps")
par(mar=c(5.1,4.6,4.1,2.1), cex.lab=1.25, cex.axis=1.25, cex=1.25)
with(data, plot(log(vesselSize) ~ tseas.975, data=data, pch=21, bg=ifelse(phenology, "black", "white"), xlab="Temperature Seasonality", ylab=expression(paste("Log Conduit Area (", mm^2, ")")), type="n"))
abline(h=log(0.001520531), col="grey50", lwd=3, lty=2)
with(data, points(log(vesselSize) ~ tseas.975, data=data, pch=20, col=colours))
text(1000, -10, "(e)", font=2)
dev.off()
setEPS()
postscript("../output/figures/figure_2_precip_seas.eps")
par(mar=c(5.1,4.6,4.1,2.1), cex.lab=1.25, cex.axis=1.25, cex=1.25)
with(data, plot(log(vesselSize) ~ pseas.975, data=data, pch=21, bg=ifelse(phenology, "black", "white"), xlab="Precipitation Seasonality", ylab=expression(paste("Log Conduit Area (", mm^2, ")")), type="n"))
abline(h=log(0.001520531), col="grey50", lwd=3, lty=2)
with(data, points(log(vesselSize) ~ pseas.975, data=data, pch=20, col=colours))
text(20, -10, "(f)", font=2)
dev.off()
setEPS()
postscript("../output/figures/figure_2_lat_lower.eps")
par(mar=c(5.1,4.6,4.1,2.1), cex.lab=1.25, cex.axis=1.25, cex=1.25)
with(data, plot(log(vesselSize) ~ decimallatitude.025, data=data, xlab=expression(paste("Minimum Latitude (",degree,")")), ylab=expression(paste("Log Conduit Area (", mm^2, ")")), type="n", xlim=c(-60,0)))
legend("topleft", legend=c("Deciduous/Freezing", "Deciduous/Non-freezing"), pch=20, col=c("orange","darkorange","skyblue2","forestgreen"), bty="n")
legend("bottomright", legend=c("Evergreen/Freezing","Evergreen/Non-freezing"), pch=20, col=c("skyblue2","forestgreen"), bty="n")
abline(h=log(0.001520531), col="grey50", lwd=3, lty=2)
with(data, points(log(vesselSize) ~ decimallatitude.025, data=data, pch=20, col=colours))
text(-55, -10, "(a)", font=2)
dev.off()
setEPS()
postscript("../output/figures/figure_2_lat_upper.eps")
par(mar=c(5.1,4.6,4.1,2.1), cex.lab=1.25, cex.axis=1.25, cex=1.25)
with(data, plot(log(vesselSize) ~ decimallatitude.975, data=data, xlab=expression(paste("Maximum Latitude (",degree,")")), ylab=expression(paste("Log Conduit Area (", mm^2, ")")), type="n", xlim=c(0,70)))
abline(h=log(0.001520531), col="grey50", lwd=3, lty=2)
with(data, points(log(vesselSize) ~ decimallatitude.975, data=data, pch=20, col=colours))
text(5, -10, "(b)", font=2)
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

# ...variance partitioning (ish)...
library(relaimpo)
model <- lm(log(vesselSize) ~ factor(phenology) + tmin.025 + pmin.025 + pseas.975 + tseas.975 + decimallatitude.025 + decimallatitude.025.sq + decimallatitude.975 + decimallatitude.975.sq, data=s.data, na.action="na.pass")
calc.relimp(model)
