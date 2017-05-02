# Non-phylogenetic plots and models - vessel stuff
# Will Pearse - 2017-04-30

# Headers
library(MuMIn)
library(relaimpo)

# Load data
growth.form <- read.csv("../speciesTraitDataAEZ3.csv", as.is=TRUE)
growth.form <- growth.form[growth.form$woodiness=="W",]
growth.form <- growth.form[,c("gs","phenology","vesselSize")]
env <- read.csv("../datasets/species_summaries.csv", as.is=TRUE)
data <- merge(growth.form, env, by.x="gs", by.y="X")
data <- data[!duplicated(data$gs),]

# Subset data to only the needed variables; rescale/re-code
data <- data[,c("gs","tmin.025","pmin.025","tseas.975","pseas.975","gbif.records","vesselSize","phenology")]
data$tmin.025 <- data$tmin.025*10
data$phenology <- as.numeric(data$phenology=="D")

# Conduit diameter modelling
conduit.model <- lm(log(vesselSize) ~ tmin.025 + pmin.025 + pseas.975 + tseas.975 + phenology, data=data, na.action="na.pass")
conduit.model.set <- dredge(conduit.model)
conduit.model.avg <- model.avg(conduit.model.set, subset=delta<4)
summary(conduit.model.avg)
#...phenology, precip seasonality, and minimum temperature (phew). Let's plot...

with(data, plot(log(vesselSize) ~ tmin.025))
with(data, plot(log(vesselSize) ~ pseas.975))
with(data, boxplot(log(vesselSize) ~ phenology))
#...so something could be significant but not important. Do a bit of scaling and relaimpo

# Scaled conduit diameter models
orig.data <- data
data$tmin.025 <- as.numeric(scale(data$tmin.025))
data$pmin.025 <- as.numeric(scale(data$pmin.025))
data$tseas.975 <- as.numeric(scale(data$tseas.975))
data$pseas.975 <- as.numeric(scale(data$pseas.975))
conduit.model <- lm(log(vesselSize) ~ tmin.025 + pmin.025 + pseas.975 + tseas.975 + phenology, data=data, na.action="na.pass")
conduit.model.set <- dredge(conduit.model)
conduit.model.avg <- model.avg(conduit.model.set, subset=delta<4)
summary(conduit.model.avg)
#...effect of minimum temperature is ~2x the effect of seasonality of precipitation

# Scaled conduit diameter models
conduit.model <- calc.relimp(log(vesselSize) ~ tmin.025 + pmin.025 + pseas.975 + tseas.975 + phenology, data=data)
conduit.model
#...the effect isn't as big here, but it's still something. Clearly temperature matters more.

# Plots
pdf("../output/figures/figure_2_temp.pdf")
with(data, plot(log(vesselSize) ~ tmin.025, data=data, pch=21, bg=ifelse(phenology, "black", "white"), xlab=expression(paste("Minimum temperature (",degree,"C)")), ylab="Log Conduit Diameter", cex=1.5))
dev.off()
pdf("../output/figures/figure_2_precip.pdf")
with(data, plot(log(vesselSize) ~ pmin.025, data=data, pch=21, bg=ifelse(phenology, "black", "white"), xlab="Minimum precipitation (mm)", ylab="Log Conduit Diameter", cex=1.5))
dev.off()
pdf("../output/figures/figure_2_temp_seas.pdf")
with(data, plot(log(vesselSize) ~ tseas.975, data=data, pch=21, bg=ifelse(phenology, "black", "white"), xlab="Temperature seasonality", ylab="Log Conduit Diameter", cex=1.5))
dev.off()
pdf("../output/figures/figure_2_precip_seas.pdf")
with(data, plot(log(vesselSize) ~ pseas.975, data=data, pch=21, bg=ifelse(phenology, "black", "white"), xlab="Precipitation seasonality", ylab="Log Conduit Diameter", cex=1.5))
dev.off()
pdf("../output/figures/figure_2_lat_lower.pdf")
with(data, plot(log(vesselSize) ~ tmin.025, data=data, pch=21, bg=ifelse(phenology, "black", "white"), xlab=expression(paste("Minimum temperature (",degree,"C)")), ylab="Log Conduit Diameter", cex=1.5))
dev.off()
pdf("../output/figures/figure_2_lat_upper.pdf")
with(data, plot(log(vesselSize) ~ tmin.025, data=data, pch=21, bg=ifelse(phenology, "black", "white"), xlab=expression(paste("Minimum temperature (",degree,"C)")), ylab="Log Conduit Diameter", cex=1.5))
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
