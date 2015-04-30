#Headers; load workspace
require(phylolm)
source("workspace.R")

#Silly models!
c.data$data$evergreen <- ifelse(c.data$data$phenology == "EV", 1, 0)
model <- phyloglm(evergreen ~ log(vesselSize) * pole.lim, data=c.data$data, phy=c.data$phy, start.alpha=0.001)
summmary(model)
#...I'm not sure I trust this model...
