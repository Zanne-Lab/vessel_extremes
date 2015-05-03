## 3 strategies paper matching tree and data
require(bayou)
require(OUwie)
require(aRbor)
setwd("~/repos/growthforms/")

tree <- read.tree("./Vascular_Plants_rooted.dated.tre")
AEZ3 <- read.csv("./speciesTraitDataAEZ3.csv")
climate <- read.csv("./WoodyClimateHemi.12May.csv")
AEZ3$species <- gsub(" ", "_", AEZ3$gs)
climate$species <- gsub(" ", "_", climate$species)
dat <- left_join(climate, AEZ3)

#td <- make.treedata(tree, dat)
#saveRDS(td, file="./output/cleandata/matchedTankAEZ3Clim.rds")
td <- readRDS(file="./output/cleandata/matchedTankAEZ3Clim.rds")

tdOUwie <- filter(td, !is.na(vesselSize), !is.na(phenology),!is.na(tmin.lo))
tdOUwie <- mutate(tdOUwie, freeze= (tmin.lo < 0), lnVs=log(vesselSize))
tdOUwie <- mutate(tdOUwie,OU1= "global", OU2phenology= phenology, OU2freeze= c("NF", "F")[as.numeric(freeze)+1],  OU4= paste(phenology, c("NF", "F")[as.numeric(freeze)+1], sep=""))
tdOUwie <- treeply(tdOUwie, reorder, "postorder")
tdOUwie <- select(tdOUwie, lnVs,OU1, OU2phenology, OU2freeze, OU4, tmin.lo)
boxplot(tdOUwie$dat$lnVs~ tdOUwie$dat$OU4)
boxplot(tdOUwie$dat$lnVs~ tdOUwie$dat$OU2phenology)
boxplot(tdOUwie$dat$lnVs~ tdOUwie$dat$OU2freeze)
plot(tdOUwie$dat$tmin.lo, tdOUwie$dat$lnVs, pch=21, bg=as.numeric(tdOUwie$dat$OU4))

require(diversitree)
require(OUwie)
tree <- tdOUwie$phy
dat <- tdOUwie$dat

lik <- list()
lik[["OU2ph"]] <- make.mk2(tree, states=setNames(as.numeric(factor(dat$OU2phenology)), tree$tip.label)-1)
lik[["OU2fr"]] <- make.mk2(tree, states=setNames(as.numeric(factor(dat$OU2freeze)), tree$tip.label)-1)
lik[["OU4"]] <- make.mkn(tree, states=setNames(as.numeric(factor(dat$OU4)), tree$tip.label), k=4)
lik[["cOU2ph"]] <- constrain(lik$OU2ph, q01 ~ q10)
lik[["cOU2fr"]] <- constrain(lik$OU2fr, q01 ~ q10)
lik[["sOU4"]] <- constrain(lik$OU4, q12 ~ q21, q13~q31, q14~q41, q23~q32, q24~q42, q34~q43)
lik[["cOU4"]] <- constrain(lik$OU4, q12 ~ q21, q13~q31, q14~q41, q23~q32, q24~q42, q34~q43, q12~q31, q12~q41, q12~q32, q12~q42, q12~q43, q21~q31, q32~q41, q42~q43, q31~q41, q31~q43, q41~q43)
mkpars <- list(c(0.1, 0.1), c(0.1, 0.1), rep(0.1, 12), c(0.1), c(0.1), rep(0.1, 6), 0.1) 
npars <- sapply(mkpars, length)

fits <- lapply(1:6,function(x) find.mle(lik[[x]], mkpars[[x]]))
fits[[3]]
aics <- 2*npars[1:6]-2*sapply(fits, function(x) x$lnLik)
names(aics) <- names(lik)[1:6]
aics

#tmp <- mcmc(lik$OU4, fits[[3]]$par, nsteps=100, w=.1)
#w <- median(diff(sapply(tmp[2:12], quantile, c(.05, .95))))
#postOU4 <- mcmc(lik$OU4, x.init=fits[[3]]$par, nsteps = 1000, w=w)
#Qs <- postOU4[sample(40:115, 10),-c(1,14)]
makeQ <- function(qpars){
  nQ <- names(qpars)
  inds <- lapply(strsplit(sapply(strsplit(nQ, 'q'), function(x)x[2]),""), as.numeric)
  Q <- matrix(NA, ncol=4, nrow=4)
  for(i in 1:length(inds)){
    Q[inds[[i]][1], inds[[i]][2]] <- as.vector(unlist(qpars[i]))
  }
  diag(Q) <--1* apply(Q, 1, sum, na.rm=TRUE)
  colnames(Q) <- c("DF", "DNF", "EVF", "EVNF")
  rownames(Q) <- c("DF", "DNF", "EVF", "EVNF")
  Q
}
Q <- makeQ(fits[[3]]$par)


## ARD model better for all...
smtrees <- list()
tree$node.label <- NULL
smtrees[["OU2ph"]] <- make.simmap(tree, setNames(dat$OU2phenology, tree$tip.label), model="ARD", nsim=10)
smtrees[["OU2fr"]] <- make.simmap(tree, setNames(dat$OU2freeze, tree$tip.label), model="ARD", nsim=10)
smtrees[["OU4"]] <- make.simmap(tree, setNames(dat$OU4, tree$tip.label), model="ARD", nsim=10, Q=Q)
plotSimmap(smtrees$OU4)

models <- c("BM1","BMS","OU1","OUM","OUMV")

#plot(tree, show.tip.label=FALSE)
#tiplabels(pch=21, bg=as.numeric(factor(dat$OU4)), col=as.numeric(factor(dat$OU4)), cex=0.5)

#phenogram(tree, setNames(dat$lnVs, tree$tip.label))
#phenogram(tdOUwie$phy,setNames( tdOUwie$dat$tmin.lo, tdOUwie$phy$tip.label))
#ouFitTmin <- fitContinuous(td$phy, setNames(td$dat$tmin.lo, td$phy$tip.label), model="OU")

#bmtmin.lo <- fitContinuous(tdOUwie$phy, setNames(tdOUwie$dat$tmin.lo, tdOUwie$phy$tip.label), model="BM")
#outmin.lo <- fitContinuous(tdOUwie$phy, setNames(tdOUwie$dat$tmin.lo, tdOUwie$phy$tip.label), model="OU")
#bmtmin.lo$opt$aic
#outmin.lo$opt$aic
#outmin.lo
#pars.ouTemp <- list(alpha=outmin.lo$opt$alpha, sig2=outmin.lo$opt$sigsq,k=0, ntheta=1, theta=outmin.lo$opt$z0, sb=numeric(0), loc=numeric(0), t2=numeric(0))

#tdcw <- treeply(tdOUwie,reorder, "postorder")
#ntaxa <- length(tdcw$phy$tip.label)
#nknown <- length(tdcw$phy$tip.label)
#asrdat <- rep(NA, ntaxa+ntaxa-1)
#asrdat[1:ntaxa] <- tdcw$dat$tmin.lo
#asrdat[ntaxa+1] <- pars.ouTemp$theta[1]

## I don't think this is working yet....stochastic OU reconstruction of ancestral temperature states
vcvOU <- function (phy, dat, pars, SE, internal = TRUE) {
  new.pars <- pars
  ntips <- length(phy$tip.label)
  sig2 <- new.pars$sig2
  alpha <- new.pars$alpha
  D <- dist.nodes(phy)
  Cii <- D[ntips + 1, ]
  C <- D
  C[, ] <- 0
  for (i in 1:nrow(D)) for (j in 1:ncol(D)) C[i, j] <- sig2/(2 * 
                                                               alpha) * exp(-alpha * D[i, j]) * (1 - exp(-2 * alpha * 
                                                                                                           (Cii[i] + Cii[j] - D[i, j])))
  diag(C) <- diag(C) + 1e-10
  if(new.pars$ntheta > 1){
    W <- .allnodes.W(cache, new.pars)
    ExpV <- W %*% new.pars$theta
  } else {
    ExpV <- rep(new.pars$theta, length(C[,1]))
  }
  return(list(ExpV = ExpV, VCV = C))
}
simConditionalOU <- function(tree, dat, pars, plot=TRUE, ...){
  ntaxa <- length(tree$tip.label)
  nnodes <- ntaxa+tree$Nnode
  V <- vcvOU(tree, dat, pars, SE=0)$VCV
  mu <- rep(pars$theta[1], nnodes)
  root <- dat[ntaxa+1]
  X <- dat[-(ntaxa+1)]
  unknown <- which(is.na(X))
  known <- which(!is.na(X))
  muk <- mu[known]
  muu <- mu[unknown]
  Vkk <- V[known, known]
  Vuu <- V[unknown, unknown]
  Vku <- V[known, unknown]
  Vuk <- V[unknown, known]
  
  mubar <- muu + Vuk%*%solve(Vkk)%*%(X[known]-muk)
  sigmabar <- Vuu - Vuk%*%solve(Vkk)%*%Vku
  res <- MASS::mvrnorm(1, mubar, sigmabar)
  X[unknown] <- res
  Sim <- c(X[1:ntaxa], root, X[(ntaxa+1):length(X)])
  names(Sim) <- c(tree$tip.label, ntaxa+1, (ntaxa+2):length(dat))
  if(plot){
    phenogram(tree, Sim,...)
  }
  return(list(dat=Sim, mu=mubar, sigma2=sigmabar))
}
#ancTemps <- list()
#for(i in 1:10){
#  pars <- pars.ouTemp
  #pars$theta <- rnorm(1, pars.ouTemp$theta, sqrt(pars$sig2/(2*pars$alpha)))
#  ancTemps[[i]] <- simConditionalOU(tdcw$phy, asrdat, pars, plot=TRUE)
#}
#ancTemps <- lapply(ancTemps, function(x) x$dat)
#ancFreeze <- lapply(ancTemps, function(x) as.numeric(x <= 0))
#plot(tree, show.tip.label=FALSE, type="fan")
#ntips <- length(tdOUwie$phy$tip.label)
#nodelabels(pch=21, bg=ancFreeze[[1]][-(1:ntips)]+1)
#tiplabels(pch=21, bg=ancFreeze[[1]][1:ntips]+1, cex=0.5, col =ancFreeze[[1]][1:ntips]+1)


#asrTemp <- ace(setNames(tdOUwie$dat$tmin.lo, tdOUwie$phy$tip.label), tree, "continuous")
OUwie.OU4 <- data.frame('Genus_species'=tdOUwie$phy$tip.label, Reg=tdOUwie$dat$OU4, X=tdOUwie$dat$lnVs)
OU4.fits <- lapply(models, function(x) lapply(smtrees$OU4, function(y) OUwie(y, OUwie.OU4, model = x, simmap.tree=TRUE)) )

OUwie.OU2ph <- data.frame('Genus_species'=tdOUwie$phy$tip.label, Reg=tdOUwie$dat$OU2ph, X=tdOUwie$dat$lnVs)
OU2ph.fits <- lapply(models, function(x) lapply(smtrees$OU2ph, function(y) OUwie(y, OUwie.OU2ph, model = x, simmap.tree=TRUE)) )

OUwie.OU2fr <- data.frame('Genus_species'=tdOUwie$phy$tip.label, Reg=tdOUwie$dat$OU2fr, X=tdOUwie$dat$lnVs)
OU2fr.fits <- lapply(models, function(x) lapply(smtrees$OU2fr, function(y) OUwie(y, OUwie.OU2fr, model = x, simmap.tree=TRUE)) )

pdf("../output/figures/OUwieFreezeAnalysis.pdf")
fits <- list(OU4=OU4.fits, OU2ph=OU2ph.fits, OU2fr=OU2fr.fits)
aics <- lapply(fits, function(y) lapply(y, function(x) sapply(x, function(z) z$AIC)))
aics <- lapply(aics, function(x) do.call(cbind, x))
aictable <- do.call(cbind, aics)
colnames(aictable) <- unlist(lapply(1:length(fits), function(x) paste(names(fits)[x], models, sep=".")))
par(mar=c(8, 4, 2, 1))
boxplot(aictable[,order(apply(aictable,2,mean))], las=2)

extractModel <- function(model, fits){
  aics <-  lapply(fits, function(x) sapply(x, function(y) y$AIC))
  aics <- do.call(cbind, aics)
  fits <- fits[!is.na(aics[,model])]
  thetaests <- lapply(fits, function(x) lapply(x, function(y) y$theta))
  parests <- lapply(fits, function(x) lapply(x, function(y) y$solution))
  oumvaests <- parests[[model]]
  oumvaests <- lapply(1:length(oumvaests), function(x) rbind(oumvaests[[x]], "Vy"=oumvaests[[x]][2,]/(2*oumvaests[[x]][1,])))
  oumvatheta <- thetaests[[model]]
  oumvaests2 <- do.call(cbind, lapply(oumvaests, melt) %>% lapply(function(x) x[,3]))
  #Vyests <- do.call(rbind, lapply(seq(1, 10, 2), function(x) sqrt(oumvaests2[x+1,]/(2*oumvaests2[x,1]))))
  #oumvaests2 <- rbind(oumvaests2, Vyests)
  medianparests <- apply(oumvaests2, 1, median)
  separests <- apply(oumvaests2, 1, sd)/(sqrt(ncol(oumvaests2)))
  oumvatheta2 <- sapply(oumvatheta, function(x) x[,1])
  medianthetaests <- apply(oumvatheta2, 1, median)
  sethetaests <- apply(oumvatheta2, 1, sd)/sqrt(ncol(oumvatheta2))
  parnames <- do.call(cbind, lapply(oumvaests, melt))[,1:2]
  sumOUMVA <- rbind(cbind(parnames, medianparests, separests), data.frame("Var1"="theta", "Var2"=parnames[seq(1, length(parnames[,2]), 3),2], "medianparests"=medianthetaests, "separests"=sethetaests))
  summaryOUMVA <- dcast(sumOUMVA, Var2 ~ Var1, value.var=c("medianparests"))
  summarySEOUMVA <- dcast(sumOUMVA, Var2 ~ Var1, value.var="separests")
  colnames(summarySEOUMVA) <- paste("SE.", colnames(summarySEOUMVA), sep="")
  SummaryOUMVA <- cbind(summaryOUMVA, summarySEOUMVA[,2:4])
  SummaryOUMVA
}

require(reshape2)
parests <- extractModel(5, fits$OU4)

par(mfrow=c(2,1))
plot(tdOUwie$dat$tmin.lo, tdOUwie$dat$lnVs, pch=21, bg=as.numeric(factor(tdOUwie$dat$OU4)))
library(denstrip)
densregion <- denstrip::densregion
theta <- parests$theta
Vy <- parests$Vy
cols <- 1:4
transparency <- 100
for(i in 1:4){
  ylim <- par("usr")[3:4]
  x <- seq(-50, 50, length = 10)
  y <- seq(ylim[1], ylim[2], length = 100)
  Z <- matrix(nrow = length(x), ncol = length(y))
  for (j in 1:length(x)) {
    Z[j, ] <- dnorm(y, theta[i], sqrt(Vy[i])/1)
  }
  if (sum(Z) != 0) {
    densregion(x, y, Z, colmax = makeTransparent(cols[i], 
                                               transparency), colmin = "transparent")
  }
}
abline(h=theta, lwd=2, lty=2, col=cols)
plot(theta, pch=21, cex=2, bg=cols, ylim=ylim, xlim=c(0.5,4.5))
lapply(1:4, function(x) lines(rep(x,2), c(theta[x]-2*sqrt(Vy[x]), theta[x]+2*sqrt(Vy[x])), col=x, lty=2, lwd=2))

dev.off()