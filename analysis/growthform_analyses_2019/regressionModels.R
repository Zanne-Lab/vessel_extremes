## Preliminary data preparation
td <- readRDS("../output/cleandata/matchednewdata3.rds")
td <- filter(td, !is.na(growth.form.cornwell.unique), growth.form.cornwell.unique %in% c("F", "C"))
td <- filter(td, !is.na(precip), !is.na(temperature))
#write.csv(td$dat, "../output/cleandata/noclimate.csv")

tree <- td$phy
dat <- td$dat
tree <- multi2di(tree, random=FALSE) # Resolve polytomies
tree$edge.length[tree$edge.length==0] <- .Machine$double.eps # Bayou doesn't like 0 length branches
pred <- dat
pred$growthC <- as.numeric(pred$growth.form.cornwell.unique=="C")
pred <- mutate(pred, lnVs =log(vesselSize), lnVn= log(vesselNumber), lnVsR = log(vesselSize/vesselNumber), lnLF=log(vesselSize*vesselNumber), tempK = (273.15+temperature), precip = precip, freeze=as.numeric(tmin.lo<0))
.pred <- dplyr::select(pred, lnVn, lnLF, growthC, tempK, precip, freeze)
.pred[,-c(2:3, 6)] <- as.data.frame(apply(.pred[,-c(2:3,6)], 2, scale))
#pred <- select(pred, lnVs, lnVn, lnVsR, lnLs, growthC, humidity, elevation,  tempK, tcoldq.me, precip,pdryq.me, alt.me)
# Int/lnVn/growthC/tempK/precip/growthCxlnVn
#g <- nrow(.pred)
#g*solve(t(as.matrix(cbind(1, .pred, .pred$growthC*.pred$tempK)))%*%as.matrix(cbind(1,.pred, .pred$growthC*.pred$tempK)))

cache <- bayou:::.prepare.ou.univariate(tree, setNames(log(dat$vesselSize), tree$tip.label), pred = .pred)
cacheLF <- bayou:::.prepare.ou.univariate(tree, setNames(.pred$lnLF, tree$tip.label), pred = .pred)
tmp <- lm(cache$dat ~ lnVn + growthC + growthC*tempK + precip + freeze, data=.pred)
tmp2 <- lm(cacheLF$dat ~ lnVn + growthC + growthC*tempK + precip + freeze, data=.pred)
#tmp2 <- lm(cache$dat ~ lnVn + growthC + precip, data=.pred)

summary(tmp)

## Growth form models
## Priors for different parameters, so that these are easily changed for all models
param.alpha <- list(scale=1)
param.sig2 <- list(scale=1)
param.beta_lnVn <- list(mean=0, sd=1)
param.beta_growthC <- list(mean=0, sd=1)
param.beta_tempK <- list(mean=0, sd=1)
param.beta_growthCxtempK <- list(mean=0, sd=1)
param.beta_precip <- list(mean=0, sd=1)
param.beta_freeze <- list(mean=0, sd=1)
param.k <- list(lambda=10, kmax=65.5)
param.sb <- list(bmax=1, prob=1)
param.theta <- list(mean=mean(cache$dat), sd=2.5*sd(cache$dat))
param.thetaLF <- list(mean=mean(cacheLF$dat), sd=2.5*sd(cacheLF$dat))


## Proposal widths
D.XX0000 <- function(nrj) list(alpha=0.75, sig2= 0.5, beta_lnVn=0.2, k=rep(1,nrj), theta=0.25, slide=1)
D.XXXXXX <- function(nrj) list(alpha=0.75, sig2= 0.5, beta_lnVn=0.2, beta_growthC=1, beta_tempK=0.2,
                               beta_growthCxtempK=0.2,beta_precip=0.2,k=rep(1,nrj), theta=0.25, slide=1)

## Models with LnVs
{
## Simplest regression; 110000
prior.110000 <- make.prior(tree, plot.prior = FALSE, 
                    dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_lnVn="dnorm",
                               #dbeta_growthC="dnorm", dbeta_tempK="dnorm", dbeta_growthCxtempK="dnorm", dbeta_precip="dnorm",
                               dsb="fixed", dk="fixed", dtheta="dnorm"), 
                    param=list(dalpha=param.alpha, dsig2=param.sig2, dbeta_lnVn=param.beta_lnVn,
                               dk="fixed", dsb="fixed", 
                               dtheta=param.theta),
                    fixed=list(k=0, ntheta=1, sb=numeric(0), loc=numeric(0))
)

model.110000 <- makeBayouModelDev(lnVs ~ lnVn, rjpars = c(), cache, prior.110000, D=D.XX0000(1))
prior.110000(model.110000$startpar)
model.110000$model$lik.fn(model.110000$startpar, cache, cache$dat)$loglik

## Rj intercept regression; N10000
prior.N10000 <- make.prior(tree, plot.prior = FALSE, 
                           dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_lnVn="dnorm",
                                      #dbeta_growthC="dnorm", dbeta_tempK="dnorm", dbeta_growthCxtempK="dnorm", dbeta_precip="dnorm",
                                      dsb="dsb", dk="cdpois", dtheta="dnorm"), 
                           param=list(dalpha=param.alpha, dsig2=param.sig2, dbeta_lnVn=param.beta_lnVn,
                                      #dbeta_growthC=list(mean=0, sd=1), dbeta_tempK=list(mean=0, sd=0.5), dbeta_growthCxtempK=list(mean=0, sd=0.002),dbeta_precip=list(mean=0, sd=0.005),
                                      dk=param.k, dsb=param.sb, 
                                      dtheta=param.theta)
)

model.N10000 <- makeBayouModelDev(lnVs ~ lnVn, rjpars = c("theta"), cache, prior.N10000, D=D.XX0000(1))
prior.N10000(model.N10000$startpar, cache)
model.N10000$model$lik.fn(model.N10000$startpar, cache, cache$dat)$loglik

## Rj intercept + slope regression; NN0000
prior.NN0000 <- make.prior(tree, plot.prior = FALSE, 
                           dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_lnVn="dnorm",
                                      #dbeta_growthC="dnorm", dbeta_tempK="dnorm", dbeta_growthCxtempK="dnorm", dbeta_precip="dnorm",
                                      dsb="dsb", dk="cdpois", dtheta="dnorm"), 
                           param=list(dalpha=param.alpha, dsig2=param.sig2, dbeta_lnVn=param.beta_lnVn,
                                      #dbeta_growthC=list(mean=0, sd=1), dbeta_tempK=list(mean=0, sd=0.5), dbeta_growthCxtempK=list(mean=0, sd=0.002),dbeta_precip=list(mean=0, sd=0.005),
                                      dk=param.k, dsb=param.sb, dtheta=param.theta)
)

model.NN0000 <- makeBayouModelDev(lnVs ~ lnVn, rjpars = c("theta", "lnVn"), cache, prior.NN0000, D=D.XX0000(2))
prior.NN0000(model.NN0000$startpar, cache)
model.NN0000$model$lik.fn(model.NN0000$startpar, cache, cache$dat)$loglik


## Full model simple; 111111
prior.111111 <- make.prior(tree, plot.prior = FALSE, 
                           dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_lnVn="dnorm",
                                      dbeta_growthC="dnorm", dbeta_tempK="dnorm", dbeta_growthCxtempK="dnorm", dbeta_precip="dnorm",
                                      dsb="fixed", dk="fixed", dtheta="dnorm"), 
                           param=list(dalpha=param.alpha, dsig2=param.sig2, dbeta_lnVn=param.beta_lnVn,
                                      dbeta_growthC=param.beta_growthC, dbeta_tempK=param.beta_tempK, dbeta_growthCxtempK=param.beta_growthCxtempK,dbeta_precip=param.beta_precip,
                                      dk="fixed", dsb="fixed", 
                                      dtheta=param.theta),
                            fixed=list(k=0, sb=numeric(0), loc=numeric(0))
)

model.111111 <- makeBayouModelDev(lnVs ~ lnVn + growthC + tempK + precip + tempK*growthC, rjpars = c(), cache, prior.111111, 
                               D = D.XXXXXX(1))
prior.111111(model.111111$startpar, cache)
model.111111$model$lik.fn(model.111111$startpar, cache, cache$dat)$loglik

## Rj intercept full model; N11111
prior.N11111 <- make.prior(tree, plot.prior = FALSE, 
                           dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_lnVn="dnorm",
                                      dbeta_growthC="dnorm", dbeta_tempK="dnorm", dbeta_growthCxtempK="dnorm", dbeta_precip="dnorm",
                                      dsb="dsb", dk="cdpois", dtheta="dnorm"), 
                           param=list(dalpha=param.alpha, dsig2=param.sig2, dbeta_lnVn=param.beta_lnVn,
                                      dbeta_growthC=param.beta_growthC, dbeta_tempK=param.beta_tempK, dbeta_growthCxtempK=param.beta_growthCxtempK,dbeta_precip=param.beta_precip,
                                      dk=param.k, dsb=param.sb, dtheta=param.theta)
                           )


model.N11111 <- makeBayouModelDev(lnVs ~ lnVn + growthC + tempK + precip + tempK*growthC, rjpars = c("theta"), cache, prior.N11111, 
                               D = D.XXXXXX(1))
prior.N11111(model.N11111$startpar, cache)
model.N11111$model$lik.fn(model.N11111$startpar, cache, cache$dat)$loglik


## Rj intercept + slope full model; NN1111
prior.NN1111 <- make.prior(tree, plot.prior = FALSE, 
                           dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_lnVn="dnorm",
                                      dbeta_growthC="dnorm", dbeta_tempK="dnorm", dbeta_growthCxtempK="dnorm", dbeta_precip="dnorm",
                                      dsb="dsb", dk="cdpois", dtheta="dnorm"), 
                           param=list(dalpha=param.alpha, dsig2=param.sig2, dbeta_lnVn=param.beta_lnVn,
                                      dbeta_growthC=param.beta_growthC, dbeta_tempK=param.beta_tempK, dbeta_growthCxtempK=param.beta_growthCxtempK,dbeta_precip=param.beta_precip,
                                      dk=param.k, dsb=param.sb, dtheta=param.theta)
)

model.NN1111 <- makeBayouModelDev(lnVs ~ lnVn + growthC + tempK + precip + tempK*growthC, rjpars = c("theta", "lnVn"), cache, prior.NN1111, 
                               D = D.XXXXXX(2))
prior.NN1111(model.NN1111$startpar, cache)
model.NN1111$model$lik.fn(model.NN1111$startpar, cache, cache$dat)$loglik

## 111110; simple full model - interaction
prior.111110 <- make.prior(tree, plot.prior = FALSE, 
                           dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_lnVn="dnorm",
                                      dbeta_growthC="dnorm", dbeta_tempK="dnorm", 
                                      #dbeta_growthCxtempK="dnorm", 
                                      dbeta_precip="dnorm",
                                      dsb="fixed", dk="fixed", dtheta="dnorm"), 
                           param=list(dalpha=param.alpha, dsig2=param.sig2, dbeta_lnVn=param.beta_lnVn,
                                      dbeta_growthC=param.beta_growthC, dbeta_tempK=param.beta_tempK, 
                                      #dbeta_growthCxtempK=param.beta_growthCxtempK,
                                      dbeta_precip=param.beta_precip,
                                      dk="fixed", dsb="fixed", dtheta=param.theta),
                           fixed=list(k=0, ntheta=1, sb=numeric(0), loc=numeric(0))
)

model.111110 <- makeBayouModelDev(lnVs ~ lnVn + growthC + tempK + precip, rjpars = c(), cache, prior.111110, 
                               D = D.XXXXXX(1))
prior.111110(model.111110$startpar, cache)
model.111110$model$lik.fn(model.111110$startpar, cache, cache$dat)$loglik

## 111001; simple model - climate - interaction
prior.111000 <- make.prior(tree, plot.prior = FALSE, 
                           dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_lnVn="dnorm",
                                      dbeta_growthC="dnorm", 
                                      #dbeta_tempK="dnorm", 
                                      #dbeta_growthCxtempK="dnorm", 
                                      #dbeta_precip="dnorm",
                                      dsb="fixed", dk="fixed", dtheta="dnorm"), 
                           param=list(dalpha=param.alpha, dsig2=param.sig2, dbeta_lnVn=param.beta_lnVn,
                                      dbeta_growthC=param.beta_growthC, 
                                      #dbeta_tempK=param.beta_tempK, 
                                      #dbeta_growthCxtempK=param.beta_growthCxtempK,
                                      #dbeta_precip=param.beta_precip,
                                      dk="fixed", dsb="fixed", dtheta=param.theta),
                           fixed=list(k=0, ntheta=1, sb=numeric(0), loc=numeric(0))
)

model.111000 <- makeBayouModelDev(lnVs ~ lnVn + growthC, rjpars = c(), cache, prior.111000, 
                               D = D.XXXXXX(1))
prior.111000(model.111000$startpar, cache)
model.111000$model$lik.fn(model.111000$startpar, cache, cache$dat)$loglik

## N11001; intercept shift model - climate - interaction
prior.N11000 <- make.prior(tree, plot.prior = FALSE, 
                           dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_lnVn="dnorm",
                                      dbeta_growthC="dnorm", 
                                      #dbeta_tempK="dnorm", 
                                      #dbeta_growthCxtempK="dnorm", 
                                      #dbeta_precip="dnorm",
                                      dsb="dsb", dk="cdpois", dtheta="dnorm"), 
                           param=list(dalpha=param.alpha, dsig2=param.sig2, dbeta_lnVn=param.beta_lnVn,
                                      dbeta_growthC=param.beta_growthC, 
                                      #dbeta_tempK=param.beta_tempK, 
                                      #dbeta_growthCxtempK=param.beta_growthCxtempK,
                                      #dbeta_precip=param.beta_precip,
                                      dk=param.k, dsb=param.sb, dtheta=param.theta)
)

model.N11000 <- makeBayouModelDev(lnVs ~ lnVn + growthC, rjpars = c("theta"), cache, prior.N11000, 
                               D = D.XXXXXX(1))
prior.N11000(model.N11000$startpar, cache)
model.N11000$model$lik.fn(model.N11000$startpar, cache, cache$dat)$loglik

## NN1000; intercept shift model - climate - interaction
prior.NN1000 <- make.prior(tree, plot.prior = FALSE, 
                           dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_lnVn="dnorm",
                                      dbeta_growthC="dnorm", 
                                      #dbeta_tempK="dnorm", 
                                      #dbeta_growthCxtempK="dnorm", 
                                      #dbeta_precip="dnorm",
                                      dsb="dsb", dk="cdpois", dtheta="dnorm"), 
                           param=list(dalpha=param.alpha, dsig2=param.sig2, dbeta_lnVn=param.beta_lnVn,
                                      dbeta_growthC=param.beta_growthC, 
                                      #dbeta_tempK=param.beta_tempK, 
                                      #dbeta_growthCxtempK=param.beta_growthCxtempK,
                                      #dbeta_precip=param.beta_precip,
                                      dk=param.k, dsb=param.sb, dtheta=param.theta)
)

model.NN1000 <- makeBayouModelDev(lnVs ~ lnVn + growthC, rjpars = c("theta", "lnVn"), cache, prior.NN1000, 
                               D = D.XXXXXX(2))
prior.NN1000(model.NN1000$startpar, cache)
model.NN1000$model$lik.fn(model.NN1000$startpar, cache, cache$dat)$loglik

## 111101; full simple model - precip
prior.111101 <- make.prior(tree, plot.prior = FALSE, 
                           dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_lnVn="dnorm",
                                      dbeta_growthC="dnorm", 
                                      dbeta_tempK="dnorm", 
                                      dbeta_growthCxtempK="dnorm", 
                                      #dbeta_precip="dnorm",
                                      dsb="fixed", dk="fixed", dtheta="dnorm"), 
                           param=list(dalpha=param.alpha, dsig2=param.sig2, dbeta_lnVn=param.beta_lnVn,
                                      dbeta_growthC=param.beta_growthC, 
                                      dbeta_tempK=param.beta_tempK, 
                                      dbeta_growthCxtempK=param.beta_growthCxtempK,
                                      #dbeta_precip=param.beta_precip,
                                      dk="fixed", dsb="fixed", dtheta=param.theta),
                           fixed=list(k=0, ntheta=1, sb=numeric(0), loc=numeric(0))
)

model.111101 <- makeBayouModelDev(lnVs ~ lnVn + growthC+tempK+growthC*tempK, rjpars = c(), cache, prior.111101, 
                               D = D.XXXXXX(1))
prior.111101(model.111101$startpar, cache)
model.111101$model$lik.fn(model.111101$startpar, cache, cache$dat)$loglik
}

## Models with LF
## Code: theta + growthC + temp + precip + growthC*tempK + freeze + freeze*growthC
{
## Proposal widths
D.LFXX00000 <- function(nrj) list(alpha=0.75, sig2= 0.5, beta_growthC=1, k=rep(1,nrj), theta=0.5, slide=1,
                                  beta_tempK=0.2, beta_growthCxtempK=0.2,beta_precip=0.2,k=rep(1,nrj))
D.LFNXXXXXX <- function(nrj) list(alpha=0.75, sig2= 0.5, beta_growthC=1, beta_tempK=0.2,
                                  beta_growthCxtempK=0.2,beta_precip=0.2,k=rep(1,nrj), theta=1, slide=1)

## Growth form only
prior.LF1100000 <- make.prior(tree, plot.prior = FALSE, 
                           dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", 
                                      #dbeta_lnVn="dnorm",
                                      dbeta_growthC="dnorm", 
                                      #dbeta_tempK="dnorm", dbeta_growthCxtempK="dnorm", dbeta_precip="dnorm",
                                      dsb="fixed", dk="fixed", dtheta="dnorm"), 
                           param=list(dalpha=param.alpha, dsig2=param.sig2, 
                                      #dbeta_lnVn=param.beta_lnVn,
                                      dbeta_growthC=param.beta_growthC,
                                      dk="fixed", dsb="fixed", 
                                      dtheta=param.theta),
                           fixed=list(k=0, ntheta=1, sb=numeric(0), loc=numeric(0))
)

model.LF1100000 <- makeBayouModelDev(lnLF ~ growthC, rjpars = c(), cacheLF, prior.LF1100000, D=D.LFXX00000(1))
prior.LF1100000(model.LF1100000$startpar)
model.LF1100000$model$lik.fn(model.LF1100000$startpar, cacheLF, cacheLF$dat)$loglik

## Growth form & temperature
prior.LF1110000 <- make.prior(tree, plot.prior = FALSE, 
                              dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", 
                                         #dbeta_lnVn="dnorm",
                                         dbeta_growthC="dnorm", 
                                         dbeta_tempK="dnorm", 
                                         #dbeta_growthCxtempK="dnorm", dbeta_precip="dnorm",
                                         dsb="fixed", dk="fixed", dtheta="dnorm"), 
                              param=list(dalpha=param.alpha, dsig2=param.sig2, 
                                         #dbeta_lnVn=param.beta_lnVn,
                                         dbeta_growthC=param.beta_growthC,
                                         dbeta_tempK = param.beta_tempK,
                                         dk="fixed", dsb="fixed", 
                                         dtheta=param.theta),
                              fixed=list(k=0, ntheta=1, sb=numeric(0), loc=numeric(0))
)

model.LF1110000 <- makeBayouModelDev(lnLF ~ growthC + tempK, rjpars = c(), cacheLF, prior.LF1110000, D=D.LFXX00000(1))
prior.LF1110000(model.LF1110000$startpar)
model.LF1110000$model$lik.fn(model.LF1110000$startpar, cacheLF, cacheLF$dat)$loglik

## Growth form & temperature & precip
prior.LF1111000 <- make.prior(tree, plot.prior = FALSE, 
                              dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", 
                                         #dbeta_lnVn="dnorm",
                                         dbeta_growthC="dnorm", 
                                         dbeta_tempK="dnorm", 
                                         #dbeta_growthCxtempK="dnorm", 
                                         dbeta_precip="dnorm",
                                         dsb="fixed", dk="fixed", dtheta="dnorm"), 
                              param=list(dalpha=param.alpha, dsig2=param.sig2, 
                                         #dbeta_lnVn=param.beta_lnVn,
                                         dbeta_growthC=param.beta_growthC,
                                         dbeta_tempK = param.beta_tempK,
                                         dbeta_precip = param.beta_precip,
                                         dk="fixed", dsb="fixed", 
                                         dtheta=param.theta),
                              fixed=list(k=0, ntheta=1, sb=numeric(0), loc=numeric(0))
)

model.LF1111000 <- makeBayouModelDev(lnLF ~ growthC + tempK + precip, rjpars = c(), cacheLF, prior.LF1111000, D=D.LFXX00000(1))
prior.LF1111000(model.LF1111000$startpar)
model.LF1111000$model$lik.fn(model.LF1111000$startpar, cacheLF, cacheLF$dat)$loglik

## Growth form & temperature & precip & growthC*tempK
prior.LF1111100 <- make.prior(tree, plot.prior = FALSE, 
                              dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", 
                                         #dbeta_lnVn="dnorm",
                                         dbeta_growthC="dnorm", 
                                         dbeta_tempK="dnorm", 
                                         dbeta_growthCxtempK="dnorm", 
                                         dbeta_precip="dnorm",
                                         dsb="fixed", dk="fixed", dtheta="dnorm"), 
                              param=list(dalpha=param.alpha, dsig2=param.sig2, 
                                         #dbeta_lnVn=param.beta_lnVn,
                                         dbeta_growthC=param.beta_growthC,
                                         dbeta_tempK = param.beta_tempK,
                                         dbeta_precip = param.beta_precip,
                                         dbeta_growthCxtempK = param.beta_growthCxtempK,
                                         dk="fixed", dsb="fixed", 
                                         dtheta=param.theta),
                              fixed=list(k=0, ntheta=1, sb=numeric(0), loc=numeric(0))
)

model.LF1111100 <- makeBayouModelDev(lnLF ~ growthC + tempK + precip +  growthC*tempK, rjpars = c(), cacheLF, prior.LF1111100, D=D.LFXX00000(1))
prior.LF1111100(model.LF1111100$startpar)
model.LF1111100$model$lik.fn(model.LF1111100$startpar, cacheLF, cacheLF$dat)$loglik

## 
## RJ & Growth form 
prior.LFN100000 <- make.prior(tree, plot.prior = FALSE, 
                              dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", 
                                         #dbeta_lnVn="dnorm",
                                         dbeta_growthC="dnorm", 
                                         #dbeta_tempK="dnorm", 
                                         #dbeta_growthCxtempK="dnorm", 
                                         #dbeta_precip="dnorm",
                                         dsb="dsb", dk="cdpois", dtheta="dnorm"), 
                              param=list(dalpha=param.alpha, dsig2=param.sig2, 
                                         #dbeta_lnVn=param.beta_lnVn,
                                         dbeta_growthC=param.beta_growthC,
                                         #dbeta_tempK = param.beta_tempK,
                                         #dbeta_precip = param.beta_precip,
                                         #dbeta_growthCxtempK = param.beta_growthCxtempK,
                                         dk=param.k, dsb=param.sb, 
                                         dtheta=param.theta))

model.LFN100000 <- makeBayouModelDev(lnLF ~ growthC, rjpars = c("theta"), cacheLF, prior.LFN100000, D=D.LFNXXXXXX(1))
prior.LFN100000(model.LFN100000$startpar)
model.LFN100000$model$lik.fn(model.LFN100000$startpar, cacheLF, cacheLF$dat)$loglik

## RJ & growthC + tempK
prior.LFN110000 <- make.prior(tree, plot.prior = FALSE, 
                              dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", 
                                         #dbeta_lnVn="dnorm",
                                         dbeta_growthC="dnorm", 
                                         dbeta_tempK="dnorm", 
                                         #dbeta_growthCxtempK="dnorm", 
                                         #dbeta_precip="dnorm",
                                         dsb="dsb", dk="cdpois", dtheta="dnorm"), 
                              param=list(dalpha=param.alpha, dsig2=param.sig2, 
                                         #dbeta_lnVn=param.beta_lnVn,
                                         dbeta_growthC=param.beta_growthC,
                                         dbeta_tempK = param.beta_tempK,
                                         #dbeta_precip = param.beta_precip,
                                         #dbeta_growthCxtempK = param.beta_growthCxtempK,
                                         dk=param.k, dsb=param.sb, 
                                         dtheta=param.theta))

model.LFN110000 <- makeBayouModelDev(lnLF ~ growthC + tempK, rjpars = c("theta"), cacheLF, prior.LFN110000, D=D.LFNXXXXXX(1))
prior.LFN110000(model.LFN110000$startpar)
model.LFN110000$model$lik.fn(model.LFN110000$startpar, cacheLF, cacheLF$dat)$loglik


## RJ & growthC + tempK + precip
prior.LFN111000 <- make.prior(tree, plot.prior = FALSE, 
                              dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", 
                                         #dbeta_lnVn="dnorm",
                                         dbeta_growthC="dnorm", 
                                         dbeta_tempK="dnorm", 
                                         #dbeta_growthCxtempK="dnorm", 
                                         dbeta_precip="dnorm",
                                         dsb="dsb", dk="cdpois", dtheta="dnorm"), 
                              param=list(dalpha=param.alpha, dsig2=param.sig2, 
                                         #dbeta_lnVn=param.beta_lnVn,
                                         dbeta_growthC=param.beta_growthC,
                                         dbeta_tempK = param.beta_tempK,
                                         dbeta_precip = param.beta_precip,
                                         #dbeta_growthCxtempK = param.beta_growthCxtempK,
                                         dk=param.k, dsb=param.sb, 
                                         dtheta=param.theta))

model.LFN111000 <- makeBayouModelDev(lnLF ~ growthC + tempK + precip, rjpars = c("theta"), cacheLF, prior.LFN111000, D=D.LFNXXXXXX(1))
prior.LFN111000(model.LFN111000$startpar)
model.LFN111000$model$lik.fn(model.LFN111000$startpar, cacheLF, cacheLF$dat)$loglik

## RJ & growthC + tempK + precip + growthC*tempK
prior.LFN111100 <- make.prior(tree, plot.prior = FALSE, 
                              dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", 
                                         #dbeta_lnVn="dnorm",
                                         dbeta_growthC="dnorm", 
                                         dbeta_tempK="dnorm", 
                                         dbeta_growthCxtempK="dnorm", 
                                         dbeta_precip="dnorm",
                                         dsb="dsb", dk="cdpois", dtheta="dnorm"), 
                              param=list(dalpha=param.alpha, dsig2=param.sig2, 
                                         #dbeta_lnVn=param.beta_lnVn,
                                         dbeta_growthC=param.beta_growthC,
                                         dbeta_tempK = param.beta_tempK,
                                         dbeta_precip = param.beta_precip,
                                         dbeta_growthCxtempK = param.beta_growthCxtempK,
                                         dk=param.k, dsb=param.sb, 
                                         dtheta=param.theta))

model.LFN111100 <- makeBayouModelDev(lnLF ~ growthC + tempK + precip + growthC*tempK, rjpars = c("theta"), cacheLF, prior.LFN111100, D=D.LFNXXXXXX(1))
prior.LFN111100(model.LFN111100$startpar)
model.LFN111100$model$lik.fn(model.LFN111100$startpar, cacheLF, cacheLF$dat)$loglik
}
priors <- lapply(objects()[grep("prior.", objects(), fixed=TRUE)], function(x) get(x)); names(priors) <- objects()[grep("prior.", objects(), fixed=TRUE)]
models <- lapply(objects()[grep("model.", objects(), fixed=TRUE)], function(x) get(x)$model); names(models) <- objects()[grep("model.", objects(), fixed=TRUE)]
startpars <- lapply(objects()[grep("model.", objects(), fixed=TRUE)], function(x) get(x)$startpar); names(startpars) <- objects()[grep("model.", objects(), fixed=TRUE)]


