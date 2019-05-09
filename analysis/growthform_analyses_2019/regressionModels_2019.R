## OLS analyses 
#ols <- list()
#ols[["NCTx"]] <- lm(dat ~ logN + growthC + growthC*tempK, data=.pred)
#ols[["NCTPx"]] <- lm(dat ~ logN + growthC + growthC*tempK + precip, data=.pred)
#ols[["NCTPxx"]] <- lm(dat ~ logN + growthC*tempK + growthC*precip, data=.pred)
#ols[["NCTP"]] <- lm(dat ~ logN + growthC + tempK + precip, data=.pred)
#ols[["NCT"]] <- lm(dat ~ logN + growthC + tempK, data=.pred)
#ols[["NCP"]] <- lm(dat ~ logN + growthC + precip, data=.pred)
#ols[["NC"]] <- lm(dat ~ logN + growthC, data=.pred)
#ols[["NTPx"]] <- lm(dat ~ logN + precip*tempK, data=.pred)
#sort(sapply(ols, AIC))

## Growth form models
## Priors for different parameters, so that these are easily changed for all models
param.alpha <- list(scale=0.05)
param.sig2 <- list(scale=0.1)
param.beta_logN <- list(mean=-1, sd=0.25)
param.beta_growthC <- list(mean=0, sd=1)
param.beta_tempK <- list(mean=0, sd=1)
param.beta_growthCxtempK <- list(mean=0, sd=1)
param.beta_precip <- list(mean=0, sd=1)
#param.beta_freeze <- list(mean=0, sd=1)
param.k <- list(lambda=1, kmax=65.5)
param.sb <- list(bmax=1, prob=1)
param.sb2 <- list(bmax=Inf, prob=tree$edge.length)
param.theta <- list(mean=mean(dat), sd=2.5*sd(dat))
#param.thetaLF <- list(mean=mean(cacheLF$dat), sd=2.5*sd(cacheLF$dat))

dist.names <- list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_logN="dnorm",dbeta_growthC="dnorm", 
                   dbeta_tempK="dnorm", dbeta_precip="dnorm", dbeta_growthCxtempK="dnorm", 
                   dk="fixed", dsb="fixed",  dtheta="dnorm")

dist.names.rj <- list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_logN="dnorm",dbeta_growthC="dnorm", 
                   dbeta_tempK="dnorm", dbeta_precip="dnorm", dbeta_growthCxtempK="dnorm", 
                   dk="cdpois", dsb="dsb",  dtheta="dnorm")

dist.params.fixed <- list(dalpha=param.alpha, dsig2=param.sig2, dbeta_logN=param.beta_logN, dbeta_growthC=param.beta_growthC, dbeta_tempK=param.beta_tempK,  
                      dbeta_precip=param.beta_precip,  dbeta_growthCxtempK=param.beta_growthCxtempK, dk="fixed", dsb="fixed", dtheta=param.theta)

dist.params.rj <- list(dalpha=param.alpha, dsig2=param.sig2, dbeta_logN=param.beta_logN, dbeta_growthC=param.beta_growthC, dbeta_tempK=param.beta_tempK,  
                          dbeta_precip=param.beta_precip,  dbeta_growthCxtempK=param.beta_growthCxtempK, dk=param.k, dsb=param.sb, dtheta=param.theta)

## Proposal widths
D.XX000 <- function(nrj) list(alpha=0.75, sig2= 0.5, beta_logN=0.2, k=rep(1,nrj), theta=0.25, slide=1)
D.XXXXX <- function(nrj) list(alpha=0.75, sig2= 0.5, beta_logN=0.25, beta_growthC=1, beta_tempK=0.25,
                               beta_growthCxtempK=0.25, beta_precip=0.25,k=rep(1,nrj), theta=0.25, slide=1)

priors <- list()
models <- list()
## Models with LogS
{
## Simplest regression; 
  
priors[['m11000']] <- make.prior(tree, plot.prior = FALSE, 
                    dists=dist.names[c(1:3, 8:10)], 
                    param=dist.params.fixed[c(1:3, 8:10)],
                    fixed=list(k=0, ntheta=1, sb=numeric(0), loc=numeric(0))
)

models[["m11000"]] <- makeBayouModel(dat ~ logN, rjpars = c(), tree, dat, pred=.pred, priors$m11000, D=D.XX000(1))
priors[["m11000"]](models$m11000$startpar)
#model.11000$model$lik.fn(model.11000$startpar, cache, cache$dat)$loglik


## GrowthC; 11100
priors[['m11100']] <- make.prior(tree, plot.prior = FALSE, 
                          dists=dist.names[c(1:3, 4, 8:10)], 
                          param=dist.params.fixed[c(1:3, 4, 8:10)],
                          fixed=list(k=0, ntheta=1, sb=numeric(0), loc=numeric(0))
)

models[['m11100']] <- makeBayouModel(dat ~ logN + growthC, rjpars = c(), tree, dat, pred=.pred, priors$m11100, D=D.XXXXX(1))
priors$m11100(models$m11100$startpar)
#model.11000$model$lik.fn(model.11000$startpar, cache, cache$dat)$loglik

## GrowthC & tempK; 11110
priors[['m11110']] <- make.prior(tree, plot.prior = FALSE, 
                                 dists=dist.names[c(1:3, 4:5, 8:10)], 
                                 param=dist.params.fixed[c(1:3, 4:5, 8:10)],
                                 fixed=list(k=0, ntheta=1, sb=numeric(0), loc=numeric(0))
)

models[['m11110']] <- makeBayouModel(dat ~ logN + growthC+tempK, rjpars = c(), tree, dat, pred=.pred, priors$m11110, D=D.XXXXX(1))
priors$m11110(models$m11110$startpar)
#model.11000$model$lik.fn(model.11000$startpar, cache, cache$dat)$loglik


## GrowthC & tempK & precip; 11111
priors[['m11111']] <- make.prior(tree, plot.prior = FALSE, 
                                 dists=dist.names[c(1:3, 4:6, 8:10)], 
                                 param=dist.params.fixed[c(1:3, 4:6, 8:10)],
                                 fixed=list(k=0, ntheta=1, sb=numeric(0), loc=numeric(0))
)

models[['m11111']] <- makeBayouModel(dat ~ logN + growthC + tempK + precip, rjpars = c(), tree, dat, pred=.pred, priors$m11111, D=D.XXXXX(1))
priors$m11111(models$m11111$startpar)
#model.11000$model$lik.fn(model.11000$startpar, cache, cache$dat)$loglik



## GrowthC x tempK & precip; 11111x
priors[['m11111x']] <- make.prior(tree, plot.prior = FALSE, 
                                 dists=dist.names[c(1:3, 4:7, 8:10)], 
                                 param=dist.params.fixed[c(1:3, 4:7, 8:10)],
                                 fixed=list(k=0, ntheta=1, sb=numeric(0), loc=numeric(0))
)

models[['m11111x']] <- makeBayouModel(dat ~ logN + growthC*tempK + precip, rjpars = c(), tree, dat, pred=.pred, priors$m11111x, D=D.XXXXX(1))
priors$m11111x(models$m11111x$startpar)
#model.11000$model$lik.fn(model.11000$startpar, cache, cache$dat)$loglik

#####Reversible Jump Models########
#########Intercept only############

priors[['mN1000']] <- make.prior(tree, plot.prior = FALSE, 
                                 dists=dist.names.rj[c(1:3, 8:10)], 
                                 param=dist.params.rj[c(1:3, 8:10)])


models[["mN1000"]] <- makeBayouModel(dat ~ logN, rjpars = c("theta"), tree, dat, pred=.pred, priors$mN1000, D=D.XXXXX(1))
priors[["mN1000"]](models$mN1000$startpar)
#model.11000$model$lik.fn(model.11000$startpar, cache, cache$dat)$loglik


## GrowthC; N1100
priors[['mN1100']] <- make.prior(tree, plot.prior = FALSE, 
                                 dists=dist.names.rj[c(1:3, 4, 8:10)], 
                                 param=dist.params.rj[c(1:3, 4, 8:10)])

models[['mN1100']] <- makeBayouModel(dat ~ logN + growthC, rjpars = c("theta"), tree, dat, pred=.pred, priors$mN1100, D=D.XXXXX(1))
priors$mN1100(models$mN1100$startpar)
#model.11000$model$lik.fn(model.11000$startpar, cache, cache$dat)$loglik

## GrowthC & tempK; N1110
priors[['mN1110']] <- make.prior(tree, plot.prior = FALSE, 
                                 dists=dist.names.rj[c(1:3, 4:5, 8:10)], 
                                 param=dist.params.rj[c(1:3, 4:5, 8:10)]
)

models[['mN1110']] <- makeBayouModel(dat ~ logN + growthC+tempK, rjpars = c("theta"), tree, dat, pred=.pred, priors$mN1110, D=D.XXXXX(1))
priors$mN1110(models$mN1110$startpar)
#model.11000$model$lik.fn(model.11000$startpar, cache, cache$dat)$loglik


## GrowthC & tempK & precip; N1111
priors[['mN1111']] <- make.prior(tree, plot.prior = FALSE, 
                                 dists=dist.names.rj[c(1:3, 4:6, 8:10)], 
                                 param=dist.params.rj[c(1:3, 4:6, 8:10)])

models[['mN1111']] <- makeBayouModel(dat ~ logN + growthC + tempK + precip, rjpars = c("theta"), tree, dat, pred=.pred, priors$mN1111, D=D.XXXXX(1))
priors$mN1111(models$mN1111$startpar)
#model.11000$model$lik.fn(model.11000$startpar, cache, cache$dat)$loglik



## GrowthC x tempK & precip; N1111x
priors[['mN1111x']] <- make.prior(tree, plot.prior = FALSE, 
                                  dists=dist.names.rj[c(1:3, 4:7, 8:10)], 
                                  param=dist.params.rj[c(1:3, 4:7, 8:10)]
)

models[['mN1111x']] <- makeBayouModel(dat ~ logN + growthC*tempK + precip, rjpars = c("theta"), tree, dat, pred=.pred, priors$mN1111x, D=D.XXXXX(1))
priors$mN1111x(models$mN1111x$startpar)


#########Intercept & Slope############

priors[['mNN000']] <- make.prior(tree, plot.prior = FALSE, 
                                 dists=dist.names.rj[c(1:3, 8:10)], 
                                 param=dist.params.rj[c(1:3, 8:10)])


models[["mNN000"]] <- makeBayouModel(dat ~ logN, rjpars = c("theta", "beta_logN"), tree, dat, pred=.pred, priors$mNN000, D=D.XXXXX(2))
priors[["mNN000"]](models$mNN000$startpar)
#model.11000$model$lik.fn(model.11000$startpar, cache, cache$dat)$loglik


## GrowthC; NN100
priors[['mNN100']] <- make.prior(tree, plot.prior = FALSE, 
                                 dists=dist.names.rj[c(1:3, 4, 8:10)], 
                                 param=dist.params.rj[c(1:3, 4, 8:10)])

models[['mNN100']] <- makeBayouModel(dat ~ logN + growthC, rjpars = c("theta", "logN"), tree, dat, pred=.pred, priors$mNN100, D=D.XXXXX(2))
priors$mNN100(models$mNN100$startpar)
#model.11000$model$lik.fn(model.11000$startpar, cache, cache$dat)$loglik

## GrowthC & tempK; N1110
priors[['mNN110']] <- make.prior(tree, plot.prior = FALSE, 
                                 dists=dist.names.rj[c(1:3, 4:5, 8:10)], 
                                 param=dist.params.rj[c(1:3, 4:5, 8:10)]
)

models[['mNN110']] <- makeBayouModel(dat ~ logN + growthC+tempK, rjpars = c("theta", "logN"), tree, dat, pred=.pred, priors$mNN110, D=D.XXXXX(2))
priors$mNN110(models$mNN110$startpar)
#model.11000$model$lik.fn(model.11000$startpar, cache, cache$dat)$loglik


## GrowthC & tempK & precip; N1111
priors[['mNN111']] <- make.prior(tree, plot.prior = FALSE, 
                                 dists=dist.names.rj[c(1:3, 4:6, 8:10)], 
                                 param=dist.params.rj[c(1:3, 4:6, 8:10)])

models[['mNN111']] <- makeBayouModel(dat ~ logN + growthC + tempK + precip, rjpars = c("theta", "logN"), tree, dat, pred=.pred, priors$mNN111, D=D.XXXXX(2))
priors$mNN111(models$mNN111$startpar)
#model.11000$model$lik.fn(model.11000$startpar, cache, cache$dat)$loglik



## GrowthC x tempK & precip; N1111x
priors[['mNN111x']] <- make.prior(tree, plot.prior = FALSE, 
                                  dists=dist.names.rj[c(1:3, 4:7, 8:10)], 
                                  param=dist.params.rj[c(1:3, 4:7, 8:10)]
)

models[['mNN111x']] <- makeBayouModel(dat ~ logN + growthC*tempK + precip, rjpars = c("theta", "logN"), tree, dat, pred=.pred, priors$mNN111x, D=D.XXXXX(2))
priors$mNN111x(models$mNN111x$startpar)

}