require(bayou)
require(devtools)
source("./analysis/betaBayouFunctions.R")
load_all("~/repos/bayou/bayou_1.0")
args <- list("rjgrowthforms", 10000, "r1")

td <- readRDS("./datasets/matcheddata.rds")
tree <- td$phy
dat <- td$data
tree <- multi2di(tree, random=FALSE) # Resolve polytomies
tree$edge.length[tree$edge.length==0] <- .Machine$double.eps # Bayou doesn't like 0 length branches
pred <- dat
pred$growthT <- as.numeric(pred$growthForm=="T")
pred$growthS_T <- as.numeric(pred$growthForm=="S_T")
pred$growthS <- as.numeric(pred$growthForm=="S")
pred$growthH <- as.numeric(pred$growthForm=="H")
pred$growthC <- as.numeric(pred$growthForm=="C")

cache <- bayou:::.prepare.ou.univariate(tree, setNames(log(dat$vesselSize), tree$tip.label), pred = pred)

custom.lik <- function(pars, cache, X, model="Custom"){
  n <- cache$n
  X <- cache$dat
  pred <- cache$pred
  betaID <- getTipMap(pars, cache)
  ## Specify the model here
  X = X - pars$beta1[betaID]*pred$log.no - pars$betaC*pred$growthC - pars$betaH*pred$growthH - pars$betaT*pred$growthT
  cache$dat <- X
  ### The part below mostly does not change
  X.c <- bayou:::C_weightmatrix(cache, pars)$resid
  transf.phy <- bayou:::C_transf_branch_lengths(cache, 1, X.c, pars$alpha)
  transf.phy$edge.length[cache$externalEdge] <- transf.phy$edge[cache$externalEdge] + cache$SE[cache$phy$edge[cache$externalEdge, 2]]^2*(2*pars$alpha)/pars$sig2
  comp <- bayou:::C_threepoint(list(n=n, N=cache$N, anc=cache$phy$edge[, 1], des=cache$phy$edge[, 2], diagMatrix=transf.phy$diagMatrix, P=X.c, root=transf.phy$root.edge, len=transf.phy$edge.length))
  if(pars$alpha==0){
    inv.yVy <- comp$PP
    detV <- comp$logd
  } else {
    inv.yVy <- comp$PP*(2*pars$alpha)/(pars$sig2)
    detV <- comp$logd+n*log(pars$sig2/(2*pars$alpha))
  }
  llh <- -0.5*(n*log(2*pi)+detV+inv.yVy)
  return(list(loglik=llh, theta=pars$theta,resid=X.c, comp=comp, transf.phy=transf.phy))
}

## Define starting parameters 
k <- rpois(1, lambda=10)
sb <- sample(1:length(cache$bdesc), k, replace=FALSE, prob = sapply(cache$bdesc, length))
startpar <- list(alpha=0.1, sig2=3, beta1=rnorm(k+1, -1, 1),betaC=0, betaH=0, betaT=0,k=k, ntheta=k+1, theta=rnorm(k+1, -1.25, 0.25), sb=sb, loc=rep(0, k), t2=2:(k+1))

monitor = function(i, lik, pr, pars, accept, accept.type, j){
  names <- c("gen", "lnL", "prior", "alpha","sig2", "rbeta1","betaC", "betaH","betaT","rtheta", "k")
  string <- "%-8i%-8.2f%-8.2f%-8.2f%-8.2f%-8.2f%-8.2f%-8.2f%-8.2f%-8.2f%-8i"
  acceptratios <- tapply(accept, accept.type, mean)
  names <- c(names, names(acceptratios))
  if(j==0){
    cat(sprintf("%-7.7s", names), "\n", sep=" ")                           
  }
  cat(sprintf(string, i, lik, pr, pars$alpha, pars$sig2, pars$beta1[1], pars$theta[1], pars$betaC, pars$betaH, pars$betaT, pars$k), sprintf("%-8.2f", acceptratios),"\n", sep="")
}

## We're going to define a custom model with variable slopes (beta1)
model.WA <- list(moves = list(alpha=".multiplierProposal", sig2=".multiplierProposal", 
                                   beta1=".vectorMultiplier", betaC=".slidingWindowProposal", betaH=".slidingWindowProposal",
                              betaT=".slidingWindowProposal",
                              theta=".adjustTheta", k=".splitmergebd", slide=".slide"
          ),
control.weights = list(alpha=4, sig2=2, beta1=10, betaC=4, betaH=4, betaT=4,
                       k=10, theta=10, slide=2 
          ),
D = list(alpha=1, sig2= 0.75, beta1=0.75,betaC=0.4, betaH=0.4,betaT=0.4, k=c(1,1), theta=2, slide=1
          ),
parorder = names(startpar)[-which(names(startpar) %in% c("sb", "loc", "t2"))],
rjpars = c("beta1","theta"),
shiftpars = c("sb", "loc", "t2"),
monitor.fn = monitor,
lik.fn = custom.lik)

prior <- make.prior(tree, plot.prior = FALSE, 
                    dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta1="dnorm",
                               dbetaC="dnorm", dbetaH="dnorm", dbetaT="dnorm",
                               dsb="dsb", dk="cdpois", dtheta="dnorm"
                    ), 
                    param=list(dalpha=list(scale=1), dsig2=list(scale=1), dbeta1=list(mean=-2, sd=1),
                               dbetaC=list(mean=0, sd=1),dbetaH=list(mean=0, sd=1),dbetaT=list(mean=0, sd=1),
                               dk=list(lambda=15, kmax=50), dsb=list(prob=1, bmax=1), 
                               dtheta=list(mean=-3, sd=1)
                    )
)

prior(startpar)
custom.lik(startpar, cache, cache$dat)$loglik


mymcmc <- bayou.makeMCMC(cache$phy, cache$dat, pred=cache$pred, SE=0, model=model.WA, prior=prior, startpar=startpar, new.dir=paste("./output/runs/",args[[1]],sep=""), outname=paste(args[[1]],args[[3]],sep=""), plot.freq=NULL, ticker.freq=1000, samp = 100)
mymcmc$run(50000)
chain <- mymcmc$load()
chain <- set.burnin(chain, 0.4)
plot(chain)
sumstats <- summary(chain)


cutoff <- 0.1
sumpars <- list(sb = which(sumstats$branch.posteriors$pp > cutoff))
sumpars$k <- length(sumpars$sb)
sumpars$ntheta <- length(sumpars$sb)+1
sumpars$loc <- rep(0, sumpars$k)
sumpars$t2 <- 2:sumpars$ntheta
sumpars$betaC <- sumstats$statistics['betaC', 1]
sumpars$betaH <- median(chain$betaH[200:500])
sumpars$betaT <- median(chain$betaT[200:500])
sumpars$theta <- sumstats$statistics['root.theta',1]
sumpars$beta1 <- sumstats$statistics['root.beta1',1]
tr <- pars2simmap(sumpars, tree)

plotSimmap(tr$tree, colors=tr$col, fsize=0.25)
summarizeDerivedState <- function(branch, chain){
  if(branch==0){
    Th <- sapply(chain$theta, function(x) x[1])
    B1 <- sapply(chain$beta1, function(x) x[1])
    B2 <- sapply(chain$beta2, function(x) x[1])
    #B3 <- unlist(chain$beta3)
  } else {
    SB <- unlist(chain$sb)
    gen <- unlist(lapply(1:length(chain$sb), function(x) rep(x, length(chain$sb[[x]]))))
    ind <- which(SB==branch)
    gen <- gen[ind]
    T2 <- unlist(chain$t2)[ind]
    B1 <- sapply(1:length(T2), function(x) chain$beta1[[gen[x]]][T2[x]])
    Th <- sapply(1:length(T2), function(x) chain$theta[[gen[x]]][T2[x]])
    B2 <- unlist(chain$beta2)[gen]
    #B3 <- unlist(chain$beta3)[gen]
  }
  medians = list(theta=median(Th), beta1=median(B1), beta2=median(B2))
  densities = list(theta=density(Th), beta1=density(B1), beta2=density(B2))
  return(list(medians=medians, densities=densities))
}


tipregs <- bayou:::.tipregime(sumpars, tree)
descendents <- lapply(1:(length(sb)+1), function(x) names(tipregs[tipregs==x])) 

par(mfrow=c(1,2))
plot(cache$pred$log.no, cache$dat, pch=21, bg=cache$pred$growthForm, col=cache$pred$growthForm, main="Regressions", ylab="log Size")
curve(sumpars$theta+sumpars$betaC+sumpars$beta1*x, add=TRUE, lty=2, lwd=2, col=1)
curve(sumpars$theta+sumpars$beta1*x, add=TRUE, lty=2, lwd=2, col=3)
curve(sumpars$theta+sumpars$betaH+sumpars$beta1*x, add=TRUE, lty=2, lwd=2, col=2)
curve(sumpars$theta+sumpars$betaT+sumpars$beta1*x, add=TRUE, lty=2, lwd=2, col="blue")


plot(density(sapply(chain$theta[200:500], function(x) x[1])), col="green", lwd=3, xlim=c(-6, 0), main="Densities", xlab="Intercept")
lines(density(chain$betaH[200:500]+sapply(chain$theta[200:500], function(x) x[1])), col="red", lwd=3)
lines(density(chain$betaC[200:500]+sapply(chain$theta[200:500], function(x) x[1])), col="black", lwd=3)
lines(density(chain$betaT[200:500]+sapply(chain$theta[200:500], function(x) x[1])), col="blue", lwd=3)

