require(bayou)
require(devtools)
require(dplyr)
require(aRbor)
source("./betaBayouFunctions.R")
args <- commandArgs(TRUE)
args <- as.list(args)
args[[2]] <- as.numeric(args[[2]])
load_all("~/repos/bayou/bayou_1.0")
#args <- list("rjsimpleregression", 10000, "r3")
print(args)

td <- readRDS("../output/cleandata/matchednewdata.rds")
td <- treeply(td, reorder, "postorder")
tree <- td$phy
dat <- td$dat
tree <- multi2di(tree, random=FALSE) # Resolve polytomies
tree$edge.length[tree$edge.length==0] <- .Machine$double.eps # Bayou doesn't like 0 length branches
pred <- dat
pred$growthC <- as.numeric(pred$support=="C")
pred <- mutate(pred, lnVs =log(vesselSize), lnVn= log(vesselNumber), lnVsR = log(vesselSize/vesselNumber), lnLs=log(vesselSize*vesselNumber), tempK = 273.15+temperature)
pred <- select(pred, lnVs, lnVn, lnVsR, lnLs, growthC, humidity, elevation,  tempK, tcoldq.me, precip,pdryq.me, alt.me)
trait <- pred$lnVs

tmp <- cbind(td$phy$tip.label, select(pred, lnVs, lnVn, growthC))
tmp2 <- tmp[order(apply(tmp[,2:3], 1, sum)),]
subset(tmp2, tmp2[,4]==1)

cache <- bayou:::.prepare.ou.univariate(tree, setNames(trait, tree$tip.label), pred = pred)

custom.lik <- function(pars, cache, X, model="Custom"){
  n <- cache$n
  X <- cache$dat
  pred <- cache$pred
  betaID <- getTipMap(pars, cache)
  ## Specify the model here
  X = X - pars$beta1[betaID]*pred$lnVn
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
k <- rpois(1, lambda=25)
sb <- sample(1:length(cache$bdesc), k, replace=FALSE, prob = sapply(cache$bdesc, length))
startpar <- list(alpha=0.1, sig2=3, beta1=rnorm(k+1, -1, 1), k=k, ntheta=k+1, theta=rnorm(k+1, -1.25, 0.25), sb=sb, loc=rep(0, k), t2=2:(k+1))
plotBayoupars(startpar, cache$phy, col=setNames(rainbow(startpar$ntheta), 1:startpar$ntheta))

monitor = function(i, lik, pr, pars, accept, accept.type, j){
  names <- c("gen", "lnL", "prior", "alpha","sig2", "rbeta1", "rtheta", "k")
  string <- "%-8i%-8.2f%-8.2f%-8.2f%-8.2f%-8.2f%-8.2f%-8i"
  acceptratios <- tapply(accept, accept.type, mean)
  names <- c(names, names(acceptratios))
  if(j==0){
    cat(sprintf("%-7.7s", names), "\n", sep=" ")                           
  }
  cat(sprintf(string, i, lik, pr, pars$alpha, pars$sig2, pars$beta1[1], pars$theta[1], pars$k), sprintf("%-8.2f", acceptratios),"\n", sep="")
}

## We're going to define a custom model with variable slopes (beta1)
model.WA <- list(moves = list(alpha=".multiplierProposal", sig2=".multiplierProposal", 
                                   beta1=".vectorMultiplier", theta=".adjustTheta", k=".splitmergebd", slide=".slide"
          ),
control.weights = list(alpha=4, sig2=2, beta1=10, 
                       k=10, theta=10, slide=2 
          ),
D = list(alpha=1, sig2= 0.75, beta1=0.75, k=c(1,1), theta=2, slide=1
          ),
parorder = names(startpar)[-which(names(startpar) %in% c("sb", "loc", "t2"))],
rjpars = c("beta1","theta"),
shiftpars = c("sb", "loc", "t2"),
monitor.fn = monitor,
lik.fn = custom.lik)

prior <- make.prior(tree, plot.prior = FALSE, 
                    dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta1="dnorm",
                               dsb="dsb", dk="cdpois", dtheta="dnorm"
                    ), 
                    param=list(dalpha=list(scale=1), dsig2=list(scale=1), dbeta1=list(mean=-2, sd=1), 
                               dk=list(lambda=15, kmax=50), dsb=list(prob=1, bmax=1), 
                               dtheta=list(mean=-3, sd=1)
                    )
)

prior(startpar)
custom.lik(startpar, cache, cache$dat)$loglik


mymcmc <- bayou.makeMCMC(cache$phy, cache$dat, pred=cache$pred, SE=0, model=model.WA, prior=prior, startpar=startpar, new.dir=paste("../output/runs/",args[[1]],sep=""), outname=paste(args[[1]],args[[3]],sep=""), plot.freq=NULL, ticker.freq=1000, samp = 100)
mymcmc$run(200000)
chain <- mymcmc$load()
chain <- set.burnin(chain, 0.4)
plot(chain)

save(chain, file=paste("../output/runs/",args[[1]],"/",args[[1]],"_chain_",args[[3]],sep=""))
save(mymcmc, file=paste("../output/runs/",args[[1]],"/",args[[1]],"_mcmc_",args[[3]],sep=""))

load(file=paste("../output/runs/",args[[1]],"/",args[[1]],"_chain_",args[[3]],".Rdata",sep=""))
load(file=paste("../output/runs/",args[[1]],"/",args[[1]],"_chain_",args[[3]], ".Rdata",sep=""))

sumstats <- summary(chain)
cutoff <- 0.2
postburn <- seq(floor(0.3*length(chain$gen)), length(chain$gen))
sumpars <- list(sb = which(sumstats$branch.posteriors$pp > cutoff))
sumpars$k <- length(sumpars$sb)
sumpars$ntheta <- length(sumpars$sb)+1
sumpars$loc <- rep(0, sumpars$k)
sumpars$t2 <- 2:sumpars$ntheta


summarizeDerivedState <- function(branch, chain){
  if(branch==0){
    Th <- sapply(chain$theta, function(x) x[1])
    B1 <- sapply(chain$beta1, function(x) x[1])
    #B2 <- sapply(chain$beta2, function(x) x[1])
    #B3 <- unlist(chain$beta3)
  } else {
    SB <- unlist(chain$sb)
    gen <- unlist(lapply(1:length(chain$sb), function(x) rep(x, length(chain$sb[[x]]))))
    ind <- which(SB==branch)
    gen <- gen[ind]
    T2 <- unlist(chain$t2)[ind]
    B1 <- sapply(1:length(T2), function(x) chain$beta1[[gen[x]]][T2[x]])
    Th <- sapply(1:length(T2), function(x) chain$theta[[gen[x]]][T2[x]])
    #B2 <- unlist(chain$beta2)[gen]
    #B3 <- unlist(chain$beta3)[gen]
  }
  medians = list(theta=median(Th), beta1=median(B1))
  densities = list(theta=density(Th), beta1=density(B1))
  return(list(medians=medians, densities=densities))
}

sb <- sumpars$sb
cladesummaries <- lapply(c(0, sb), function(x) summarizeDerivedState(x, chain))
regressions <- t(sapply(cladesummaries, function(x) unlist(x$medians)))
rownames(regressions) = c("root", sumpars$sb[sb])
#cache <- bayou:::.prepare.ou.univariate(, dat, SE=0, pred)
tipregs <- bayou:::.tipregime(sumpars, tree)
descendents <- lapply(1:(length(sb)+1), function(x) names(tipregs[tipregs==x])) 
#descendents <- lapply(sb, function(x) na.exclude(cache$tip.label[cache$edge[c(cache$bdesc[[x]], x), 2]] ))
nodesc <- which(sapply(descendents, length)==0)


pal <- heat.colors
tr <- pars2simmap(sumpars, cache$phy)
pdf("../output/figures/rjsimpleregressionmap.pdf", height=20, width=12)
plotSimmap(tr$tree, colors = tr$col, fsize=0.15)

par(mfrow=c(2,2), mar=c(3, 3, 1, 1), bg="white")
plot(c(-7, 0), c(0, 4), type="n", main="Theta", ylab="density", xlab="Theta")
lapply(1:length(cladesummaries), function(x) lines(cladesummaries[[x]]$densities$theta, col=x))
plot(c(-2.5, 0), c(0, 15), type="n", main="Beta1", ylab="density", xlab="Beta1")
lapply(1:length(cladesummaries), function(x) lines(cladesummaries[[x]]$densities$beta1, col=x))

par(mfrow=c(2,2), mar=c(3,3,5,1), bg="black", ask=FALSE, col.axis="white", col.lab="white", col.main="white")
c1 <- "pink"
palx <- rainbow
#beta3 <- median(chains$beta3)
#missing.pred <- do.call(rbind, chainsmp$missing.pred)
#mps <- apply(missing.pred, 2, median)
#impPred <- pred
#impPred[is.na(pred[,3]),3] <- mps
for(i in (1:nrow(regressions))){
  plotBayoupars(sumpars, tree, colors=setNames(c(palx(nrow(regressions))[i], rep("gray80", nrow(regressions)-1)), c(i, (1:nrow(regressions))[-i])), fsize=0.2)
  plot(pred$lnVn, dat, xlab="lnMass", ylab="lnBMR", pch=21, bg=makeTransparent("gray20", 100), col =makeTransparent("gray20", 100) )
  include <- which(names(dat) %in% descendents[[i]])
  text(pred[include, 'lnVn'], dat[include], labels=names(dat[include]), col="white", cex=0.4, pos = 4)
  points(pred[include,'lnVn'], dat[include], pch=21, bg=palx(nrow(regressions))[i])
  print(descendents[[i]])
  expected <- regressions[i,1]+regressions[i,2]*pred[include,'lnVn']
  o <- order(pred[include,1])
  lines(pred[include,'lnVn'][o], expected[o], col=palx(nrow(regressions))[i], lwd=2)
  plot(cladesummaries[[i]]$densities$beta1, col=palx(nrow(regressions))[i], xlim=c(-2.5, 0), lwd=3, main="Beta2")
  plot(cladesummaries[[i]]$densities$theta, col=palx(nrow(regressions))[i], xlim=c(-7,0), lwd=3, main="Theta")
  
}

dev.off()
#par(mfrow=c(1,2), bg="black")




#plotBranchHeatMap(cache$phy, chain, "theta", burnin=0.3, nn=1000, pal, show.tip.label=FALSE, legend_ticks=seq(0.65, 0.85, 0.05), edge.width=2)
#plotBranchHeatMap(cache$phy, chain, "beta1", burnin=0.3, nn=1000, pal, show.tip.label=FALSE, legend_ticks=seq(-2, 0, 0.5), edge.width=2)
#dev.off()














#pdf("../output/figures/lnVnlnvsLabeled.pdf")
#plot(pred$lnVn, pred$lnVs,  pch=21,col=makeTransparent("gray50", 100), bg=makeTransparent("gray50", 100), xlab="ln Vessel Number", ylab="ln Vessel Size")
#points(pred$lnVn[pred$growthC==1],pred$lnVs[pred$growthC==1], pch=21, col=makeTransparent("red",100), bg=makeTransparent("red",100))
#text(pred$lnVn[which(pred$growthC==1)],pred$lnVs[which(pred$growthC==1)], label=cache$phy$tip.label[which(pred$growthC==1)],pos=4, cex=0.25, col="darkblue")
#dev.off()
