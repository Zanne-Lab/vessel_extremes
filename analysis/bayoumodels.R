growthformModels <- 
  list(
    "VsNR ~ rj(theta) + rj(tcoldq.me) + rj(pdryq.me) + growthC" =
      list(
        lik = function(pars, cache, X, model="Custom"){
          n <- cache$n
          X <- cache$dat
          pred <- cache$pred
          betaID <- getTipMap(pars, cache)
          ## Specify the model here
          X = X -pars$beta1[betaID]*pred$tcoldq.me-pars$beta2[betaID]*pred$pdryq.me- pars$betaC*pred$growthC
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
      )
    "VsNR ~ rj(theta) + temp + lshumid + elevation + growthC + growthC*temperature = 

    )