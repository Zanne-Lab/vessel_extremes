source("workspace.R")
transition.calc <- function(tree, continuous, discrete, simmap.model="ER", simmap.nsim=1000, simmap.pi="estimated", anc.ML.maxit=100000){
  #Assertions and argument checking
  if(!inherits(tree, "phylo")) stop("Error: '", deparse(substitute(simmap)), "' must be of class 'phylo'")
  if(!is.factor(discrete)) stop("Error: '", deparse(substitute(discrete)), "' must be a factor; preferably a discrete character!")
  if(is.null(names(discrete))) stop("Error: '", deparse(substitute(discrete)), "' must be named")
  if(!is.numeric(continuous)) stop("Error: '", deparse(substitute(continuous)), "' must be a numeric; preferably a continuous character!")
  if(is.null(names(continuous))) stop("Error: '", deparse(substitute(continuous)), "' must be named")
  if(!identical(sort(tree$tip.label), sort(names(discrete)))) stop("Error: mismatch between'", deparse(substitute(discrete)), "' and phylogeny")
  if(!identical(sort(tree$tip.label), sort(names(continuous)))) stop("Error: mismatch between'", deparse(substitute(continuous)), "' and phylogeny")
  
  #Make simmap
  simmap <- make.simmap(tree, discrete, model=simmap.model, nsim=simmap.nsim, pi=simmap.pi)
  t.ltt <- ltt(simmap[[1]], plot=FALSE, gamma=FALSE)
  
  #Find transitions (can be multiple per branch); pre-allocation could make too long a data.frame (CHECK!)
  states <- levels(discrete)
  tmp <- rep(NA, length(simmap)*length(simmap[[1]]$maps))
  transitions <- data.frame(from=tmp, to=tmp, end.node=tmp)
  x <- 1
  for(i in seq_along(simmap)){
      for(j in seq_along(simmap[[i]]$maps)){
          if(tree$edge[j,2]>length(tree$tip.label) & length(simmap[[i]]$maps[[j]]) > 1){
              for(k in seq(from=1, to=length(simmap[[i]]$maps[[j]])-1)){
                  transitions$from[x] <- states[which(states==names(simmap[[i]]$maps[[j]])[k])]
                  transitions$to[x] <- states[which(states==names(simmap[[i]]$maps[[j]])[k+1])]
                  transitions$end.node[x] <- tree$edge[j,2]
                  x <- x + 1
                  if(x >= nrow(transitions))
                      transitions <- rbind(transitions, data.frame(from=tmp, to=tmp, end.node=tmp))
              }
          } else {
              transitions$from[x] <- states[which(states==names(simmap[[i]]$maps[[j]])[1])]
              transitions$to[x] <- states[which(states==names(simmap[[i]]$maps[[j]])[1])]
              transitions$end.node[x] <- tree$edge[j,2]
              x <- x + 1
              if(x == nrow(transitions))
                  transitions <- rbind(transitions, data.frame(from=tmp, to=tmp, end.node=tmp))
          }
      }
  }
  transitions <- transitions[!is.na(transitions$from),]
  #Age the transitions
  transitions$age <- max(t.ltt$times) - t.ltt$times[match(transitions$end.node, names(t.ltt$times))]
  transitions$first <- c(TRUE, rep(FALSE, nrow(transitions)-1))
  
  #Reconstruct continuous state
  anc.continuous <- fastAnc(tree, continuous, CI=TRUE)
  #Get the modal reconstructed nodal value and plot against that
  transitions$transition <- with(transitions, paste(from, to, sep="_"))
  #Prepare output and return
  output <- list(transitions=transitions, cont.sim=anc.continuous)
  class(output) <- "transition.calc"
  return(output)
}




#Silly models!
c.data$data$evergreen <- ifelse(c.data$data$phenology == "EV", 1, 0)
c.data$data$interac <- with(c.data$data, log(vesselSize) * pole.lim)
model <- with(c.data, transition.calc(phy, setNames(data$interac, rownames(data)), setNames(factor(data$evergreen), rownames(data)), simmap.nsim=1000))

counts <- with(model$transitions, table(transition, end.node))
modal.trans <- setNames(rownames(counts)[unlist(apply(counts, 2, function(x) which(max(x) == x)[1]))], unique(model$transitions$end.node))
modal.trans.frac <- setNames(numeric(length(modal.trans)), names(modal.trans))
for(i in seq_along(modal.trans.frac))
    modal.trans.frac[i] <- counts[which(rownames(counts)==modal.trans[i]),i] / sum(counts[,i])
modal.trans <- modal.trans[names(modal.trans) %in% names(model$cont.sim$ace)]
t <- model$cont.sim$ace[names(model$cont.sim$ace) %in% names(modal.trans)]
identical(names(modal.trans), names(t))
save.image("trans_wip.RData")
