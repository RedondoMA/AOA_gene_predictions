phyEstimateDisc_mine<-function (phy, trait, best.state = TRUE, cutoff = 0.5, ...) 
{
  if (is.vector(trait) | is.factor(trait)) {
    trait <- data.frame(trait)
  }
  trait[, 1] <- factor(trait[, 1])
  trait.orig <- trait
  sppObs <- row.names(trait)
  sppUnobs <- phy$tip.label[!(phy$tip.label %in% sppObs)]
  trtlevels <- levels(trait[, 1])
  res <- as.data.frame(matrix(nrow = length(sppUnobs), ncol = length(trtlevels), 
                              dimnames = list(sppUnobs, trtlevels)))
  for (i in sppUnobs) {
    tree <- drop.tip(phy, subset(sppUnobs, sppUnobs != i))
    tree <- root(tree, i, resolve.root = FALSE)
    edge <- Nnode(tree) - 1 + which(tree$tip.label == i)
    bl <- tree$edge.length[edge]
    tree <- drop.tip(tree, i)
    trait <- trait.orig[tree$tip.label, ]
    est <- ace(trait, tree, type = "discrete", ...)
    val <- est$lik.anc[1, ]
    res[i, ] <- val
  }
  if (best.state) {
    beststate <- as.data.frame(matrix(nrow = dim(res)[1], 
                                      ncol = 2))
    colnames(beststate) <- c("estimated.state", "estimated.state.support")
    rownames(beststate) <- rownames(res)
    for (i in 1:dim(res)[1]) {
      res2<-as.numeric(res[i,])
      names(res2)<-colnames(res)
      best <- -sort(-(res2))[1]
      if (best >= cutoff) {
        beststate[i, 1] <- names(best)
        beststate[i, 2] <- best
      }
      else {
        beststate[i, 1] <- NA
        beststate[i, 2] <- NA
      }
    }
  }
  if (best.state) {
    return(cbind(as.matrix(res), beststate))
  }
  else {
    return(as.matrix(res))
  }
}
