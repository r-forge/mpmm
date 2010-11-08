mpmm.predict <- function(dataset,params,dims) {
  if (ncol(dataset) < 3) {
    dataset <- cbind(dataset,1)
    params$phi[[3]] <- matrix(1,params$C,1)
  }
  
  T <- nrow(dataset)
  probs <- matrix(0,T,C)
  for (k in 1:params$C) {
    probs[,k] <- params$phi[[1]][k,dataset[,1]] *
                 params$phi[[2]][k,dataset[,2]] *
                 params$phi[[3]][k,dataset[,3]] *
                 params$pi[k]
  }
  probs <- rowSums(probs)
  return(probs)
}
