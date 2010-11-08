
rownormalize <- function(x) {
  if (is.null(dim(x)))
    return(x/sum(x))
  rs <- rowSums(x)
  tmp <- do.call(cbind,rep(list(rs),ncol(x)))
  x / tmp
}

checkDataset <- function(dataset) {
  if (class(dataset)!="matrix")
    warning("Dataset should be a matrix.")
  if (! ncol(dataset) %in% c(2,3))
    warning("Dataset should be a matrix with either two or three columns.")
  if (any(is.na(dataset)))
    warning("Dataset contains NAs.  Missing data is not yet supported.")
}

mpmm.llk <- function(dataset,params,dims) {
  return(sum(log(mpmm.predict(dataset,params,dims))))
}
mpmm.predictive.dist <- function(params,dims) {
  D <- length(dims)
  vecs <- lapply(dims, function(x) {1:x})
  edgelist <- expand.grid(vecs)
  prob <- mpmm.predict(edgelist,params,dims)
  colnames(edgelist) <- paste("d",1:D,sep="")
  return(cbind(edgelist,prob=prob))
}

discrete.draw <- function(probs)
  1 + sum(runif(1) > cumsum(probs)) # or  which(rmultinom(1,1,probs) == 1)

