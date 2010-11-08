mpmm.generate <- function(T,params,dims) {
  D <- length(dims)
  K <- length(params$phi)
  dataset <- matrix(0,T,D)
  z <- rep(0,T)
  for (i in 1:T) {
    z[i] <- discrete.draw(params$pi)
    for (d in 1:D) {
      zd <- z[i]
      dataset[i,d] <- discrete.draw(params$phi[[d]][zd,])
    }
  }
  list(dataset=dataset,assignments=z)
}


