mpmm.cgs <- function(dataset,C,dims,priors,niter=10,z.init=NULL) {
  checkDataset(dataset)
    
  # C code only handles edge data with 3 vertices; add a column if dataset has fewer.
  if (ncol(dataset) < 3) {
    dataset <- cbind(dataset,1)
    dims <- c(dims,1)
    priors$phi[[3]] <- 1
  }

  s <- dataset[,1]
  r <- dataset[,2]
  a <- dataset[,3]
  alpha <- priors$pi
  beta <- priors$phi[[1]]
  gamma <- priors$phi[[2]]
  delta <- priors$phi[[3]]

  zh <- z.init  # initial state assignments for all events
  if (is.null(zh))
    zh <- sample(1:C,nrow(dataset),replace=TRUE)

  # Compute sufficient statistics
  mode_counts <- table(factor(zh,1:C))
  sender_counts <- table(factor(zh,1:C), factor(s,1:dims[1]))
  receiver_counts <- table(factor(zh,1:C), factor(r,1:dims[2]))
  action_counts <- table(factor(zh,1:C), factor(a,1:dims[3]))

  # Gibbs sample class assignments for each event
  storage.mode(s) <- storage.mode(r) <- storage.mode(a) <- "integer"
  storage.mode(niter) <- storage.mode(zh) <- "integer"
  b <- .Call ("mpmm_cgs",s,r,a,mode_counts,sender_counts,receiver_counts,action_counts,
              alpha,beta,gamma,delta,niter,zh,PACKAGE="mpmm")

  # Smoothed estimates
  spi <- as.vector((b$mode_counts+alpha)/sum(b$mode_counts+alpha))
  sphi <- list()
  sphi[[1]] <- rownormalize(b$sender_counts + beta)
  sphi[[2]] <- rownormalize(b$receiver_counts + gamma)
  counts <- list(b$sender_counts, b$receiver_counts)

  if (dims[3]>1) {
    sphi[[3]] <- rownormalize(b$action_counts + delta)
    counts <- c(counts,b$action_counts)
  } else {
    dims <- dims[1:2]
  }
  
  results <- list(assignments = b$assignments,dims = dims,
                  params=list(pi=spi,phi=sphi,C=C),
                  priors=priors, counts=counts)
                  
  cat("\n")
  return(results)
}
