
mpmm.plot <- function(fit) {
  pred <- mpmm.predictive.dist(fit$params,fit$dims)
  names(pred)[ncol(pred)] <- "value"
  plotmat(pred)
}

plotmat <- function(mat,labels=c("Sender","Receiver"),limits=c(0,max(mat[,3]))) {
  require(ggplot2)
  colnames(mat)[3] <- "value"
  qq <- ggplot(mat,aes(x=d2,y=d1)) +
       geom_tile(aes(fill=value)) + 
       scale_fill_gradient(low="white",high="black",limits=limits)+
       theme_bw() + labs(x=labels[2],y=labels[1],fill="Probability")+
       coord_equal(ratio=1) +
       opts(legend.position = "none",
            panel.grid.minor=theme_blank(),
            panel.grid.major=theme_blank())
  qq
}

mpmm.plotphi <- function(fit) {
  require(ggplot2)
  res <- c()
  for (i in 1:length(fit$phi)) {
    df <- melt(fit$phi[[i]])
    names(df) <- c("class","entity","value")
    df$dim <- i
    res <- rbind(res,df)
  }
  qq <- ggplot(res,aes(x=entity,y=factor(class))) + facet_grid(dim~.)
  qq + geom_tile(aes(fill=value)) + scale_fill_gradient(low="white",high="black") +
       theme_bw() + labs(x="Actor",y="Class",fill="Prob.") +
       opts(panel.grid.minor=theme_blank(),panel.grid.major=theme_blank()) 
}

mult.dir <- function(dataset,dims,prior=1/prod(dims)) {
  checkDataset(dataset)
  if (ncol(dataset)==2)
    e <- paste(dataset[,1],dataset[,2])
  if (ncol(dataset)==3)
    e <- paste(dataset[,1],dataset[,2],dataset[,3])
  tb <- table(e)
  prob <- (tb + prior)/(sum(tb + prior) + prod(dims)*prior)
  b <- strsplit(names(tb)," ")
  b <- lapply(b,as.numeric)
  b <- do.call(rbind,b)
  b <- cbind(b,tb,prob)
  rownames(b) <- c()
  colnames(b) <- c(paste("d",1:length(dims),sep=""),"count","prob")
  b <- as.data.frame(b)
  return(b)
}

mpmm.plotassignments <- function(dataset,fit,sort=TRUE) {
  require(ggplot2)  
  # Only handles 2-dimensional data
  d <- cbind(dataset,z=fit$assignments)
  b <- mult.dir(d,c(max(dataset[,1]),max(dataset[,2]),fit$params$C))

  if (sort) {
    mapSenderGroups <- apply(fit$params$phi[[1]],2,which.max)
    so <- order(mapSenderGroups)
    mapReceiverGroups <- apply(fit$params$phi[[2]],2,which.max)
    to <- order(mapReceiverGroups)
    b[,1] <- reorder(factor(b[,1]),match(b[,1],so))
    b[,2] <- reorder(factor(b[,2]),match(b[,2],to))
  }

  b$d3 <- factor(b$d3)
  qq <- ggplot(b,aes(x=d2,y=d1)) +
        geom_point(aes(alpha=log(count),size=log(count),colour=d3))+
        scale_size(to=c(1,4))+
        facet_grid(~d3) + labs(x="receiver",y="sender",colour="class") +
        opts(axis.text.x=theme_text(angle=-90, hjust=0, size=5),
             axis.text.y=theme_text(vjust=0, size=5))
  qq
}
