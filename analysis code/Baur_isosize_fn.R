# MRA functions following Baur and Leuenberger 2011
# code kindly provided by H. Baur
# keep this script in the working directory, it will be sourced during the 
### analysis of skull measurements

##### Arguments
# x           a numeric matrix (or data frame) that provides the data for the principal component analysis
# npc         number of shape components shown in numerical output, as usually only the first two or three components contain useful information. Default value is npc=3
# rpc         number of decimal places used for printing components and loadings. Default value is rpc=4 

##### Values
# isosize    isometric size (i.e., the centered data multiplide by the isometric size vector)


##########################################
##### Isometric size
##########################################

isosize <- function(x) {
  X <- log(x)
  X <- scale(X, center=TRUE, scale=FALSE)
  p <- dim(X)[2]
  I <- diag(1,p,p)
  a0 <- as.vector(rep(1,p)/p)
  P <- I-p*(a0%*%t(a0))
  Y <- X%*%P
  colnames(Y) <- colnames(X)
  isosize <- X%*%a0
  isosize	
}




