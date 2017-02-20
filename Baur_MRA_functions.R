# MRA functions following Baur and Leuenberger 2011
# code kindly provided by H. Baur
# keep this script in the working directory, it will be sourced during the 
### analysis of skull measurements


##########################################
##### PCA in isometry free shape space
##########################################

##### Arguments
 # x           a numeric matrix (or data frame) that provides the data for the principal component analysis
 # npc         number of shape components shown in numerical output, as usually only the first two or three components contain useful information. Default value is npc=3
 # rpc         number of decimal places used for printing components and loadings. Default value is rpc=4 

##### Values
 # shapePCA    returns the results of the shape PCA as an object of class "prcomp". Mainly used for displaying the screeplot
 # loadings    a matrix of variable loadings (i.e., a matrix whose columns contain the eigenvectors)
 # components  a matrix containing the variance, proportion of varaiance and cumulative proportion of variance
 # PCmatrix    a matrix of rotated data (i.e., the centered data multiplied by the loadings)
 # isosize     a matrix containing isometric size (i.e., the centered data multiplide by the isometric size vector)

shapePCA <- function(x, npc=3, rpc=4){
X <- log(x)
X <- scale(X, center=TRUE, scale=FALSE)
p <- dim(X)[2]
I <- diag(1,p,p)
a0 <- as.vector(rep(1,p)/p)
P <- I-p*(a0%*%t(a0))
Y <- X%*%P
colnames(Y) <- colnames(X)
isosize <- X%*%a0; colnames(isosize) <- "isosize"
PCA <- prcomp(Y,center=FALSE, scale=FALSE)
loadings <- PCA$rotation
p <- dim(loadings)[2]
PCmatrix <- PCA$x
components.a <- PCA$sdev^2; components.b <- components.a/sum(components.a); components.c <- cumsum(components.b); components <- rbind(components.a, components.b, components.c)
colnames(loadings) <- paste("shape.PC", 1:p, sep="")
colnames(PCmatrix) <- colnames(loadings)
colnames(components) <- colnames(loadings); rownames(components) <- c("Variance", "Proportion of Variance", "Cumulative Proportion")
list(PCA=PCA, PCmatrix=PCmatrix[,1:npc], isosize=as.numeric(isosize), loadings=round(loadings[,1:npc], rpc), components=round(components[,1:npc], rpc))
}





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




