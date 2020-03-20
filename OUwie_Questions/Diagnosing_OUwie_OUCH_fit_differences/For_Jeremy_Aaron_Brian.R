## The purpose of this code is to compare how OUwie and OUCH and calculate the covariance and weight matrices.

require(OUwie)
require(phytools)
require(corpcor)
require(ouch)
require(ape)
require(magrittr)

## Load trees, regimes, phenotypic data
treelist <- readRDS("treelist.RDS")
ou.tree <- treelist[["ou.tree"]] ## phylogeny in ouchtree format
ape.tree <- treelist[["ape.tree"]] ## phylogeny in phylo format
ou.reg <- treelist[["ou.reg"]] ## OU regime specification for ouch
ou.dat <- treelist[["ou.dat"]] ## phenotypic data in ouch format
ape.dat <- treelist[["ape.dat"]] ## phenotypic data in OUwie format

## Run both programs
ou.OU <- hansen(ou.dat, ou.tree, ou.reg, sqrt.alpha=1, sigma=1)
ouwie.OUM <- OUwie(ape.tree, ape.dat, model='OUM', root.station=TRUE, scaleHeight=TRUE)

## The two clearly converge on really different estimates
summary(ou.OU)
ouwie.OUM

## And they don't agree on the calculation of likelihood for the same set of parameters
hansen(ou.dat, ou.tree, ou.reg, sqrt.alpha=sqrt(ouwie.OUM$solution["alpha",1]), sigma=sqrt(ouwie.OUM$solution["sigma.sq",1]), fit=FALSE)@loglik
ouwie.OUM$loglik

#########################################################################
#                                                                       #
#	COMPUTING THE WEIGHT AND COVARIANCE MATRICES FOR OUCH           #
#                                                                       #
#########################################################################

## Create a single regime painting as well to compare between OUCH and OUwie to see if that is the source of differences
ou.reg.1 = as.factor(rep("x",length(ou.reg)))
names(ou.reg.1) <- as.character(1:length(ou.reg.1))
## rename the original two-regime painting to keep track
ou.reg.2 <- ou.reg

## Parameter estimates
sqrt.alpha <- 1
sigma <- 1

## From inside hansen.R hansen() function
source("hansen.R")
nm <- rownames(ou.dat)
data <- lapply(as.list(ou.dat),function(x){names(x)<-nm;x})
nchar <- length(ou.dat)
data <- lapply(ou.dat,function(x)x[ou.tree@nodes])
dat <- do.call(c,lapply(ou.dat,function(y)y[ou.tree@term]))
nsymargs <- nchar*(nchar+1)/2
nalpha <- length(sqrt.alpha)
nsigma <- length(sigma)

## Compute beta for the single regime case
nm.1 <- deparse(substitute(ou.reg.1))[1]
regimes.1 <- list(ou.reg.1)
names(regimes.1) <- nm.1
regimes.1 <- rep(regimes.1,nchar)
regimes.1 <- lapply(regimes.1,function(x)x[ou.tree@nodes])
beta.1 <- regime.spec(ou.tree,regimes.1)
## Compute beta for the two-regime case
nm.2 <- deparse(substitute(ou.reg.2))[1]
regimes.2 <- list(ou.reg.2)
names(regimes.2) <- nm.2
regimes.2 <- rep(regimes.2,nchar)
regimes.2 <- lapply(regimes.2,function(x)x[ou.tree@nodes])
beta.2 <- regime.spec(ou.tree,regimes.2)

## From inside hansen.R ou.lik.fn() function
alpha <- sym.par(sqrt.alpha)
sigma <- sym.par(sigma)
n <- length(ou.dat)
ev <- eigen(alpha,symmetric=TRUE)

## Compute the weight matrix for the single regime and two-regime cases
system("R CMD SHLIB weight-matrix.c")
dyn.load("weight-matrix.so")
ouch.W.1 <- .Call("ouch_weights",object=ou.tree,lambda=ev$values,S=ev$vectors,beta=beta.1)
ouch.W.2 <- .Call("ouch_weights",object=ou.tree,lambda=ev$values,S=ev$vectors,beta=beta.2)

## Compute the covariance matrix, which does not depend on regime paintings
system("R CMD SHLIB covar-matrix.c")
dyn.load("covar-matrix.so")
ouch.V <- .Call("ouch_covar",object=ou.tree,lambda=ev$values,S=ev$vectors,sigma.sq=sigma)
## This executes eqn (A8) in Butler and King 2004 for single regime and two-regime cases to compute the theta values
source("glssoln.R")
gsol.1 <- try(
    glssoln(ouch.W.1,dat,ouch.V),
    silent=FALSE
)
gsol.2 <- try(
    glssoln(ouch.W.2,dat,ouch.V),
    silent=FALSE
)
ouch.theta.1 <- gsol.1$coeff
ouch.theta.2 <- gsol.2$coeff

#########################################################################
#                                                                       #
#	COMPUTING THE WEIGHT AND COVARIANCE MATRICES FOR OUwie          #
#                                                                       #
#########################################################################

## Here are the key functions
source("varcov.ou.R")
source("weight.mat.R")
## the two-regime phylogeny and phenotypic data
phy.2 <- ape.tree
data.2 <- ape.dat
## single regime phylogeny and phenotypic data
phy.1 <- phy.2
phy.1$node.label <- rep("1",length(phy.1$node.label))
data.1 <- data.2
data.1$xMpleth[data.1$xMpleth=="2"] <- "1"

## Two-regime case first
data <- data.2
phy <- phy.2
data<-data.frame(data[,2], data[,3], row.names=data[,1])
data<-data[phy$tip.label,]
n=max(phy$edge[,1])
ntips=length(phy$tip.label)
tot.states<-factor(c(phy$node.label,as.character(data[,1])))
k<-length(levels(tot.states))
int.states<-factor(phy$node.label)
phy$node.label=as.numeric(int.states)
tip.states<-factor(data[,1])
data[,1]<-as.numeric(tip.states)

root.state <- phy$node.label[1]
int.state <- phy$node.label[-1]
##New tree matrix to be used for subsetting regimes
edges=cbind(c(1:(n-1)),phy$edge,nodeHeights(phy))
edges[,4:5]<-edges[,4:5]/max(nodeHeights(phy))
edges=edges[sort.list(edges[,3]),]
mm<-c(data[,1],int.state)
regime <- matrix(0,nrow=length(mm),ncol=length(unique(mm)))
##Generates an indicator matrix from the regime vector
for (i in 1:length(mm)) {
    regime[i,mm[i]] <- 1
}
##Finishes the edges matrix
edges=cbind(edges,regime)
##Resort the edge matrix so that it looks like the original matrix order
edges=edges[sort.list(edges[,1]),]
x<-as.matrix(data[,2])
index.mat<-matrix(0,2,k)
np=2
index<-matrix(TRUE,2,k)
index.mat[1,1:k]<-1
index.mat[2,1:k]<-2
param.count<-np+k
bool=TRUE

## parameter estimates
p <- c(1,1)
Rate.mat <- matrix(1, 2, k)
Rate.mat[] <- c(p, 1e-10)[index.mat]
N<-length(x[,1])
## covariance matrix does not depend on regime painting
ouwie.V<-varcov.ou(phy, edges, Rate.mat, root.state=1, simmap.tree=FALSE, scaleHeight=TRUE)
## but the weight matrix does
ouwie.W.2<-weight.mat(phy, edges, Rate.mat, root.state=1, simmap.tree=FALSE, scaleHeight=TRUE, assume.station=TRUE)
## here is the calculation of the thetas
ouwie.theta.2<-pseudoinverse(t(ouwie.W.2)%*%pseudoinverse(ouwie.V)%*%ouwie.W.2)%*%t(ouwie.W.2)%*%pseudoinverse(ouwie.V)%*%x

## Repeat for the single regime case
data <- data.1
phy <- phy.1
data<-data.frame(data[,2], data[,3], row.names=data[,1])
data<-data[phy$tip.label,]
n=max(phy$edge[,1])
ntips=length(phy$tip.label)
tot.states<-factor(c(phy$node.label,as.character(data[,1])))
k<-length(levels(tot.states))
int.states<-factor(phy$node.label)
phy$node.label=as.numeric(int.states)
tip.states<-factor(data[,1])
data[,1]<-as.numeric(tip.states)

root.state <- phy$node.label[1]
int.state <- phy$node.label[-1]
##New tree matrix to be used for subsetting regimes
edges=cbind(c(1:(n-1)),phy$edge,nodeHeights(phy))
edges[,4:5]<-edges[,4:5]/max(nodeHeights(phy))
edges=edges[sort.list(edges[,3]),]
mm<-c(data[,1],int.state)
regime <- matrix(0,nrow=length(mm),ncol=length(unique(mm)))
##Generates an indicator matrix from the regime vector
for (i in 1:length(mm)) {
    regime[i,mm[i]] <- 1
}
##Finishes the edges matrix
edges=cbind(edges,regime)
##Resort the edge matrix so that it looks like the original matrix order
edges=edges[sort.list(edges[,1]),]
x<-as.matrix(data[,2])
index.mat<-matrix(0,2,k)
np=2
index<-matrix(TRUE,2,k)
index.mat[1,1:k]<-1
index.mat[2,1:k]<-2
param.count<-np+k
bool=TRUE

## parameter estimates
p <- c(1,1)
Rate.mat <- matrix(1, 2, k)
Rate.mat[] <- c(p, 1e-10)[index.mat]
N<-length(x[,1])
ouwie.W.1<-weight.mat(phy, edges, Rate.mat, root.state=1, simmap.tree=FALSE, scaleHeight=TRUE, assume.station=TRUE)
## here is the calculation of the thetas
ouwie.theta.1<-pseudoinverse(t(ouwie.W.1)%*%pseudoinverse(ouwie.V)%*%ouwie.W.1)%*%t(ouwie.W.1)%*%pseudoinverse(ouwie.V)%*%x

##############################################################3
##
##               COMPARING THE RESULTS
##
##############################################################3

## If you plug the OUCH weight and covariance matrices for the single-regime case into the OUwie calcuation, you get results that are identical to the OUwie results. (Note that I have to reverse the order of the data because the tips are ordered differently between the OUCH tree and the OUwie tree: sapply(ape.tree$tip.label, function(n) which(ou.tree@nodelabels[90:195]==n)) %>% as.numeric)
pseudoinverse(t(ouch.W.1)%*%pseudoinverse(ouch.V)%*%ouch.W.1)%*%t(ouch.W.1)%*%pseudoinverse(ouch.V)%*%rev(x)
pseudoinverse(t(ouwie.W.1)%*%pseudoinverse(ouwie.V)%*%ouwie.W.1)%*%t(ouwie.W.1)%*%pseudoinverse(ouwie.V)%*%x

## This also works in the other direction - if you plug the OUwie matrices into the OUCH calculation, you get the OUCH result back out again
dat <- do.call(c,lapply(ou.dat,function(y)y[ou.tree@term]))
glssoln(ouch.W.1,dat,ouch.V)$coeff
glssoln(ouwie.W.1,rev(dat),ouwie.V)$coeff

## But these two do not agree. This suggests that the problem is with what OUwie is doing internally to the phenotypic data.
glssoln(ouwie.W.1,rev(dat),ouwie.V)$coeff
glssoln(ouwie.W.1,x,ouwie.V)$coeff 

## But this is not true for the tworegime case
pseudoinverse(t(ouch.W.2)%*%pseudoinverse(ouch.V)%*%ouch.W.2)%*%t(ouch.W.2)%*%pseudoinverse(ouch.V)%*%rev(x)
ouwie.theta.2



## Weight matrices for the single regime case are identical
ouch.W.1

## the fact that the theta values are different, however, suggests that the covariance matrices are not being calculated identically...
ouch.theta.1
ouwie.theta.1

## ... because if you plug the ouwie weight and covariance matrices into the ouch algorithm, the theta values are the same
dat <- do.call(c,lapply(ou.dat,function(y)y[ou.tree@term]))
glssoln(ouwie.W.1,rev(dat),ouwie.V)$coeff
glssoln(ouwie.W.2,rev(dat),ouwie.V)$coeff
ouwie.theta.1
ouwie.theta.2

## Weight matrices for the two-regime case appear different
ouch.W.2
ouwie.W.2
## but the tips are in a different order between OUCH and OUwie, and the regimes are also opposie (e.g., regime 1 in OUCH is regime 2 in OUwie).
## In particular, the tips are in exactly the opposite order between OUwie and OUCH

## A proper comparison shows that things are not the same, though they show similar patterns (e.g., the first 20 entries are 0 in both weight matrices, as are the last 7, and the values shift at the same places).
rev(ouch.W.2[,1])
ouwie.W.2[,2]
