## The purpose of this code is to compare how OUwie and OUCH and calculate the covariance and weight matrices.

require(OUwie)
require(phytools)
require(corpcor)
require(ouch)
require(ape)
require(magrittr)

## Load trees, regimes, phenotypic data
treelist <- readRDS("treelist_small_tree.RDS")
ou.tree <- treelist[["ou.tree"]] ## phylogeny in ouchtree format
ape.tree <- treelist[["ape.tree"]] ## phylogeny in phylo format
ou.reg <- treelist[["ou.reg"]] ## OU regime specification for ouch
ou.dat <- treelist[["ou.dat"]] ## phenotypic data in ouch format
ape.dat <- treelist[["ape.dat"]] ## phenotypic data in OUwie format

## Run both programs
ou.OU <- hansen(ou.dat, ou.tree, ou.reg, sqrt.alpha=1, sigma=1)
ouwie.OUM <- OUwie(ape.tree, ape.dat, model='OUM', root.station=TRUE, scaleHeight=TRUE)

## The two are pretty similar in what they converge on
summary(ou.OU)
ouwie.OUM

## BUT they don't agree on the calculation of likelihood for the same set of parameters
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
if(file.exists("weight-matrix.so")) {## rebuild C executable on this computer
    system("rm weight-matrix.so")
    system("rm weight-matrix.o")
}
system("R CMD SHLIB weight-matrix.c")
dyn.load("weight-matrix.so")
ouch.W.1 <- .Call("ouch_weights",object=ou.tree,lambda=ev$values,S=ev$vectors,beta=beta.1)
ouch.W.2 <- .Call("ouch_weights",object=ou.tree,lambda=ev$values,S=ev$vectors,beta=beta.2)

## Compute the covariance matrix, which does not depend on regime paintings
if(file.exists("covar-matrix.so")) {## rebuild C executable on this computer
    system("rm covar-matrix.so")
    system("rm covar-matrix.o")
}
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
data.1$reg[data.1$reg=="2"] <- "1"

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

## Recover the OUCH data
dat <- do.call(c,lapply(ou.dat,function(y)y[ou.tree@term]))

## Here, the ouch and ouwie results are identical for the single-regime case
ouwie.theta.1
ouch.theta.1

## BUT they are different for the two-regime case
ouwie.theta.2
ouch.theta.2

## AND the covariance matrices appear to be different
ouch.V
ouwie.V

## HOWEVER... even though the covariance matrices are different, plugging either of them into eq. A8 in Butler and King (2004) produces identical results, regardless of whether you use the OUwie data or ouch data. This suggests that the covariance matrices are actually the same, even though they appear to be different.
glssoln(ouwie.W.1,x,ouwie.V)$coeff
glssoln(ouch.W.1,x,ouch.V)$coeff
glssoln(ouwie.W.1,dat,ouwie.V)$coeff
glssoln(ouch.W.1,dat,ouch.V)$coeff

## And they are, up to a constant: if I subtract a value from every entry in the ouch covariance matrix, I get the OUwie covariance matrix
round(ouch.V-ouch.V[1,20], 5)==round(ouwie.V,5)

## BUT, the weight matrices for the two-regime case are different
ouch.W.2
ouwie.W.2

## since the data looks the same for both OUwie and ouch, if you use the ouch weight matrix, you get back the ouch result, regardless of whether you use the OUwie data and/or the OUwie covariance matrix
glssoln(ouch.W.2,x,ouwie.V)$coeff ## ouch weight, OUwie data, OUwie covariance
glssoln(ouch.W.2,dat,ouwie.V)$coeff ## ouch weight, ouch data, OUwie covariance
ouch.theta.2

## but if you use the OUwie weight matrix, you get back the OUwie result, regardless of whether you use the ouch data and/or the ouch covariance matrix
glssoln(ouwie.W.2,dat,ouch.V)$coeff
glssoln(ouwie.W.2,x,ouch.V)$coeff
ouwie.theta.2

## This suggests that the two programs are definitely calculating the weight matrices differently.
## Looking at the phylogeny, for Sp 20 there are 7 branches back to the root. From this
ou.tree@epochs[[1]]
## we can easily calculate by hand the entry for Sp20 in Eq. A7 of Butler and King. If epochs = branch times, then it would seem like the entry for W[20,1] should be
exp(-1) *
    (rev(exp(ou.tree@epochs[[20]]))[2]-rev(exp(ou.tree@epochs[[20]]))[1] +
        rev(exp(ou.tree@epochs[[20]]))[3]-rev(exp(ou.tree@epochs[[20]]))[2] +
            rev(exp(ou.tree@epochs[[20]]))[4]-rev(exp(ou.tree@epochs[[20]]))[3] +
                rev(exp(ou.tree@epochs[[20]]))[5]-rev(exp(ou.tree@epochs[[20]]))[4] +
                    rev(exp(ou.tree@epochs[[20]]))[6]-rev(exp(ou.tree@epochs[[20]]))[5] +
                        rev(exp(ou.tree@epochs[[20]]))[7]-rev(exp(ou.tree@epochs[[20]]))[6])
## (Note: this is = exp(-1)*(rev(exp(ou.tree@epochs[[20]]))[7]-rev(exp(ou.tree@epochs[[20]]))[1]) = exp(-1)*(exp(1)-exp(0)) = exp(0)-exp(-1)
## Thus, for species 8-20, W[8:20,1] = exp(0)-exp(-1), and W[8:20,2] = 0 because none of these species ever experience regime 2.

## For Sp 1-7, there is one uncertainty: when does the regime change from regime 1 to regime 2? It's actually ambiguous for this tree because the regimes go back to the root. Since regime 1 is designated as the ancestral regime, the regime changes somewhere along the first branch (based on the epochs, sometime between t=0 and 0.3988038). If you assume the regime switches at zero, then the entry in the weight matrix for species 1 in regime 1 would clearly be 0. But
ouch.W.2[1,1]
## is clearly non-zero, so the regime switch is happening somewhere along the branch. If you assume that the regime switch happens at the next branching point (t=0.3988038), then the entry in the weight matrix (W[1,1]) should be
exp(-1) * (exp(0.3988038) - exp(0))
## And the entry for W[1,2] should be
exp(-1) * (exp(1) - exp(0.3988038))
## because the rest of the phylogenetic history of species 1 is spent in regime 2.

## Thus the expected weight matrix is
exp.W = data.frame(
    "1"=c(rep(exp(-1) * (exp(0.3988038) - exp(0)), 7),
        rep(exp(-1) * (exp(1) - exp(0)), 13)),
    "2"=c(rep(exp(-1) * (exp(1) - exp(0.3988038)), 7),
        rep(0, 13)))
## If you restandardize by dividing by the rowSum,
exp.W/rowSums(exp.W)
## you end up with something that is not the same as *either* the ouch weight matrix OR the OUwie weight matrix (although Sp 8-20 match the ouch weight matrix
ouch.W.2
ouwie.W.2

## What if, instead of assuming the regime switches at the branch point, you assume the switch from regime 1 to 2 happens halfway down the branch?
exp.W = data.frame(
    "1"=c(rep(exp(-1) * (exp(0.3988038/2) - exp(0)), 7),
        rep(exp(-1) * (exp(1) - exp(0)), 13)),
    "2"=c(rep(exp(-1) * (exp(1) - exp(0.3988038/2)), 7),
        rep(0, 13)))
exp.W/rowSums(exp.W)
## now the weights for species 1-7 match the OUwie weight matrix, whereas the weights for species 8-20 match the ouch weight matrix!?
OUwie.W.2
ouch.W.2

## My last note here is that the OUCH weight matrix entries for Sp 1-7 in regime 1 are equal to exp(-1) and in regime 2 are equal to exp(0)-exp(-1), which seems oddly coincidental to me.
