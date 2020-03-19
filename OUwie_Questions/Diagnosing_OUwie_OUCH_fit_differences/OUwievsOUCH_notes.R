## This document goes through a lot of ways that OUwie and OUCH could be different from each other. It was created back in 2014. Here are the main findings:
## OUwie and OUCH create regime paintings in identical ways
## Though superficially different-looking, both OUwie and OUCH will calculate theta values and likelihoods identically, given the same weight and covariance matrices.


## My supposition is that OUCH and OUwie are not filling in the weight matrix used for calculating the MLE equivalently, which is why the two do not produce analogous results. Aaron gives an equation for the entries of the weight matrix (A7), so I can compare the weight matrices calculated by OUCH and OUwie and compare them against the formula in Bulter and King 2004.

###########
# Set-up for running OUwie on the salamander tree, from OUwievsOUCH.R
###########

require(OUwie)

############
# Trees  - one for each hypothesis, with the internal nodes labelled appropriately (below)
############
tree <- read.nexus("av_ultra_fulldataset.nex")
tree_xMpleth <- tree_xplethMpleth <- tree_xpleth <- tree_metamorphosis <- tree

############################################
# Regimes at Nodes for each hypothesis
############################################
nodedata <- read.csv("ouwie_nodelabels.csv")  # each hypothesis is nodedata$xpleth, nodedata$xMpleth, etc.

## Node Colors for plots
nodecol <- nodedata
nodecol[-1] <- sapply(nodecol[-1], as.character)
nodecol[nodecol=="x"] <- "black"
nodecol[nodecol=="Mpleth"] <- "red"
nodecol[nodecol=="pleth"] <- "yellow"
nodecol[nodecol=="M"] <- "green"

## code Regime names as numbers for the OUwie algorithm
nn <- nodedata
nn[-1] <- sapply(nn[-1], as.character)
nn[nn=="x"] <- 1
nn[nn=="Mpleth"] <- 2
nn[nn=="pleth"] <- 3
nn[nn=="M"] <- 4

tree_xMpleth$node.label <- nn$xMpleth
tree_metamorphosis$node.label <- nn$metamorphosis

############################################################
# Tip Data (species labels, genome size and regimes at tips)
############################################################
dat <- read.csv("tree1.dat.csv")
dat$genomesize <- log(dat$genomesize)

tdat <- dat[!is.na(dat$genomesize),]  # tip data
ndat <- dat[is.na(dat$genomesize),]   # node data

oo <- sapply(tree$tip.label, function(x) grep(x, tdat$labels))    ## reorder tdat to match order of tip labels in phylo tree

tipdat <- tdat[oo,]
tipdat <- tipdat[c("labels", "xMpleth", "metamorphosis", "genomesize")]

## Tip Colors for plots
tipcol <- tipdat
tipcol[2:3] <- sapply(tipcol[2:4], as.character)
tipcol[tipcol=="x"] <- "black"
tipcol[tipcol=="Mpleth"] <- "red"
tipcol[tipcol=="M"] <- "green"

## code Regime names as numbers for the OUwie algorithm
tt <- tipdat
tt[2:3] <- sapply(tt[2:3], as.character)
tt[tt=="x"] <- 1
tt[tt=="Mpleth"] <- 2
tt[tt=="M"] <- 4

### input dataframes - one for each hypothesis
data.xMpleth <- tt[c("labels", "xMpleth", "genomesize")]
data.metamorphosis <- tt[c("labels", "metamorphosis", "genomesize")]

###############################################################
###  xMpleth
###############################################################
plot(tree_xMpleth, cex=.6)     # plot tree with slightly reduced font size
nodelabels(pch=21, bg=nodecol$xMpleth)   # plot node labels
tiplabels(pch=21, bg=tipcol$xMpleth)   # plot tip labels
title("xMpleth")

## Running OUwie, where the model estimates separate optima for each regime, but a single alpha and sigma for the entire tree.
xMpleth.OUM <- OUwie(tree_xMpleth,data.xMpleth,model=c("OUM"),root.station=TRUE)

## Running OUwie with a single-rate BM model
xMpleth.BM1 <- OUwie(tree_xMpleth,data.xMpleth,model=c("BM1"),root.station=TRUE)


##############################
#
#	OUCH analysis
#
##############################
require(ape)
require(ouch)

# setup
dat <- read.csv("tree1.dat.csv")
tree <- read.nexus("av_ultra_fulldataset.nex")
tree <- ape2ouch(tree)
rownames(dat) <- dat$nodes
gsz <- log(dat['genomesize'])
regimes <- dat[c("xMpleth")]
regimes <- unlist(regimes)
names(regimes) <- rownames(dat)
ouch.ou <- hansen(gsz, tree, regimes, sqrt.alpha=1, sigma=1)
ouch.bm <- brown(gsz, tree)


# Plot the xMpleth hypothesis
quartz()
plot(tree, regimes=regimes['xMpleth'], cex=.35)

## One possibility, just from looking at the plotted trees, for why the two methods do not produce analogous results is because they differ in dealing with regime switches. It is unclear from the plots where the regime shifts are occurring on each tree: OUwie colors the node, whereas OUCH colors the branch.

## Regardless, for the xMpleth hypothesis, the initial weight matrix is calculated in the following way for each method
## OUwie:
#xMpleth.OUM <- OUwie(tree_xMpleth,data.xMpleth,model=c("OUM"),root.station=FALSE)
# Taken from the source code OUwie.R
phy <- tree_xMpleth
data <- data.xMpleth
model <- 'OUM'
root.station <- TRUE
simmap.tree <- FALSE
scaleHeight <- FALSE
lb=0.000001
ub <- 1000
clade <- NULL
mserr <- 'none'
diagn <- FALSE
## Make sure the data is in the same order as the tip labels
data<-data.frame(data[,2], data[,3], row.names=data[,1])
data<-data[phy$tip.label,]
## Values to be used throughout
n=max(phy$edge[,1])
ntips=length(phy$tip.label)
## Obtain a list of all the regime states. This is a solution for
## instances when tip states and the internal nodes are not of equal
## length:
tot.states <- factor(c(phy$node.label,as.character(data[,1])))
k <- length(levels(tot.states))
int.states <- factor(phy$node.label)
phy$node.label <- as.numeric(int.states)
tip.states <- factor(data[,1])
data[,1] <- as.numeric(tip.states)
## Obtain root state and internal node labels
root.state <- phy$node.label[1]
int.state <- phy$node.label[-1]
## New tree matrix to be used for subsetting regimes
require(phytools)
edges <- cbind(c(1:(n-1)),phy$edge,nodeHeights(phy))
edges <- edges[sort.list(edges[,3]),]

mm<-c(data[,1],int.state)
regime <- matrix(0,nrow=length(mm),ncol=length(unique(mm)))
## Generates an indicator matrix from the regime vector
for (i in 1:length(mm)) {
    regime[i,mm[i]] <- 1
}
## Finishes the edges matrix
edges=cbind(edges,regime)
## Resort the edge matrix so that it looks like the original matrix order
edges=edges[sort.list(edges[,1]),]

## The first column of 'edges' is an alternative numbering scheme for
## the branches that progress from the root towards the branches
## recursively.
## The second and third columns give ancestor-descendent node pairs
## (i.e. branches)
## The fourth and fifth columns give the node times for the ancestor
## and descendent (i.e. branch lengths)
## The sixth and seventh columns are indicator variables for which
## regime each branch belongs in

## From this, I can answer my question (raised above) about how OUwie
## is "painting" regimes onto the tree. In particular, looking the
## clade of the four Desmognathus spp. on the OUCH tree, the branch
## leading to that clade is in the Mpleth regime.  From the following
## plot,
quartz()
plot(tree_xMpleth, cex=0.4)
nodelabels()
## the equivalent branch is the one connecting nodes 145 and 146 (on
## the OUCH tree, this is the branch connecting nodes 52 and 51; nodes
## 145 and 52 are at time 0.603723793, and nodes 146 and 51 are at
## time 0.711690577, confirming). OUwie assigns this branch into
## regime 2 (the Mpleth regime).
## In other words, THE PAINTINGS ARE IDENTICAL FOR OUWIE AND OUCH.

bool=root.station
Rate.mat <- matrix(1, 2, k)

## Phenotypic data matrix
x<-as.matrix(data[,2])
## Generate the appropriate parameter matrix structure
index.mat<-matrix(0,2,k)
np=2 ## Two parameters, alpha and sigma
index<-matrix(TRUE,2,k) ## create a 2x2 matrix for the two parameter and two regimes
index.mat[1,1:k]<-1
index.mat[2,1:k]<-2
if(root.station==TRUE){
    param.count<-np+k ## 4 parameters
}
bool=root.station
Rate.mat <- matrix(1, 2, k)
## Make an initial parameter guess (the same as the one I will use below to test OUCH)
p <- c(1, 1)

## Compute the likelihood of these parameters (using the 'dev' function in OUwie.R)
source('~/Downloads/OUwie/R/weight.mat.R')
source('~/Downloads/OUwie/R/varcov.ou.R')
Rate.mat[] <- c(p, 1e-10)[index.mat] ## same alpha and sigma for both regimes
N <- length(x[,1]) ## number of data
## compute the variance-covariance and weight matrices from Butler and King 2004
V <- varcov.ou(phy, edges, Rate.mat, root.state=root.state)
W <- weight.mat(phy, edges, Rate.mat, root.state=root.state, assume.station=bool)
theta<-Inf
require(corpcor)
try(theta<-pseudoinverse(t(W)%*%pseudoinverse(V)%*%W)%*%t(W)%*%pseudoinverse(V)%*%x, silent=TRUE)
DET<-determinant(V, logarithm=TRUE)

logl<--.5*(t(W%*%theta-x)%*%pseudoinverse(V)%*%(W%*%theta-x))-.5*as.numeric(DET$modulus)-.5*(N*log(2*pi))

################################################
## OUCH does something different
dat <- read.csv("tree1.dat.csv")
tree <- read.nexus("av_ultra_fulldataset.nex")
tree <- ape2ouch(tree)
rownames(dat) <- dat$nodes
regimes <- dat[c("xMpleth")]
regimes <- unlist(regimes)
names(regimes) <- rownames(dat)
gsz <- log(dat['genomesize'])

#hansen(gsz, tree, regimes, sqrt.alpha=1, sigma=1)
sqrt.alpha=1
sigma=1
data <- gsz
nm <- rownames(data)
data <- lapply(as.list(data),function(x){names(x)<-nm;x})
nchar <- length(data)
data <- lapply(data,function(x)x[tree@nodes])
dat <- do.call(c,lapply(data,function(y)y[tree@term]))
nsymargs <- nchar*(nchar+1)/2
nalpha <- length(sqrt.alpha)
nsigma <- length(sigma)
nm <- deparse(substitute(regimes))[1]
regimes <- list(regimes)
names(regimes) <- nm
regime.spec <- function (object, regimes) {
    nterm <- object@nterm
    nchar <- length(regimes)
    reg <- sets.of.regimes(object,regimes)
    nreg <- sapply(reg,length)
    beta <- vector(mode='list',length=nterm)
    for (i in seq_len(nterm)) {
        p <- object@lineages[[object@term[i]]]
        np <- length(p)
        beta[[i]] <- vector(mode='list',length=nchar)
        for (n in seq_len(nchar)) {
            beta[[i]][[n]] <- matrix(data=NA,nrow=np,ncol=nreg[n])
            for (ell in seq_len(nreg[n])) {
                beta[[i]][[n]][,ell] <- ifelse(regimes[[n]][p]==reg[[n]][ell],1,0)
            }
        }
    }
    beta
}
sets.of.regimes <- function (object, regimes) {
  lapply(regimes,function(x)sort(unique(x)))
}
## for each tip species, give the history of regimes from the root to that tip. This returns a list, with one entry per tip. Each entry is a matrix with row number equal to the number of nodes between root and tip and column number equal to the number of regimes. The entries of this matrix are indicator variables (0,1) indicating which regime the branch was in.
beta <- regime.spec(tree,regimes)

optim.diagn <- vector(mode='list', length=0)
sym.par <- function (x) {
  nchar <- floor(sqrt(2*length(x)))
  if (nchar*(nchar+1)!=2*length(x)) {
    stop("a symmetric matrix is parameterized by a triangular number of parameters",call.=FALSE)
  }
  y <- matrix(0,nchar,nchar)
  y[lower.tri(y,diag=TRUE)] <- x
  y%*%t(y)
}

alpha <- sym.par(sqrt.alpha)
sigma <- sym.par(sigma)
## compute the likelihood using ou.lik.fn in hansen.R
n <- length(dat)
ev <- eigen(alpha,symmetric=TRUE)
w <- .Call('ouch_weights',object=tree,lambda=ev$values,S=ev$vectors,beta=beta,PACKAGE='ouch')
v <- .Call('ouch_covar',object=tree,lambda=ev$values,S=ev$vectors,sigma.sq=sigma,PACKAGE='ouch')
source('~/Dropbox/Research/ouch/ouch/R/glssoln.R')
gsol <- try(glssoln(w,dat,v),
            silent=FALSE)
e <- gsol$residuals
theta <- gsol$coeff
q <- e%*%solve(v,e)
det.v <- determinant(v,logarithm=TRUE)
dev <- n*log(2*pi)+as.numeric(det.v$modulus)+q[1,1]

#############################
## Okay, W != w and V != v, and logl != dev, so the two methods are doing something different. One way to spot check it would be to use the W and V matrices calculated using OUwie and plug those into the code for calculating the loglikelihood from OUCH, to see if I can pin down exactly where the two methods diverge from one another.
gsol <- try(glssoln(W,dat,V),
            silent=FALSE)
e <- gsol$residuals
theta <- gsol$coeff
q <- e%*%solve(V,e)
det.v <- determinant(V,logarithm=TRUE)
dev <- n*log(2*pi)+as.numeric(det.v$modulus)+q[1,1]
## Okay, the deviances still do not agree, so the methods disgree not only the calculation of the covariance and weight matrices, but also on how you use these matrices to calculate the thetas and likelihoods.

## Comparing the codes more carefully, they share a number of things in common:
## OUwie liklihood = -1/2*(STUFF + DET$modulus + N*log(2*pi))
## OUCH liklihood = q[1,1] + det.v$modulus + n*log(2*pi)
## so if STUFF, which = t(W%*%theta-x)%*%pseudoinverse(V)%*%(W%*%theta-x) is equal to q[1,1], which takes the residuals of the glssoln, then the two agree.
q[1,1]
t(W%*%theta-x)%*%pseudoinverse(V)%*%(W%*%theta-x)
## these do not agree, but since V, W, and x must be the same, the difference is coming through how each is calculating theta.
## OUwie calculates theta in the following way
theta<-pseudoinverse(t(W)%*%pseudoinverse(V)%*%W)%*%t(W)%*%pseudoinverse(V)%*%x
## where pseudoinverse uses singular value decomposition to find the inverse of a matrix. From the help() for pseudoinverse:
## Any rectangular real matrix M can be decomposed as
## M = U D V'
## where U and V are orthogonal, V' is t(V) and D is a diagonal matrix containing only the positive singular values (as determined by 'tol', see also 'fast.svd').
## The pseudoinverse is then obtained as
## iM = V D^(-1) U'.
## Note that tol = max(dim(M))*max(D)*.Machine$double.eps

## whereas OUCH uses glssoln, the main line being
#(s$v[,inds,drop=FALSE]%*%(svals %*% t(s$u[,inds,drop=FALSE])))%*%forwardsolve(vh,x,upper.tri=TRUE,transpose=TRUE)
## where
vh <- chol(v) ## v is the covariance matrix
s <- svd(forwardsolve(vh,w,upper.tri=TRUE,transpose=TRUE)) ## a is the weight matrix
tol <- sqrt(.Machine$double.eps)
inds <- s$d > tol*max(s$d)
svals <- s$d[inds,drop=FALSE]
r <- length(svals)
svals <- diag(1/svals,nrow=r,ncol=r)
## chol computes the Choleski factorization of a real symmetric positive-definite square matrix
## svd computes the singular-value decomposition of a rectangular matrix
## forwardsolve solves a triangular system of linear equations, where vh is an upper triangular matrix giving the coefficients for the system to be solved and a (the weight matrix) is a matrix whose columns give the right-hand sides for the equations. With transpose=TRUE, solves t(vh)%*%y==a, returning y

## There is a lot of similarity here, as pseudoinverse appears to be a function that performs many of the operations that Aaron is doing
vh <- chol(v) ## v is the covariance matrix
mat1 <- forwardsolve(vh,w,upper.tri=TRUE,transpose=TRUE)
t(vh)%*%mat1==w ## check: TRUE
round(mat1,3)==round(pseudoinverse(t(vh))%*%w,3) ## TRUE!
## so pseudoinverse(t(vh))%*%w gives the same answer as forwardsolve(vh,w,...)

## Then forwardsolve(vh,x,...)=pseudoinverse(t(vh),x)
## I wonder if V returned by OUwie is very similar to t(vh) in OUCH?? V appears in many places where t(vh) would appear in Aaron's code, if it were rewritten using pseudoinverse rather than forwardsolve.
## Taking that suggestion seriously, pseudoinverse(V)%*%x is equivalent to forwardsolve(vh,x,...) and pseudoinverse(t(W)%*%pseudoinverse(V)%*%W)%*%t(W) is equal to (s$v[,inds,drop=FALSE]%*%(svals %*% t(s$u[,inds,drop=FALSE])))

## Moreover, the singular value decomposition of a matrix X, as computed by svd, is X = U D V', and 'svd' returns a list with elements
## d: a vector containing the singular values of X
## u: a matrix whose columns contain the left singular vectors of X
## v: a matrix whose columns contain the right singular vectors of X
## OUCH then does, essentially,
s <- svd(forwardsolve(vh,a,upper.tri=TRUE,transpose=TRUE))
s$v %*% diag(s$d) %*% t(s$u)
## which should be identical to the return of the pseudoinverse, based on what I have above!
## So if s=svd(forwardsolve(t(vh),w)), then s$v %*% diag(1/s$d) %*% t(s$u)
## should be equal to pseudoinverse(forwardsolve(t(vh,w)))?
require(corpcor)
round(s$v %*% (diag(1/s$d) %*% t(s$u)),7)==round(pseudoinverse(forwardsolve(t(vh),w)),7) ## YES!

## Then
round(s$v %*% (diag(1/s$d) %*% t(s$u)),7)==round(pseudoinverse(pseudoinverse(t(vh))%*%w),7) ## YES!

## So Aaron's calculation of theta ('y' in glssoln.R) is equivalent to the following:
pseudoinverse(pseudoinverse(t(vh))%*%w)%*%(pseudoinverse(t(vh))%*%x)
## But what OUwie does for the calculation of theta is:
pseudoinverse(t(w)%*%pseudoinverse(v)%*%w)%*%t(w)%*%pseudoinverse(v)%*%x

## Plugging in the values above, I see that these are, in fact, equivalent!
## So the two methods agree on the calculation of the thetas, though they do it in different ways.
x <- dat
vh <- chol(v) ## v is the covariance matrix
s <- svd(forwardsolve(vh,w,upper.tri=TRUE,transpose=TRUE))
## Aaron's calculation of theta
theta <- (s$v%*%(diag(1/s$d)%*%t(s$u)))%*%forwardsolve(vh,x,upper.tri=TRUE,transpose=TRUE)
## The equivalent, using pseudoinverses
theta <- pseudoinverse(pseudoinverse(t(vh))%*%w)%*%(pseudoinverse(t(vh))%*%x)
## OUwie's calculation of theta is:
theta <- pseudoinverse(t(w)%*%pseudoinverse(v)%*%w)%*%t(w)%*%pseudoinverse(v)%*%x

## So the theta calculations agree. What about the likelihoods, given theta?
## Aaron takes the theta estimates and compute the residuals
e <- w%*%theta-x
dim(e) <- 106
## He then computes 'q' as
q <- e%*%solve(v,e)
## and uses q[1,1] in the calculation of the likelihood
q[1,1]

## In the equivalent spot of the likelihood calculation in OUwie, those authors uses
t(w%*%theta-x)%*%pseudoinverse(v)%*%(w%*%theta-x)

## These also agree. So, given the same weight and variance-covariance matrices, the two methods will calculate identical likelihood values.

#########
## SO IT APPEARS THAT ALL OF THE DIFFERENCES ARE DUE TO DIFFERENCES IN THE CALCULATION OF THE WEIGHT AND VARIANCE-COVARIANCE MATRICES.
#########

## To understand why this is, I need to look at Butler and King 2004, which reports the calculation of theta as
## pseduoinverse(t(W)%*%pseudoinverse(V)%*%W)%*%pseudoinverse(V)%*%x
## and the calculation of the likelihood as
## t(x-W%*%theta) %*% pseudoinverse(V) %*% (x-W%*%theta)
## so OUwie takes Butler and King at its literal word, whereas Aaron is doing something a little more cheeky with the code for hansen.

## YUCK. I was really hoping to avoid having to compare the codes calculating the variance-covariance and weight matrices.

## Easier than working with the salamander tree is to use a simpler tree, with a simple regime painting.
source('~/Dropbox/Research/OUCH_Simulation_study/treegen.R')
## Generates a tree with only 9 tips
set.seed(771201)
tree <- gen.tree(waittime=0.8, endtime=1)
labels <- rep('', nrow(tree))
labels[which(tree$time==1)] <- paste0('Sp',1:9)
tree$labels <- labels
ou.tree <- ouchtree(nodes=tree$node, ancestors=tree$ancestor, times=tree$time, labels=tree$labels)
plot(ou.tree)

## Generate a single clade in a new regime
ou.reg <- paint(ou.tree, subtree=c('1'='1','13'='2'), branch=c('1'='1','13'='2'))
plot(ou.tree, regimes=ou.reg)

## Generate phenotypic data
x <- gen.phenotypic.values(ou.tree, ou.reg, theta=c(-1,1), alpha=3, sigma=1)

require(pmc)
ape.tree <- convert(ot=ou.tree)
nodedata <- data.frame(node=seq(10,17), reg=c(rep('1',6),'2','2'))
nodedata$reg <- as.character(nodedata$reg)
ape.tree$node.label <- nodedata$reg

ape.dat <- as(tree, 'data.frame')[,c(4,1,2,3)]
ape.dat$val <- x
ape.dat$reg <- ou.reg
tdat <- ape.dat[!is.na(ape.dat$val),]
oo <- sapply(ape.tree$tip.label, function(x) grep(x, tdat$labels))    ## reorder tda
tipdat <- tdat[oo,]
tipdat <- tipdat[c('labels','reg','val')]
tipdat[2] <- as.vector(sapply(tipdat[2], as.character))
ape.dat <- tipdat

plot(ape.tree)
nodelabels(pch=21, bg=nodedata$reg)
tiplabels(pch=21, bg=ape.dat$reg)

## input dataframes
ou.dat <- data.frame(val=x)
## Check with Brownian motion
ape.BM1 <- OUwie(ape.tree, ape.dat, model='BM1')
ou.BM <- brown(ou.dat, ou.tree)

ape.OUM <- OUwie(ape.tree, ape.dat, model='OUM')
ou.OU <- hansen(ou.dat, ou.tree, ou.reg, sqrt.alpha=1, sigma=1)
## generating and fitting many different phenotypic datasets, I find that the two typically converge to similar AIC estimates, although ouch always has lower AIC, and the theta estimates are usually pretty similar, but the alpha and sigma estimates are often wildly different. This is actually pretty hard to understand, because the alpha estimate is used to determine the theta estimates, so it is unclear how the two can be so different, and yet produce similar thetas.


## Trying that
tree <- read.nexus("av_ultra_fulldataset.nex")
tree <- ape2ouch(tree)
dat <- read.csv("tree1.dat.csv")
## Put dat in the same order as the tree
dat <- dat[order(dat$nodes),]
## check that orders are identical
as.character(dat$labels[90:195])==tree@nodelabels[90:195]
rownames(dat) <- dat$nodes
gsz <- log(dat['genomesize'])
regimes <- dat[c("xMpleth")]
regimes <- unlist(regimes)
names(regimes) <- rownames(dat)
ouch.OU <- hansen(gsz, tree, regimes, sqrt.alpha=1, sigma=1)

## convert to an apetree from the ouchtree object (cannot use 'convert' from pmc because of polytomies)
ot <- tree
n <- ot@nnodes
anc <- as.integer(ot@ancestors[!is.na(ot@ancestors)])
## cannot use this, because it assumes binary tree
#internal <- sort(anc)[seq(2,n-1,by=2)]
internal <- which(!as.numeric(ot@nodes)%in%ot@term)
tmp <- integer(n)
tmp[internal] <- 1
tips <- which(tmp==0)
root <- which(is.na(ot@ancestors))
internal <- internal[internal!=root]
new_ids <- integer(n)
new_ids[tips] <- 1:length(tips)
new_ids[root] <- length(tips)+1
new_ids[internal] <- (length(tips)+2):(length(tips)+1+length(internal))
new_ancestor <- new_ids[as.integer(ot@ancestors)]
edge <- matrix(NA, n-1, 2)
edge[,1] <- new_ancestor[!is.na(new_ancestor)]
edge[,2] <- new_ids[!is.na(new_ancestor)]
anc <- as.integer(ot@ancestors[!is.na(ot@ancestors)])
lengths <- ot@times[!is.na(ot@ancestors)] - ot@times[anc]
labels <- ot@nodelabels[tips]
apetree <- list(edge=edge, Nnode=89, tip.label = labels, edge.length= lengths )
class(apetree) <- 'phylo'
## This apetree is identical to the ouchtree object, which can be confirmed by plotting:
plot(tree, cex=0.5)
quartz()
plot(apetree, cex=0.5)
## Another way to confirm
tree@nodelabels[90:195]==apetree$tip.label

ape.dat <- as(tree, 'data.frame')[,c(4,1,2,3)]
ape.dat$genomesize <- gsz
apereg <- as.character(dat$xMpleth)
apereg[apereg=='x'] <- 1
apereg[apereg=='Mpleth'] <- 2
ape.dat$xMpleth <- apereg
tdat <- ape.dat[!is.na(ape.dat$genomesize),]
ape.dat <- tdat[c('labels','xMpleth','genomesize')]
ape.dat$genomesize <- as.numeric(unlist(ape.dat$genomesize))
apetree$node.label <- ape.dat$xMpleth[1:89]

ouwie.oum <- OUwie(apetree, ape.dat, model='OUM', root.station=TRUE)
## Still doesn't work. However, I am wondering now whether it has something to do with the tree having polytomies.

## Try again with a larger simulated tree, just to see if we can get different data
set.seed(77120)
tree <- gen.tree(waittime=0.26, endtime=1)
labels <- rep('', nrow(tree))
labels[which(tree$time==1)] <- paste0('Sp',1:length(which(tree$time==1)))
tree$labels <- labels
ou.tree <- ouchtree(nodes=tree$node, ancestors=tree$ancestor, times=tree$time, labels=tree$labels)
plot(ou.tree)

## Generate a single clade in a new regime
ou.reg <- paint(ou.tree, subtree=c('1'='1','2'='2'), branch=c('1'='1','2'='2'))
plot(ou.tree, regimes=ou.reg)

## Generate phenotypic data
set.seed(1234)
x <- gen.phenotypic.values(ou.tree, ou.reg, theta=c(-1,1), alpha=1.5, sigma=0.5)

require(pmc)
ape.tree <- convert(ot=ou.tree)
nodedata <- data.frame(node=120:237, reg=c(rep('1',62),rep('2',56)))
nodedata$reg <- as.character(nodedata$reg)
ape.tree$node.label <- nodedata$reg

ape.dat <- as(ou.tree, 'data.frame')[,c(4,1,2,3)]
ape.dat$val <- x
ape.dat$reg <- ou.reg
tipdat <- ape.dat[!is.na(ape.dat$val),]
tipdat <- tipdat[c('labels','reg','val')]
tipdat[2] <- as.vector(sapply(tipdat[2], as.character))
ape.dat <- tipdat

plot(ape.tree)
nodelabels(pch=21, bg=nodedata$reg)
tiplabels(pch=21, bg=ape.data$reg)

## input dataframes
ou.dat <- data.frame(val=x)
## Check with Brownian motion
ape.BM1 <- OUwie(ape.tree, ape.dat, model='BM1')
ou.BM <- brown(ou.dat, ou.tree)

ape.OUM <- OUwie(ape.tree, ape.dat, model='OUM')
ou.OU <- hansen(ou.dat, ou.tree, ou.reg, sqrt.alpha=1, sigma=1)

## No, more likely it is because of other things. In any case, it
## explains what is happening much better than anything else. The
## algorithm is getting stuck on a boundary. The thing to do is to
## rewrite the code to use an unconstrained optimization. Since alpha
## and sigma are both positive, you can just log-transform them and
## then transform back.

## One way to confirm would be to choose different alpha and sigma values for the data generation and see if the algorithm converges with a different choice.

## Generate phenotypic data (increase alpha from 1.5 to 3 and reduce sigma to 0.25)
x <- gen.phenotypic.values(ou.tree, ou.reg, theta=c(-1,1), alpha=3, sigma=0.25)
ape.dat <- as(ou.tree, 'data.frame')[,c(4,1,2,3)]
ape.dat$val <- x
ape.dat$reg <- ou.reg
tipdat <- ape.dat[!is.na(ape.dat$val),]
tipdat <- tipdat[c('labels','reg','val')]
tipdat[2] <- as.vector(sapply(tipdat[2], as.character))
ape.dat <- tipdat
ou.dat <- data.frame(val=x)
## Check with Brownian motion
ape.BM1 <- OUwie(ape.tree, ape.dat, model='BM1')
ou.BM <- brown(ou.dat, ou.tree)
## OU
ape.OUM <- OUwie(ape.tree, ape.dat, model='OUM')
ou.OU <- hansen(ou.dat, ou.tree, ou.reg, sqrt.alpha=1, sigma=1)
## Nope, keeps getting stuck on the boundary.

## Obviously, for better or worse, OUwie doesn't work as well as OUCH. The question is whether this is because of a mistake with my code or whether it is an error with OUwie.
data(tworegime)
## Looking at the structure of the trait data, 'trait', the regimes are specified as an integer variable, not a character variable - could this be causing problems?
ape.dat[,2] <- as.integer(ape.dat[,2])
ape.OUM <- OUwie(ape.tree, ape.dat, model='OUM') ## Nope.

## Rewrite OUwie.R as OUwie_Clay.R, which allows for unconstrained optimization, but only for non-BM models with non-SIMMAP trees.
require(numDeriv)
require(phytools)
require(corpcor)
source('OUwie/R/varcov.ou.R')
source('OUwie/R/weight.mat.R')
source('OUwie/R/OUwie_Clay_unconstrained.R')
ape.OUM.1 <- OUwie.Clay(ape.tree, ape.dat, model='OUM', alpha=1, sigma.sq=1)
source('OUwie/R/OUwie_Clay_constrained.R')
ape.OUM.2 <- OUwie.Clay(ape.tree, ape.dat, model='OUM', alpha=2, sigma.sq=0.1)
## I wonder if the problem is with nloptr, as the optimization tends to run itself into a very low alpha value, regardless of whether I use the phytools-inspired initialization or I allow user-specified starts, and regardless of whether I use unconstrained or constrained optimization, and even if I start it with the true alpha and sigma values.

## Try a new version that uses 'optim' rather than 'nloptr'
source('OUwie/R/OUwie_Clay_unconstrained_optim.R')
ape.OUM.3 <- OUwie.Clay(ape.tree, ape.dat, model='OUM', alpha=1, sigma.sq=1)
source('OUwie/R/OUwie_Clay_constrained_optim.R')
ape.OUM.4 <- OUwie.Clay(ape.tree, ape.dat, model='OUM', alpha=1, sigma.sq=1)

## None of these changes make two figs worth of a difference. The problems persist. I think I am going to have to dive into the variance-covariance matrix and weight matrix calculations.

require(OUwie)
require(ouch)
## Example dataset for OUwie
data(tworegime)
select.reg<-character(length(tree$node.label))
select.reg[tree$node.label == 1] <- "black"
select.reg[tree$node.label == 2] <- "red"
ouwie.fit <- OUwie(tree,trait,model=c("OUM"),root.station=TRUE)

## Perform same analysis with ouch
outree <- ape2ouch(tree)
outree@nodelabels[1:63] <- as.character(1:63)
## Figure out which nodes subtend new regimes
plot(tree)
nodelabels(pch=21, bg=select.reg)
quartz()
plot(outree)

## I am still not entirely sure how ouwie paints regimes, so try four different possibilities
regimes1 <- paint(outree, subtree=c('1'=1, '32'=2), branch=c('1'=1, '32'=2))
regimes2 <- paint(outree, subtree=c('1'=1, '32'=2), branch=c('1'=1))
regimes3 <- paint(outree, subtree=c('1'=1, '32'=2), branch=c('32'=2))
regimes4 <- paint(outree, subtree=c('1'=1, '32'=2))


outrait = sapply(outree@nodelabels, function(x) ifelse(any(trait[,1]==x),trait[trait[,1]==x,3],NA))
names(outrait) <- outree@nodes
ouch.fit1 <- hansen(outrait, outree, regimes1, sqrt.alpha=1, sigma=1)
ouch.fit2 <- hansen(outrait, outree, regimes2, sqrt.alpha=1, sigma=1)
ouch.fit3 <- hansen(outrait, outree, regimes3, sqrt.alpha=1, sigma=1)
ouch.fit4 <- hansen(outrait, outree, regimes4, sqrt.alpha=1, sigma=1)
## estimates do not agree between ouwie and ouch.

## Can I get them to agree simply by removing the initialization step of OUwie, forcing it to use the same initial estimates as ouch is using?
source('OUwie/R/OUwie_noinit.R')
ouwie.fit.2 <- OUwie(tree, trait, model=c('OUM'), root.station=TRUE, alpha=1, sigma=1)
## No.

##

