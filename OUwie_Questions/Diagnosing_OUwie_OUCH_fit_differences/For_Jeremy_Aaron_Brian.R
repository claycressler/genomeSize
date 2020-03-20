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
ape.OUM <- OUwie(ape.tree, ape.dat, model='OUM')

#########################################################################
#                                                                       #
#	COMPUTING THE WEIGHT AND COVARIANCE MATRICES FOR OUCH           #
#                                                                       #
#########################################################################
data <- ou.dat
tree <- ou.tree
regimes <- ou.reg
sqrt.alpha <- 1
sigma <- 1

## From inside hansen.R hansen() function
source("hansen.R")
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
regimes <- rep(regimes,nchar)
regimes <- lapply(regimes,function(x)x[tree@nodes])
ouch.regimes <- regimes ## store this because it is critical for calculating the weight and covariance matrices
beta <- regime.spec(tree,regimes)

## From inside hansen.R ou.lik.fn() function
alpha <- sym.par(sqrt.alpha)
sigma <- sym.par(sigma)
n <- length(dat)
ev <- eigen(alpha,symmetric=TRUE)
system("R CMD SHLIB weight-matrix.c")
dyn.load("weight-matrix.so")
ouch.w <- .Call("ouch_weights",object=tree,lambda=ev$values,S=ev$vectors,beta=beta)
system("R CMD SHLIB covar-matrix.c")
dyn.load("covar-matrix.so")
ouch.v <- .Call("ouch_covar",object=tree,lambda=ev$values,S=ev$vectors,sigma.sq=sigma)
## This executes eqn (A8) in Butler and King 2004
source("glssoln.R")
gsol <- try(
    glssoln(ouch.w,dat,ouch.v),
    silent=FALSE
)
ouch.theta <- gsol$coeff

#########################################################################
#                                                                       #
#	COMPUTING THE WEIGHT AND COVARIANCE MATRICES FOR OUwie          #
#                                                                       #
#########################################################################

## Here are the key functions
source("varcov.ou.R")
source("weight.mat.R")
phy <- ape.tree
data <- ape.dat

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
ouwie.V<-varcov.ou(phy, edges, Rate.mat, root.state=1, simmap.tree=FALSE, scaleHeight=TRUE)
ouwie.W<-weight.mat(phy, edges, Rate.mat, root.state=1, simmap.tree=FALSE, scaleHeight=TRUE, assume.station=TRUE)
## here is the calculation of the thetas
diag(ouwie.V) <- diag(ouwie.V)+p[length(p)] ## error is estimated
ouwie.theta<-pseudoinverse(t(ouwie.W)%*%pseudoinverse(ouwie.V)%*%ouwie.W)%*%t(ouwie.W)%*%pseudoinverse(ouwie.V)%*%x

## Creates an ouchtree object.
##
## 'waittime' must be specified by the user, and is the
## expected time to a branching event, assuming the branching
## follows a Poisson process (that is, that branching events
## are exponentially distributed)
##
## 'endtime' must also be specified and determines how long to
## let the branching process proceed.
##
## The smaller the value of waittime, or the larger the value
## of endtime, the larger the tree that will be produced. If
## only these two variables are set, the tree that is produced
## will be ultrametric; specifically, it will be a Yule tree,
## corresponding to a pure birth process.
##
## 'time' can also be specified, setting the time of the root
## node.

## A birth-death branching process can be simulated by setting
## the value of extrate, which determines the rate at which
## lineages go extinct. This rate is also assumed to be
## exponentially distributed.
##
gen.tree <- function(waittime, endtime,
                     time = 0, extrate = 0,
                     anc.node = NA, tree) {
    if (missing(tree))
        tree <- matrix(NA,0,3,
                       dimnames=list(NULL,c("node","ancestor","time")))
    curr.node <- nrow(tree)+1
    if (time < endtime) {
        tree <- rbind(tree,c(curr.node,anc.node,time))
        extinct.time <- if (extrate > 0) {
            time+rexp(2,rate=extrate)
        } else {
            c(Inf,Inf)
        }
        time <- pmin(time+rexp(2,rate=1/waittime),endtime)
        for (k in 1:2) {
            if (time[k] < extinct.time[k]) {
                tree <- Recall(time=time[k],
                               waittime=waittime,endtime=endtime,extrate=extrate,anc.node=curr.node,tree=tree)
            } else {
                tree <- rbind(tree,c(nrow(tree)+1,curr.node,extinct.time[k]))
            }
        }
    } else {
        tree <- rbind(tree,c(curr.node,anc.node,time))
    }
    return(as.data.frame(tree))
}

## Generate phenotypic data according to an Ornstein-Uhlenbeck
## process, given an ouchtree object specified by tree, regimes
## specified by regimes, and OU parameters theta, alpha, and sigma.
## This is done by recursion, which gen.phenotypic.values acting as
## the "parent" function that calls set.phenotypic.values, the
## function that is recursed over.
gen.phenotypic.values <- function(tree, regimes, theta, alpha, sigma, remove.anc=TRUE) {
  if (class(tree) == 'data.frame')
    tree <- with(tree, ouchtree(node,ancestor,time))
  if (tree@nnodes != length(regimes))
    stop('Length of regime specification != number of nodes in tree; each node in tree must be assigned to a regime')
  if (length(theta) != length(levels(regimes)))
    stop('Each regime must have a unique theta value specified as the selective optima')

  root.reg <- as.numeric(regimes[1])
  root.pheno <- calc.root.pheno(theta[root.reg], alpha, sigma)
  pheno <<- vector(mode='numeric', length=tree@nnodes)
  set.phenotypic.values(curr.node=1, anc.node=0, anc.pheno=root.pheno, tree=tree, regimes=regimes, theta=theta, alpha=alpha, sigma=sigma)

  ## if T, set the phenotypic values of internal nodes to NA
  if (remove.anc==T)
    pheno[as.numeric(setdiff(tree@nodes, tree@term))] <<- NA

  pheno
}
## Calculate the phenotypic value of the ancestral node as being drawn
## from the stationary distribution of the OU process.
calc.root.pheno <- function(theta, alpha, sigma) {
  ## mean and variance of stationary distribution
  mean <- theta
  var <- sigma^2/(2*alpha)
  rnorm(1, mean=mean, sd=sqrt(var))
}
## This function recursively sets the phenotypic values of the nodes
## in an ouchtree object specified by tree.  The values are, by
## default, set for a vector called 'pheno'
set.phenotypic.values <- function(curr.node, anc.node, anc.pheno, tree, regimes, theta, alpha, sigma) {
  ## set the phenotype of the current node
  if (curr.node == 1) {
    pheno[curr.node] <<- anc.pheno

    ## get the descendents of the root node and set their phenotypes
    desc <- which(tree@ancestors==curr.node)
    Recall(curr.node=desc[1], anc.node=curr.node, anc.pheno=anc.pheno, tree=tree, regimes=regimes, theta=theta, alpha=alpha, sigma=sigma)
    Recall(curr.node=desc[2], anc.node=curr.node, anc.pheno=anc.pheno, tree=tree, regimes=regimes, theta=theta, alpha=alpha, sigma=sigma)
  }
  else {
    ## generate a new phenotypic value
    x0 <- anc.pheno
    t <- as.numeric(tree@times[curr.node])-as.numeric(tree@times[anc.node])
    th <- theta[as.numeric(levels(regimes)[regimes[curr.node]])]
    a <- alpha
    p <- x0*exp(-a*t) + th*(1-exp(-a*t)) +
      sigma*sqrt((1-exp(-2*a*t))/(2*a))*rnorm(1, 0, 1)

    ## set the value of the phenotype
    pheno[curr.node] <<- p

    ## if the curr.node is not a tip, get its ancestors and set their
    ## phenotypes
    if (!(curr.node %in% tree@term)) {
      desc <- which(tree@ancestors==curr.node)

      ## set the phenotypes for all descendents
      for (n in 1:length(desc))
        Recall(curr.node=desc[n], anc.node=curr.node, anc.pheno=p, tree=tree,
               regimes=regimes, theta=theta, alpha=alpha, sigma=sigma)
    }
  }
}
