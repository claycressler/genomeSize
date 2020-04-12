## The purpose of this code is to compare how OUwie and OUCH and calculate the covariance and weight matrices.
require(ouch)

#########################################################################
#                                                                       #
#	     GENERATE A THREE-TIP TREE WITH TWO REGIMES                 #
#                                                                       #
#########################################################################
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

set.seed(7787358)
tree <- gen.tree(waittime=0.8, endtime=1)
labels <- rep('', nrow(tree))
labels[which(tree$time==1)] <- paste0('Sp',1:sum(tree$time==1))
tree$labels <- labels
ou.tree <- ouchtree(nodes=tree$node, ancestors=tree$ancestor, times=tree$time, labels=tree$labels)
ou.reg <- paint(ou.tree, subtree=c('1'='1','2'='2'), branch=c('1'='1', '2'='2'))
plot(ou.tree, regimes=ou.reg)

#########################################################################
#                                                                       #
#	      COMPUTING THE WEIGHT AND COVARIANCE MATRICES              #
#                                                                       #
#########################################################################

## From hansen.R
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

sym.par <- function (x) {
  nchar <- floor(sqrt(2*length(x)))
  if (nchar*(nchar+1)!=2*length(x)) {
    stop("a symmetric matrix is parameterized by a triangular number of parameters",call.=FALSE)
  }
  y <- matrix(0,nchar,nchar)
  y[lower.tri(y,diag=TRUE)] <- x
  y%*%t(y)
}

sets.of.regimes <- function (object, regimes) {
  lapply(regimes,function(x)sort(unique(x)))
}

## Parameter estimates
sqrt.alpha <- 1
sigma <- 1

## From inside hansen.R hansen() function
## Compute beta for the two-regime case
nm <- deparse(substitute(ou.reg))[1]
regimes <- list(ou.reg)
names(regimes) <- nm
regimes <- rep(regimes,1)
regimes <- lapply(regimes,function(x)x[ou.tree@nodes])
beta <- regime.spec(ou.tree,regimes)


## From inside hansen.R ou.lik.fn() function
alpha <- sym.par(sqrt.alpha)
sigma <- sym.par(sigma)
ev <- eigen(alpha,symmetric=TRUE)

## Compute the weight matrix for the single regime and two-regime cases
if(file.exists("weight-matrix.so")) {## rebuild C executable on this computer
    system("rm weight-matrix.so")
    system("rm weight-matrix.o")
}
system("R CMD SHLIB weight-matrix.c")
dyn.load("weight-matrix.so")
ouch.W <- .Call("ouch_weights",object=ou.tree,lambda=ev$values,S=ev$vectors,beta=beta)

## Compute the covariance matrix, which does not depend on regime paintings
if(file.exists("covar-matrix.so")) {## rebuild C executable on this computer
    system("rm covar-matrix.so")
    system("rm covar-matrix.o")
}
system("R CMD SHLIB covar-matrix.c")
dyn.load("covar-matrix.so")
ouch.V <- .Call("ouch_covar",object=ou.tree,lambda=ev$values,S=ev$vectors,sigma.sq=sigma)

print("Weight matrix:")
print(ouch.W)
print("Covariance matrix:")
print(ouch.V)

exp.W <- matrix(0, ncol=2, nrow=3)
## Okay, from Aaron I understand this now. The regime at the root is critical and must be taken into account, because regardless of what regime a tip ends up in, we assume it has spent an infinite amount of time in the root regime. So, for this case where the root is in regime 1, and Sp1 and Sp2 have spent T=1 time units evolving under regime 2, the weight of regime 1 is:
exp.W[1:2,1] <- exp(-1) + exp(-1) * (exp(0)-exp(0))
## This can be understood in the following way: the infinite amount of time that these species spent in regime 1 is captured by exp(-1). From t=0, though, the influence of regime 1 wanes. This waning is captured by exp(-1) * (exp(0)-exp(0)).
## The weight for regime 2, on the other hand is just:
exp.W[1:2,2] <- exp(-1) * (exp(1)-exp(0))
## reflecting the increasing influence of regime 2
## Species 3, on the other hand has never been in any other regime, so the weight of regime 1 is
exp.W[3,1] <- exp(-1) + exp(-1) * (exp(1)-exp(0))
print("Expected weight matrix based on Eq A7")
print(exp.W)


## But what is OUwie doing?
library(OUwie)
library(phytools)
source("convert.R") ## for converting ouchtree to 'phylo' format
ape.tree <- convert(ot=ou.tree)
## Label the internal nodes with their regime
## figure out which clade is in regime 2
## root is ntips+1, find its daughters
daughters <- which(ape.tree$edge[,1]==ou.tree@nterm+1)
## which daughter corresponds to node 2 in the ouchtree?
node2 <- which(round(ape.tree$edge.length[daughters],4)==round(ou.tree@times[2],4))
## what is the nodelabel for this daughter
start <- ape.tree$edge[daughters[node2],2]
## what is the nodelabel for the other daughter?
finish <- ape.tree$edge[daughters[-node2],2]
## create node.labels
nodes <- (ou.tree@nterm+1):(ape.tree$Nnode+ou.tree@nterm)
regs <- rep('1', length(nodes))
reg2 <- start:(finish-1)
regs[which(nodes%in%reg2)] <- '2'
ape.tree$node.label <- regs

## Verifies that the ape.tree looks like the ou.tree
plot(ape.tree)
nodelabels(pch=21, bg=ape.tree$node.label)

## Generate some fake data
data <- as(ou.tree, 'data.frame')[,c(4,1,2,3)]
data$val <- c(NA,NA,-1,-1,1)
data$reg <- ou.reg
data <- data[!is.na(data$val),c('labels','reg','val')]

## the files you need (I've included a print() statement in OUwie.R to print the weight matrix)
source("OUwie.R")
source("weight.mat.R")
source("varcov.ou.R")
library(corpcor)

OUwie(ape.tree, data, model='OUM', root.station=TRUE, scaleHeight=TRUE)

## The initial weight matrix is
## [1] "weight matrix"
##          [,1]       [,2]
##[1,] 0.9420175 0.05798245
##[2,] 0.9420175 0.05798245
##[3,] 1.0000000 0.00000000

## Does this make sense?
## For OUwie, the amount of time in regime 2 is not from 0 but halfway between t=0 and the branching point between Sp1 and Sp2. Then the weight of regime 1 is:
exp.W[1:2,1] <- exp(-1) + exp(-1) * (exp(ou.tree@times[2]/2)-exp(0))
## This can be understood in the following way: the infinite amount of time that these species spent in regime 1 is captured by exp(-1). From t=0, though, the influence of regime 1 wanes. This waning is captured by exp(-1) * (exp(0)-exp(0)).
## The weight for regime 2, on the other hand is just:
exp.W[1:2,2] <- exp(-1) * (exp(1)-exp(ou.tree@times[2]/2))
## reflecting the increasing influence of regime 2
## Species 3, on the other hand has never been in any other regime, so the weight of regime 1 is
exp.W[3,1] <- exp(-1) + exp(-1) * (exp(1)-exp(0))
print(exp.W)
