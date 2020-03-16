##############################
#
#	OUwie analysis
#
##############################

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


## Setting root.station=FALSE causes OUwie to estimate theta_0; setting root.station=TRUE drops theta_0 from the model and assumes that it is distributed according to the stationary distribution of the OU process
xMpleth.OUMVA <- OUwie(tree_xMpleth,data.xMpleth,model=c("OUMVA"),root.station=FALSE)
xMpleth.OUMV <- OUwie(tree_xMpleth,data.xMpleth,model=c("OUMV"),root.station= FALSE)
xMpleth.OUM <- OUwie(tree_xMpleth,data.xMpleth,model=c("OUM"),root.station=FALSE)
xMpleth.OUM.2 <- OUwie(tree_xMpleth,data.xMpleth,model=c("OUM"),root.station=TRUE)

xMpleth.BMS <- OUwie(tree_xMpleth,data.xMpleth,model=c("BMS"),root.station= FALSE)
xMpleth.BM1 <- OUwie(tree_xMpleth,data.xMpleth,model=c("BM1"),root.station= FALSE)
xMpleth.OUMV.T <- OUwie(tree_xMpleth,data.xMpleth,model=c("OUMV"),root.station= TRUE)

xMpleth.OUMV    ## Display results of OUwie OUM model
xMpleth.OUMV.T  ## Display results root station=T

###############################################################
###  metamorphosis
###############################################################
quartz()
plot(tree_metamorphosis, cex=.6)     # plot tree with slightly reduced font size
nodelabels(pch=21, bg=nodecol$metamorphosis)   # plot node lables
tiplabels(pch=21, bg=tipcol$metamorphosis)   # plot tip labels
title("metamorphosis")

metamorphosis.OUMVA <- OUwie(tree_metamorphosis,data.metamorphosis,model=c("OUMVA"),root.station=FALSE)
metamorphosis.OUMV <- OUwie(tree_metamorphosis,data.metamorphosis,model=c("OUMV"),root.station= FALSE)
metamorphosis.OUM <- OUwie(tree_metamorphosis,data.metamorphosis,model=c("OUM"),root.station= FALSE)
metamorphosis.BMS <- OUwie(tree_metamorphosis,data.metamorphosis,model=c("BMS"),root.station= FALSE)
metamorphosis.BM1 <- OUwie(tree_metamorphosis,data.metamorphosis,model=c("BM1"),root.station= FALSE)

metamorphosis.OUM <- OUwie(tree_metamorphosis,data.metamorphosis,model=c("OUM"),root.station= TRUE)

xMpleth.OUM.2
# Fit
#       -lnL      AIC     AICc model ntax
#  0.2255174 7.548965 7.945005   OUM  106


# Rates
#                  1         2
# alpha    0.0000010 0.0000010
# sigma.sq 0.1215066 0.1215066

# Optima
#                  1         2
# estimate 3.7383082 3.0147907
# se       0.1094239 0.2487064

# Arrived at a reliable solution


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
regimes <- dat[c("xMpleth", "metamorphosis",   "OU1")]
gsz <- log(dat['genomesize'])

# Plot the xMpleth hypothesis
plot(tree, regimes=regimes['xMpleth'], cex=.5)

# fit models. Results in res.bm, res.ou
res.bm <- brown(gsz,  tree)
res.ou <- apply(regimes, 2, function(x) hansen(gsz, tree, factor(x), sqrt.alpha=1, sigma=1))

##############################
#  Display results of model fits individually instead of all together as in tables below
##############################
res.ou[['xMpleth']]
res.ou[['OU1']]
res.bm

## Make nice tables of results
aic.dat <- as.data.frame(t(sapply(c(BM=res.bm, res.ou),function(x)unlist(summary(x)[c('aic','aic.c','sic', 'dof')]))))
param.names <- c('sigma.sq.matrix', 'alpha.matrix', 'theta.genomesize.M', 'theta.genomesize.Mpleth', 'theta.genomesize.x', 'theta.genomesize.global')
ou.param <- as.data.frame(t(sapply(res.ou, function(x) unlist(coef(x))[param.names])))
bm.param <- as.data.frame(t(unlist(coef(res.bm))[param.names]))
rownames(bm.param) <- "bm"
names(ou.param) <- names(bm.param) <- param.names
param <- rbind(bm.param, ou.param)

aic.dat
#                    aic    aic.c      sic dof
# BM            12.61833 12.73483 17.94521   2
# xMpleth       12.81293 13.20897 23.46669   4
# metamorphosis 13.94093 14.33697 24.59469   4
# OU1           20.08993 20.32522 28.08024   3

param
#               sigma.sq.matrix alpha.matrix theta.genomesize.M theta.genomesize.Mpleth theta.genomesize.x
# bm                  0.1323883           NA                 NA                      NA                 NA
# xMpleth             0.1272212   0.07223458                 NA               -5.965962           3.732793
# metamorphosis       0.1285745   0.07216152           3.528691                      NA          10.058045
# OU1                 0.1375711   0.05788527                 NA                      NA                 NA
#               theta.genomesize.global
# bm                                 NA
# xMpleth                            NA
# metamorphosis                      NA
# OU1                          3.700393


# Fit OUCH xMpleth model to alpha and sigma optimized by OUwie

xMpleth <- res.ou[["xMpleth"]]
xMpleth.OUwieparams <- update ( xMpleth, sqrt.alpha = sqrt(xMpleth.OUM.2$solution[1]), sigma = sqrt(xMpleth.OUM.2$solution[2]), fit=FALSE)
summary(xMpleth.OUwieparams)
xMpleth.OUM.2

### NOTE these results are very different from the OUwie results.
# $call
# hansen(data = data, tree = object, regimes = regimes, sqrt.alpha = sqrt.alpha,
#     sigma = sigma, fit = FALSE)

# $conv.code
# NULL

# $optimizer.message
# NULL

# $alpha
#       [,1]
# [1,] 1e-06

# $sigma.squared
#           [,1]
# [1,] 0.1215066

# $optima
# $optima$genomesize
#        Mpleth             x
# -6.794763e+05  3.735338e+00


# $loglik
# [1] -7.510268

# $deviance
# [1] 15.02054

# $aic
# [1] 23.02054

# $aic.c
# [1] 23.41658

# $sic
# [1] 33.67429

# $dof
# [1] 4

### OUwie results  -- both optima and IC are different!
#> xMpleth.OUM

# Fit
#       -lnL      AIC     AICc model ntax
#  0.2255174 7.548965 7.945005   OUM  106


# Rates
#                  1         2
# alpha    0.0000010 0.0000010
# sigma.sq 0.1215066 0.1215066

# Optima
#                  1         2
# estimate 3.7383082 3.0147907
# se       0.1094239 0.2487064

# Arrived at a reliable solution

# The question is why? It must be that
