fitComparison <- function(data) {
  
  #########################################################
  ## Pairwise comparision 1. mpdMpleth (OUMV) versus BM1 ##
  #########################################################
  
  ## Fit the oumv.mpdMpleth data to the mpdMpleth hypothesis under an OUMV model (the true data-generating model)
  fit.oumv.mpdMpleth.to.oumv.mpdMpleth.data <- OUwie(tree_mpdMpleth,data[["oumv.mpdMpleth.data"]],model=c("OUMV"),root.station=FALSE, scaleHeight=TRUE, shift.point=1, quiet=TRUE)
  
  ## Fit the oumv.mpdMpleth data to the BM1 hypothesis
  fit.bm1.to.oumv.mpdMpleth.data <- OUwie(tree_mpdMpleth,data[["oumv.mpdMpleth.data"]],model=c("BM1"),root.station=FALSE, scaleHeight=TRUE, shift.point=1, quiet=TRUE)
  
  ## Fit the bm1 data to the mpdMpleth hypothesis under an OUMV model
  bm1.as.mpdMpleth <- data[["bm1.data"]]
  bm1.as.mpdMpleth[,2] <- data[["oumv.mpdMpleth.data"]][,2]
  fit.oumv.mpdMpleth.to.bm1.data<- OUwie(tree_mpdMpleth,bm1.as.mpdMpleth,model=c("OUMV"),root.station=FALSE, scaleHeight=TRUE, shift.point=1, quiet=TRUE)
  
  ## Fit the bm1 data to the bm1 hypothesis (the true data-generating model)
  fit.bm1.to.bm1.data <- OUwie(tree_mpdMpleth, data[["bm1.data"]], model=c("BM1"),root.station=FALSE, scaleHeight=TRUE, shift.point=1, quiet=TRUE)
  
  #####################################################################
  ## Pairwise comparision 2. mpdMpleth (OUMV) versus mpdMpleth (OUM) ##
  #####################################################################
 
  ## Fit the mpdMpleth hypothesis (OUM) model to the oumv.mpdMpleth data
  fit.oum.mpdMpleth.to.oumv.mpdMpleth.data <- OUwie(tree_mpdMpleth,data[["oumv.mpdMpleth.data"]],model=c("OUM"),root.station=FALSE, scaleHeight=TRUE, shift.point=1, quiet=TRUE)
  
  ## Fit the mpdMpleth hypothesis (OUMV) model to the oum.mpdMpleth data
  fit.oumv.mpdMpleth.to.oum.mpdMpleth.data <- OUwie(tree_mpdMpleth,data[["oum.mpdMpleth.data"]],model=c("OUMV"),root.station=FALSE, scaleHeight=TRUE, shift.point=1, quiet=TRUE)
  
  ## Fit the mpdMpleth hypothesis (OUM) model to the oum.mpdMpleth data (the true data-generating model)
  fit.oum.mpdMpleth.to.oum.mpdMpleth.data <- OUwie(tree_mpdMpleth,data[["oum.mpdMpleth.data"]],model=c("OUM"),root.station=FALSE, scaleHeight=TRUE, shift.point=1, quiet=TRUE)
  
  ###############################################################
  ## Pairwise comparision 3. mpdMpleth (OUMV) versus mpd (OUM) ##
  ###############################################################
  
  ## Fit the mpd (OUM) model to the oumv.mpdMpleth data
  mpdMpleth.as.mpd <- data[["oumv.mpdMpleth.data"]]
  mpdMpleth.as.mpd[,2] <- data[["oum.mpd.data"]][,2]
  fit.oum.mpd.to.oumv.mpdMpleth.data <- OUwie(tree_mpd, mpdMpleth.as.mpd, model=c("OUM"),root.station=FALSE, scaleHeight=TRUE, shift.point=1, quiet=TRUE)
  
  ## Fit the mpdMpleth (OUMV) model to the oum.mpd data
  mpd.as.mpdMpleth <- data[["oum.mpd.data"]]
  mpd.as.mpdMpleth[,2] <- data[["oumv.mpdMpleth.data"]][,2]
  fit.oumv.mpdMpleth.to.oum.mpd.data <- OUwie(tree_mpdMpleth, mpd.as.mpdMpleth, model=c("OUMV"),root.station=FALSE, scaleHeight=TRUE, shift.point=1, quiet=TRUE)
  
  ## Fit the mpd (OUM) model to the oum.mpd data (the true data-generating model)
  fit.oum.mpd.to.oum.mpd.data <- OUwie(tree_mpd, data[["oum.mpd.data"]], model=c("OUM"),root.station=FALSE, scaleHeight=TRUE, shift.point=1, quiet=TRUE)
  
  #########################################################################
  ## Pairwise comparision 4. mpdMpleth (OUMV) versus metamorphosis (OUM) ##
  #########################################################################
 
  ## Fit the metamorphosis (OUM) model to the oumv.mpdMpleth data
  mpdMpleth.as.meta <- data[["oumv.mpdMpleth.data"]]
  mpdMpleth.as.meta[,2] <- data[["oum.meta.data"]][,2]
  fit.oum.meta.to.oumv.mpdMpleth.data <- OUwie(tree_metamorphosis, mpdMpleth.as.meta, model=c("OUM"),root.station=FALSE, scaleHeight=TRUE, shift.point=1, quiet=TRUE)
  
  ## Fit the mpdMpleth (OUMV) model to the oum.meta data
  meta.as.mpdMpleth <- data[["oum.meta.data"]]
  meta.as.mpdMpleth[,2] <- data[["oumv.mpdMpleth.data"]][,2]
  fit.oumv.mpdMpleth.to.oum.meta.data <- OUwie(tree_mpdMpleth, meta.as.mpdMpleth, model=c("OUMV"),root.station=FALSE, scaleHeight=TRUE, shift.point=1, quiet=TRUE)
  
  ## Fit the metamorphosis (OUM) model to the oum.meta data (the true data-generating model)
  fit.oum.meta.to.oum.meta.data <- OUwie(tree_metamorphosis, data[["oum.meta.data"]], model=c("OUM"),root.station=FALSE, scaleHeight=TRUE, shift.point=1, quiet=TRUE)
  
  #########################################################################
  
  list(fit.oumv.mpdMpleth.to.oumv.mpdMpleth.data=fit.oumv.mpdMpleth.to.oumv.mpdMpleth.data,
       fit.bm1.to.oumv.mpdMpleth.data=fit.bm1.to.oumv.mpdMpleth.data,
       fit.oumv.mpdMpleth.to.bm1.data=fit.oumv.mpdMpleth.to.bm1.data,
       fit.bm1.to.bm1.data=fit.bm1.to.bm1.data,
       fit.oum.mpdMpleth.to.oumv.mpdMpleth.data=fit.oum.mpdMpleth.to.oumv.mpdMpleth.data,
       fit.oumv.mpdMpleth.to.oum.mpdMpleth.data=fit.oumv.mpdMpleth.to.oum.mpdMpleth.data,
       fit.oum.mpdMpleth.to.oum.mpdMpleth.data=fit.oum.mpdMpleth.to.oum.mpdMpleth.data,
       fit.oum.mpd.to.oumv.mpdMpleth.data=fit.oum.mpd.to.oumv.mpdMpleth.data,
       fit.oumv.mpdMpleth.to.oum.mpd.data=fit.oumv.mpdMpleth.to.oum.mpd.data,
       fit.oum.mpd.to.oum.mpd.data=fit.oum.mpd.to.oum.mpd.data,
       fit.oum.meta.to.oumv.mpdMpleth.data=fit.oum.meta.to.oumv.mpdMpleth.data,
       fit.oumv.mpdMpleth.to.oum.meta.data=fit.oumv.mpdMpleth.to.oum.meta.data,
       fit.oum.meta.to.oum.meta.data=fit.oum.meta.to.oum.meta.data)
}

library(OUwie)
library(tidyverse)
library(parallel)
tree <- read.nexus("av_ultra_fulldataset.nex")
tree_metamorphosis <- tree_mpd <- tree_mpdMpleth <- tree
## Load the hypotheses
nodedata <- read.csv('ouwie_nodelabels.csv')
nodedata <- nodedata[,c('node','metamorphosis','mpd','mpdMpleth')]
## metamorphosis has two regimes: M, x
## mpd has three regimes: D, M, P
## mpdMpleth has four regimes: D, M, Mpleth, P
## Node Colors for plots
nodecol <- nodedata
nodecol[-1] <- sapply(nodecol[-1], as.character)
nodecol[nodecol=="x"] <- "brown"
nodecol[nodecol=="M"] <- "red"
nodecol[nodecol=="D"] <- "black"
nodecol[nodecol=="P"] <- "dodgerblue1"
nodecol[nodecol=="Mpleth"] <- "purple"
## code Regime names as numbers for the OUwie algorithm
nn <- nodedata
nn[-1] <- sapply(nn[-1], as.character)
nn[nn=="D"] <- 1
nn[nn=="M"] <- 2
nn[nn=="Mpleth"] <- 3
nn[nn=="P"] <- 4
nn[nn=="x"] <- 5
## Add node labels to trees
tree_metamorphosis$node.label <- nn$metamorphosis
tree_mpd$node.label <- nn$mpd
tree_mpdMpleth$node.label <- nn$mpdMpleth
## Tip data
dat <- read.csv("tree1.dat.csv")
dat$genomesize <- log(dat$genomesize)
tdat <- dat[!is.na(dat$genomesize),]  # tip data
ndat <- dat[is.na(dat$genomesize),]   # node data
## reorder tdat to match order of tip labels in phylo tree
oo <- sapply(tree$tip.label, function(x) grep(x, tdat$labels))
tipdat <- tdat[oo,]
tipdat <- tipdat[c("labels", "metamorphosis", "mpd", "mpdMpleth", "genomesize")]
## Tip Colors for plots
tipcol <- tipdat
tipcol[,1:(ncol(tipdat)-1)] <- sapply(tipcol[,1:(ncol(tipdat)-1)], as.character)
tipcol[tipcol=="x"] <- "brown"
tipcol[tipcol=="M"] <- "red"
tipcol[tipcol=="D"] <- "black"
tipcol[tipcol=="P"] <- "dodgerblue1"
tipcol[tipcol=="Mpleth"] <- "purple"
## code Regime names as numbers for the OUwie algorithm
tt <- tipdat
tt[,1:(ncol(tipdat)-1)] <- sapply(tt[,1:(ncol(tipdat)-1)], as.character)
tt[tt=="D"] <- 1
tt[tt=="M"] <- 2
tt[tt=="Mpleth"] <- 3
tt[tt=="P"] <- 4
tt[tt=="x"] <- 5

### input dataframes - one for each hypothesis
data.metamorphosis <- tt[c("labels", "metamorphosis", "genomesize")]
data.mpd <- tt[c("labels", "mpd", "genomesize")]
data.mpdMpleth <- tt[c("labels", "mpdMpleth", "genomesize")]

PMCdata <- readRDS("PMCdata.RDS")

PMCfits <- list(length=1000)

mclapply(PMCdata[1:100],
         fitComparison,
         mc.cores=12) -> PMCfits[1:100]



