library(ape)
library(ouch)
library(OUwie)
library(tidyverse)
library(parallel)

tree <- read.nexus("av_ultra_fulldataset.nex")
tree_meta <- tree_meta.nf <- tree_meta.dd.paed <- tree_meta.nf.dd.paed <- tree
## Load the hypotheses
nodedata <- read.csv('ouwie_nodelabels.csv',stringsAsFactors=FALSE)
nodedata <- nodedata[,c('node','meta','meta.nf','meta.dd.paed','meta.nf.dd.paed')]
## meta has two regimes: M, x
## meta.dd.paed has three regimes: D, M, P
## meta.nf has three regimes, x, M, Mpleth
## meta.nf.dd.paed has four regimes: D, M, Mpleth, P
## code Regime names as numbers for the OUwie algorithm
nn <- nodedata
nn[-1] <- sapply(nn[-1], as.character)
nn[nn=="D"] <- 1
nn[nn=="M"] <- 2
nn[nn=="Mpleth"] <- 3
nn[nn=="P"] <- 4
nn[nn=="x"] <- 5
## Add node labels to trees
tree_meta$node.label <- nn$meta
tree_meta.nf$node.label <- nn$meta.nf
tree_meta.dd.paed$node.label <- nn$meta.dd.paed
tree_meta.nf.dd.paed$node.label <- nn$meta.nf.dd.paed
## Tip data
dat <- read.csv("tree1.dat.csv")
dat$genomesize <- log(dat$genomesize)
tdat <- dat[!is.na(dat$genomesize),]  # tip data
ndat <- dat[is.na(dat$genomesize),]   # node data
## reorder tdat to match order of tip labels in phylo tree
oo <- sapply(tree$tip.label, function(x) grep(x, tdat$labels))
tipdat <- tdat[oo,]
tipdat <- tipdat[c("labels", "meta", "meta.nf", "meta.dd.paed","meta.nf.dd.paed", "genomesize")]
## code Regime names as numbers for the OUwie algorithm
tt <- tipdat
tt[,1:(ncol(tipdat)-1)] <- sapply(tt[,1:(ncol(tipdat)-1)], as.character)
tt[tt=="D"] <- 1
tt[tt=="M"] <- 2
tt[tt=="Mpleth"] <- 3
tt[tt=="P"] <- 4
tt[tt=="x"] <- 5
### input dataframes - one for each hypothesis
data.meta <- tt[c("labels", "meta", "genomesize")]
data.meta.nf <- tt[c("labels", "meta.nf", "genomesize")]
data.meta.dd.paed <- tt[c("labels","meta.dd.paed","genomesize")]
data.meta.nf.dd.paed <- tt[c("labels", "meta.nf.dd.paed", "genomesize")]

fitComparison <- function(data) {
  
  #########################################################
  ## Pairwise comparison 1. meta.nf.dd.paed (OUMV) versus BM1 ##
  #########################################################
  # tIn <- Sys.time()
  
  ## Fit the metamorphosis (OUM) model to the oum.meta data (the true data-generating model)
  fit.oum.meta.to.oum.meta.data <- OUwie(tree_meta, data[["oum.meta.data"]], model=c("OUM"),root.station=FALSE, scaleHeight=TRUE, shift.point=1, quiet=TRUE)
  
  ## Fit the oum.metamorphosis data to the BM1 hypothesis
  fit.bm1.to.oum.meta.data <- OUwie(tree_meta,data[["oum.meta.data"]],model=c("BM1"),root.station=FALSE, scaleHeight=TRUE, shift.point=1, quiet=TRUE)
  
  ## Fit the bm1 data to the metamorphosis hypothesis under an OUM model
  bm1.as.meta <- data[["bm1.data"]]
  bm1.as.meta[,2] <- data[["oum.meta.data"]][,2]
  fit.oum.meta.to.bm1.data<- OUwie(tree_meta,bm1.as.meta,model=c("OUM"),root.station=FALSE, scaleHeight=TRUE, shift.point=1, quiet=TRUE)
  
  ## Fit the bm1 data to the bm1 hypothesis (the true data-generating model)
  fit.bm1.to.bm1.data <- OUwie(tree_meta, data[["bm1.data"]], model=c("BM1"),root.station=FALSE, scaleHeight=TRUE, shift.point=1, quiet=TRUE)
  
  # print(paste("Comparison #1 took", signif(Sys.time()-tIn,3)))
  
  ##################################################################
  ## Pairwise comparison 2. meta.dd.paed (OUM) versus meta (OUM) ##
  ##################################################################
  
  # tIn = Sys.time()
  
  ## Fit the metamorphosis (OUM) model to the oum.meta.dd.paed data
  meta.dd.paed.as.meta <- data[["oum.meta.dd.paed.data"]]
  meta.dd.paed.as.meta[,2] <- data[["oum.meta.data"]][,2]
  fit.oum.meta.to.oum.meta.dd.paed.data <- OUwie(tree_meta, meta.dd.paed.as.meta, model=c("OUM"),root.station=FALSE, scaleHeight=TRUE, shift.point=1, quiet=TRUE)
  
  ## Fit the meta.dd.paed (OUM) model to the oum.meta data
  meta.as.meta.dd.paed <- data[["oum.meta.data"]]
  meta.as.meta.dd.paed[,2] <- data[["oum.meta.dd.paed.data"]][,2]
  fit.oum.meta.dd.paed.to.oum.meta.data <- OUwie(tree_meta.dd.paed, meta.as.meta.dd.paed, model=c("OUM"),root.station=FALSE, scaleHeight=TRUE, shift.point=1, quiet=TRUE)
  
  ## Fit the meta.dd.paed (OUM) model to the oum.meta.dd.paed data (the true data-generating model)
  fit.oum.meta.dd.paed.to.oum.meta.dd.paed.data <- OUwie(tree_meta.dd.paed, data[["oum.meta.dd.paed.data"]], model=c("OUM"),root.station=FALSE, scaleHeight=TRUE, shift.point=1, quiet=TRUE)
  
  # print(paste("Comparison #2 took", signif(Sys.time()-tIn,3)))
  
  ##############################################################
  ## Pairwise comparison 3. meta.nf (OUMV) versus meta (OUM) ##
  ##############################################################
  # tIn <- Sys.time()
  ## Fit the metamorphosis (OUM) model to the oumv.meta.nf data
  meta.nf.as.meta <- data[["oumv.meta.nf.data"]]
  meta.nf.as.meta[,2] <- data[["oum.meta.data"]][,2]
  fit.oum.meta.to.oumv.meta.nf.data <- OUwie(tree_meta, meta.nf.as.meta, model=c("OUM"),root.station=FALSE, scaleHeight=TRUE, shift.point=1, quiet=TRUE)
  
  ## Fit the meta.nf (OUMV) model to the oum.meta data
  meta.as.meta.nf <- data[["oum.meta.data"]]
  meta.as.meta.nf[,2] <- data[["oumv.meta.nf.data"]][,2]
  fit.oumv.meta.nf.to.oum.meta.data <- OUwie(tree_meta.nf, meta.as.meta.nf, model=c("OUMV"),root.station=FALSE, scaleHeight=TRUE, shift.point=1, quiet=TRUE)
  
  ## Fit the meta.nf (OUMV) model to the oumv.meta.nf data (the true data-generating model)
  fit.oumv.meta.nf.to.oumv.meta.nf.data <- OUwie(tree_meta.nf, data[["oumv.meta.nf.data"]], model=c("OUMV"),root.station=FALSE, scaleHeight=TRUE, shift.point=1, quiet=TRUE)
  # print(paste("Comparison #3 took", signif(Sys.time()-tIn,3)))
  
  ###############################################################
  ## Pairwise comparison 4. meta.nf.dd.paed (OUMV) versus meta.dd.paed (OUM) ##
  ###############################################################
  # tIn <- Sys.time()
  ## Fit the meta.dd.paed (OUM) model to the oumv.meta.nf.dd.paed data
  meta.nf.dd.paed.as.meta.dd.paed <- data[["oumv.meta.nf.dd.paed.data"]]
  meta.nf.dd.paed.as.meta.dd.paed[,2] <- data[["oum.meta.dd.paed.data"]][,2]
  fit.oum.meta.dd.paed.to.oumv.meta.nf.dd.paed.data <- OUwie(tree_meta.dd.paed, meta.nf.dd.paed.as.meta.dd.paed, model=c("OUM"),root.station=FALSE, scaleHeight=TRUE, shift.point=1, quiet=TRUE)
  
  ## Fit the meta.nf.dd.paed (OUMV) model to the oum.meta.dd.paed data
  meta.dd.paed.as.meta.nf.dd.paed <- data[["oum.meta.dd.paed.data"]]
  meta.dd.paed.as.meta.nf.dd.paed[,2] <- data[["oumv.meta.nf.dd.paed.data"]][,2]
  fit.oumv.meta.nf.dd.paed.to.oum.meta.dd.paed.data <- OUwie(tree_meta.nf.dd.paed, meta.dd.paed.as.meta.nf.dd.paed, model=c("OUMV"),root.station=FALSE, scaleHeight=TRUE, shift.point=1, quiet=TRUE)
  
  ## Fit the meta.nf.dd.paed model to the oumv.meta.nf.dd.paed data
  fit.oumv.meta.nf.dd.paed.to.oumv.meta.nf.dd.paed.data <- OUwie(tree_meta.nf.dd.paed,data[["oumv.meta.nf.dd.paed.data"]],model=c("OUMV"),root.station=FALSE, scaleHeight=TRUE, shift.point=1, quiet=TRUE)
  # print(paste("Comparison #4 took", signif(Sys.time()-tIn,3)))
  
  ##################################################################
  ## Pairwise comparison 5. meta.nf.dd.paed (OUMV) versus meta.nf (OUMV) ##
  ##################################################################
  # tIn <- Sys.time()
  ## Fit the meta.nf.dd.paed model to the oumv.meta.nf data
  meta.nf.as.meta.nf.dd.paed <- data[["oumv.meta.nf.data"]]
  meta.nf.as.meta.nf.dd.paed[,2] <- data[["oumv.meta.nf.dd.paed.data"]][,2]
  fit.oumv.meta.nf.dd.paed.to.oumv.meta.nf.data <- OUwie(tree_meta.nf.dd.paed, meta.nf.as.meta.nf.dd.paed, model=c("OUMV"),root.station=FALSE, scaleHeight=TRUE, shift.point=1, quiet=TRUE)
  
  ## Fit the meta.nf model to the oumv.meta.nf.dd.paed data
  meta.nf.dd.paed.as.meta.nf <- data[["oumv.meta.nf.dd.paed.data"]]
  meta.nf.dd.paed.as.meta.nf[,2] <- data[["oumv.meta.nf.data"]][,2]
  fit.oumv.meta.nf.to.oumv.meta.nf.dd.paed.data <- OUwie(tree_meta.nf, meta.nf.dd.paed.as.meta.nf, model=c("OUMV"),root.station=FALSE, scaleHeight=TRUE, shift.point=1, quiet=TRUE)
  # print(paste("Comparison #5 took", signif(Sys.time()-tIn,3)))
  
  #####################################################################
  ## Pairwise comparison 6. meta.nf.dd.paed (OUMV) versus meta.nf.dd.paed (OUM) ##
  #####################################################################
  # tIn <- Sys.time()
  ## Fit the meta.nf.dd.paed hypothesis (OUM) model to the oumv.meta.nf.dd.paed data
  fit.oum.meta.nf.dd.paed.to.oumv.meta.nf.dd.paed.data <- OUwie(tree_meta.nf.dd.paed,data[["oumv.meta.nf.dd.paed.data"]],model=c("OUM"),root.station=FALSE, scaleHeight=TRUE, shift.point=1, quiet=TRUE)
  
  ## Fit the meta.nf.dd.paed hypothesis (OUMV) model to the oum.meta.nf.dd.paed data
  fit.oumv.meta.nf.dd.paed.to.oum.meta.nf.dd.paed.data <- OUwie(tree_meta.nf.dd.paed,data[["oum.meta.nf.dd.paed.data"]],model=c("OUMV"),root.station=FALSE, scaleHeight=TRUE, shift.point=1, quiet=TRUE)
  
  ## Fit the meta.nf.dd.paed hypothesis (OUM) model to the oum.meta.nf.dd.paed data (the true data-generating model)
  fit.oum.meta.nf.dd.paed.to.oum.meta.nf.dd.paed.data <- OUwie(tree_meta.nf.dd.paed,data[["oum.meta.nf.dd.paed.data"]],model=c("OUM"),root.station=FALSE, scaleHeight=TRUE, shift.point=1, quiet=TRUE)
  # print(paste("Comparison #6 took", signif(Sys.time()-tIn,3)))
  
  #########################################################################
  
  list(## Comparison #1
    fit.oum.meta.to.oum.meta.data=fit.oum.meta.to.oum.meta.data,
    fit.bm1.to.oum.meta.data=fit.bm1.to.oum.meta.data,
    fit.oum.meta.to.bm1.data=fit.oum.meta.to.bm1.data,
    fit.bm1.to.bm1.data=fit.bm1.to.bm1.data,
    ## Comparison #2
    fit.oum.meta.to.oum.meta.dd.paed.data=fit.oum.meta.to.oum.meta.dd.paed.data,
    fit.oum.meta.dd.paed.to.oum.meta.data=fit.oum.meta.dd.paed.to.oum.meta.data,
    fit.oum.meta.dd.paed.to.oum.meta.dd.paed.data=fit.oum.meta.dd.paed.to.oum.meta.dd.paed.data,
    ## Comparison #3
    fit.oum.meta.to.oumv.meta.nf.data=fit.oum.meta.to.oumv.meta.nf.data,
    fit.oumv.meta.nf.to.oum.meta.data=fit.oumv.meta.nf.to.oum.meta.data,
    fit.oumv.meta.nf.to.oumv.meta.nf.data=fit.oumv.meta.nf.to.oumv.meta.nf.data,
    ## Comparison #4
    fit.oum.meta.dd.paed.to.oumv.meta.nf.dd.paed.data=fit.oum.meta.dd.paed.to.oumv.meta.nf.dd.paed.data,
    fit.oumv.meta.nf.dd.paed.to.oum.meta.dd.paed.data=fit.oumv.meta.nf.dd.paed.to.oum.meta.dd.paed.data,
    fit.oumv.meta.nf.dd.paed.to.oumv.meta.nf.dd.paed.data=fit.oumv.meta.nf.dd.paed.to.oumv.meta.nf.dd.paed.data,
    ## Comparison #5
    fit.oumv.meta.nf.dd.paed.to.oumv.meta.nf.data=fit.oumv.meta.nf.dd.paed.to.oumv.meta.nf.data,
    fit.oumv.meta.nf.to.oumv.meta.nf.dd.paed.data=fit.oumv.meta.nf.to.oumv.meta.nf.dd.paed.data,
    ## Comparison #6
    fit.oum.meta.nf.dd.paed.to.oumv.meta.nf.dd.paed.data=fit.oum.meta.nf.dd.paed.to.oumv.meta.nf.dd.paed.data,
    fit.oumv.meta.nf.dd.paed.to.oum.meta.nf.dd.paed.data=fit.oumv.meta.nf.dd.paed.to.oum.meta.nf.dd.paed.data,
    fit.oum.meta.nf.dd.paed.to.oum.meta.nf.dd.paed.data=fit.oum.meta.nf.dd.paed.to.oum.meta.nf.dd.paed.data)
}

results <- readRDS("all_model_fits.RDS")

PMCdata <- readRDS("PMCdata.RDS")

numCores <- detectCores()
PMCfits <- vector(mode='list', length=500)
for (i in 1:31) {
  print(i)
  tIn <- Sys.time()
  mclapply(PMCdata[((i-1)*16+1):(16*i)],
           fitComparison,
           mc.cores=numCores) -> fits
  PMCfits[((i-1)*16+1):(16*i)] <- fits
  saveRDS(PMCfits, file="PMCfits.RDS")
  system("git add .")
  system("git commit -m 'updating PMCfits'")
  system("git push origin master")
  print(paste("Rep", i, "took", tIn-Sys.time()))
}

mclapply(PMCdata[497:500],
	 fitComparison,
	 mc.cores=4) -> PMCfits[497:500]
saveRDS(PMCfits, file="PMCfits.RDS")
system("git add .")
system("git commit -m 'updating PMCfits'")
system("git push origin master")

