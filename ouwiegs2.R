require(OUwie)

############
# Trees  - one for each hypothesis, with the internal nodes labelled appropriately (below)
#  tree is the individual tree,
#  ttree is the list of trees, one for each hypothesis
############
tree <- read.nexus("av_ultra_fulldataset.nex")

hyp <- c("metamorphosis", "mpd", "mpdMpleth", "MxMpleth", "xMpleth", "xpleth", "xplethMpleth", "xplethMplethPpleth")
ttree <- list( tree, tree, tree, tree, tree, tree, tree, tree)
names(ttree) <- hyp


############
# Plot nodelabels to identify which nodes to associate with which regimes
############
plot(tree, cex=.6)
nodelabels()
tiplabels()


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
nodecol[nodecol=="P"] <- "blue"
nodecol[nodecol=="Ppleth"] <- "purple"
nodecol[nodecol=="D"] <- "yellow"

## code Regime names as numbers for the OUwie algorithm  (it doesn't deal with factors, apparently)
nn <- nodedata
nn[-1] <- sapply(nn[-1], as.character) 
nn[nn=="x"] <- 1
nn[nn=="Mpleth"] <- 2
nn[nn=="pleth"] <- 3
nn[nn=="M"] <- 4
nn[nn=="P"] <- 5
nn[nn=="Ppleth"] <- 6
nn[nn=="D"] <- 7

## add node labels to each tree in ttree list
for (i in hyp)  { ttree[[i]]$node.label <- nn[[i]] }

############################################################
# Tip Data (species labels, genome size and regimes at tips)
############################################################
dat <- read.csv("tree1.dat.csv")
dat$genomesize <- log(dat$genomesize)

tdat <- dat[!is.na(dat$genomesize),]  # tip data
ndat <- dat[is.na(dat$genomesize),]   # node data

oo <- sapply(tree$tip.label, function(x) grep(x, tdat$labels))    ## reorder tdat to match order of tip labels in phylo tree

tipdat <- tdat[oo,]
tipdat <- tipdat[c("labels", hyp, "genomesize")]

## Tip Colors for plots
tipcol <- tipdat
tipcol[2:(ncol(tipdat)-1)] <- sapply(tipcol[2:(ncol(tipdat)-1)], as.character)
tipcol[tipcol=="x"] <- "black"
tipcol[tipcol=="Mpleth"] <- "red"
tipcol[tipcol=="pleth"] <- "yellow"
tipcol[tipcol=="M"] <- "green"
tipcol[tipcol=="P"] <- "blue"
tipcol[tipcol=="Ppleth"] <- "purple"
tipcol[tipcol=="D"] <- "yellow"


## code Regime names as numbers for the OUwie algorithm
tt <- tipdat
tt[2:(ncol(tipdat)-1)] <- sapply(tt[2:(ncol(tipdat)-1)], as.character) 
tt[tt=="x"] <- 1
tt[tt=="Mpleth"] <- 2
tt[tt=="pleth"] <- 3
tt[tt=="M"] <- 4
tt[tt=="P"] <- 5
tt[tt=="Ppleth"] <- 6
tt[tt=="D"] <- 7

### list input dataframes - one for each hypothesis
ddat <- vector("list", length=length(hyp))
ddat <- lapply( hyp, function(x) tt[c("labels", x, "genomesize")])
names(ddat) <- hyp


###############################################################
###  plottree function   -- takes tree list and the specific hypothesis as input
###############################################################
plottree <- function( treelist, hypothesis) {
	plot(treelist[[hypothesis]], cex=0.6)      						# plot tree with slightly reduced font size
	nodelabels(pch=21, bg=nodecol[[hypothesis]])	# plot node labels
	tiplabels(pch=21, bg=tipcol[[hypothesis]])		# plot tip labels
	title(hypothesis)
}

pdf("OUwieplots_allregimes.pdf", width=7, height=11)
plottree (ttree, 'metamorphosis')
plottree (ttree, 'mpd')
plottree (ttree, 'mpdMpleth')
plottree (ttree, 'MxMpleth')
plottree (ttree, 'xMpleth')
plottree (ttree, 'xpleth')
plottree (ttree, 'xplethMpleth')
plottree (ttree, 'xplethMplethPpleth')
dev.off()

###############################################################
###  Run OUwie models
###############################################################
ouwie_models <- c("BM1", "BMS", "OUM", "OUMV", "OUMVA")

metamorphosis.out <- lapply(ouwie_models, function(x) OUwie(ttree[['metamorphosis']], ddat[['metamorphosis']], model=c(x), root.station=TRUE))
mpd.out <- lapply(ouwie_models, function(x) OUwie(ttree[['mpd']], ddat[['mpd']], model=c(x), root.station=TRUE))
mpdMpleth.out <- lapply(ouwie_models, function(x) OUwie(ttree[['mpdMpleth']], ddat[['mpdMpleth']], model=c(x), root.station=TRUE))
MxMpleth.out <- lapply(ouwie_models, function(x) OUwie(ttree[['MxMpleth']], ddat[['MxMpleth']], model=c(x), root.station=TRUE))
xMpleth.out <- lapply(ouwie_models, function(x) OUwie(ttree[['xMpleth']], ddat[['xMpleth']], model=c(x), root.station=TRUE))
xpleth.out <- lapply(ouwie_models, function(x) OUwie(ttree[['xpleth']], ddat[['xpleth']], model=c(x), root.station=TRUE))
xplethMpleth.out <- lapply(ouwie_models, function(x) OUwie(ttree[['xplethMpleth']], ddat[['xplethMpleth']], model=c(x), root.station=TRUE))
xplethMplethPpleth.out <- lapply(ouwie_models, function(x) OUwie(ttree[['xplethMplethPpleth']], ddat[['xplethMplethPpleth']], model=c(x), root.station=TRUE))

###############################################################
###  Gather Results to Compare Models
###############################################################

aic.metamorphosis <- sapply(metamorphosis.out, function(x) x$AIC)
aic.mpd <- sapply(mpd.out, function(x) x$AIC)
aic.mpdMpleth <- sapply(mpdMpleth.out, function(x) x$AIC)
aic.MxMpleth <- sapply(MxMpleth.out, function(x) x$AIC)
aic.xMpleth <- sapply(xMpleth.out, function(x) x$AIC)
aic.xpleth <- sapply(xpleth.out, function(x) x$AIC)
aic.xplethMpleth <- sapply(xplethMpleth.out, function(x) x$AIC)
aic.xplethMplethPpleth <- sapply(xplethMplethPpleth.out, function(x) x$AIC)

names(aic.metamorphosis) <- names(aic.mpd) <- names(aic.mpdMpleth) <- names(aic.MxMpleth) <- names(aic.xMpleth) <- names(aic.xpleth) <- names(aic.xplethMpleth) <- names(aic.xplethMplethPpleth) <- ouwie_models

aics <- cbind(aic.metamorphosis, aic.mpd, aic.mpdMpleth, aic.MxMpleth, aic.xMpleth, aic.xpleth, aic.xplethMpleth, aic.xplethMplethPpleth)
write.csv(aics, "OUwie_modelfits2_tree1.csv")

save( metamorphosis.out, mpd.out, mpdMpleth.out, MxMpleth.out, xMpleth.out, xpleth.out, xplethMpleth.out, xplethMplethPpleth.out, file="OUwie_results_tree1.Rda" )
























