require(OUwie)

############
# Trees  - one for each hypothesis, with the internal nodes labelled appropriately (below)
############
tree <- read.nexus("av_ultra_fulldataset.nex")
tree_xMpleth <- tree_xplethMpleth <- tree_xpleth <- tree_metamorphosis <- tree_mpd <- tree_mpdMpleth <- tree_MxMpleth <- tree_xplethMplethPpleth <- tree

############
# Plot nodelabels to identify which nodes to associate with which regimes
############
plot(tree)
nodelabels()
tiplabels()

############################################
# Regimes at Nodes for each hypothesis
############################################
nodedata <- read.csv("ouwie_nodelabels.csv")  # each hypothesis is nodedata$xpleth, nodedata$xMpleth, etc.
nodedata <- nodedata[,c('node','metamorphosis','mpd','mpdMpleth')]

## Node Colors for plots
## All direct developers are plethodons, but not all plethodons are direct developers
## metamorphosis
nodecol <- nodedata
nodecol[-1] <- sapply(nodecol[-1], as.character)
nodecol[nodecol=="x"] <- "black"
nodecol[nodecol=="Mpleth"] <- "red"
nodecol[nodecol=="pleth"] <- "yellow"
nodecol[nodecol=="M"] <- "green"
nodecol[nodecol=="P"] <- "blue"
nodecol[nodecol=="Ppleth"] <- "purple"

## code Regime names as numbers for the OUwie algorithm
nn <- nodedata
nn[-1] <- sapply(nn[-1], as.character)
nn[nn=="x"] <- 1
nn[nn=="Mpleth"] <- 2
nn[nn=="pleth"] <- 3
nn[nn=="M"] <- 4
nn[nn=="P"] <- 5
nn[nn=="Ppleth"] <- 6


tree_xMpleth$node.label <- nn$xMpleth
tree_xplethMpleth$node.label <- nn$xplethMpleth
tree_xpleth$node.label <- nn$xpleth
tree_metamorphosis$node.label <- nn$metamorphosis
tree_mpd$node.label <- nn$mpd
tree_mpdMpleth$node.label <- nn$mpdMpleth
tree_MxMpleth$node.label <- nn$MxMpleth
tree_xplethMplethPpleth$node.label <- nn$xplethMplethPpleth


############################################################
# Tip Data (species labels, genome size and regimes at tips)
############################################################
dat <- read.csv("tree1.dat.csv")
dat$genomesize <- log(dat$genomesize)
# make the pleth regime consistent with xpleth
#names(dat)[names(dat)=="pleth"] <- "xpleth"
#dat$xpleth <- as.character(dat$xpleth)
#dat$xpleth[dat$xpleth=="NP"] <- "x"

tdat <- dat[!is.na(dat$genomesize),]  # tip data
ndat <- dat[is.na(dat$genomesize),]   # node data

oo <- sapply(tree$tip.label, function(x) grep(x, tdat$labels))    ## reorder tdat to match order of tip labels in phylo tree

tipdat <- tdat[oo,]
tipdat <- tipdat[c("labels", "MxMpleth", "xMpleth", "xpleth", "xplethMpleth", "xplethMplethPpleth", "metamorphosis", "mpd", "mpdMpleth", "genomesize")]

## Tip Colors for plots
tipcol <- tipdat
tipcol[1:ncol(tipdat)-1] <- sapply(tipcol[1:ncol(tipdat)-1], as.character)
tipcol[tipcol=="x"] <- "black"
tipcol[tipcol=="Mpleth"] <- "red"
tipcol[tipcol=="pleth"] <- "yellow"
tipcol[tipcol=="M"] <- "green"
tipcol[tipcol=="P"] <- "blue"
tipcol[tipcol=="Ppleth"] <- "purple"


## code Regime names as numbers for the OUwie algorithm
tt <- tipdat
tt[1:ncol(tipdat)-1] <- sapply(tt[1:ncol(tipdat)-1], as.character)
tt[tt=="x"] <- 1
tt[tt=="Mpleth"] <- 2
tt[tt=="pleth"] <- 3
tt[tt=="M"] <- 4
tt[tt=="P"] <- 5
tt[tt=="Ppleth"] <- 6

### input dataframes - one for each hypothesis
data.xMpleth <- tt[c("labels", "xMpleth", "genomesize")]
data.xplethMpleth <- tt[c("labels", "xplethMpleth", "genomesize")]
data.xpleth <- tt[c("labels", "xpleth", "genomesize")]
data.metamorphosis <- tt[c("labels", "metamorphosis", "genomesize")]
data.mpd <- tt[c("labels", "mpd", "genomesize")]
data.mpdMpleth <- tt[c("labels", "mpdMpleth", "genomesize")]
data.MxMpleth <- tt[c("labels", "MxMpleth", "genomesize")]
data.xplethMplethPpleth <- tt[c("labels", "xplethMplethPpleth", "genomesize")]

###############################################################
###  xMpleth
###############################################################
plot(tree_xMpleth, cex=.6)     # plot tree with slightly reduced font size
nodelabels(pch=21, bg=nodecol$xMpleth)   # plot node labels
tiplabels(pch=21, bg=tipcol$xMpleth)   # plot tip labels
title("xMpleth")


xMpleth.OUMVA <- OUwie(tree_xMpleth,data.xMpleth,model=c("OUMVA"),root.station=TRUE)
xMpleth.OUMV <- OUwie(tree_xMpleth,data.xMpleth,model=c("OUMV"),root.station=TRUE)
xMpleth.OUM <- OUwie(tree_xMpleth,data.xMpleth,model=c("OUM"),root.station=TRUE)
xMpleth.BMS <- OUwie(tree_xMpleth,data.xMpleth,model=c("BMS"),root.station=TRUE)
xMpleth.BM1 <- OUwie(tree_xMpleth,data.xMpleth,model=c("BM1"),root.station=TRUE)

###############################################################
###  xplethMpleth
###############################################################
quartz()
plot(tree_xplethMpleth, cex=.6)     # plot tree with slightly reduced font size
nodelabels(pch=21, bg=nodecol$xplethMpleth)   # plot node lables
tiplabels(pch=21, bg=tipcol$xplethMpleth)   # plot tip labels
title("xplethMpleth")


xplethMpleth.OUMVA <- OUwie(tree_xplethMpleth,data.xplethMpleth,model=c("OUMVA"),root.station=TRUE)
xplethMpleth.OUMV <- OUwie(tree_xplethMpleth,data.xplethMpleth,model=c("OUMV"),root.station=TRUE)
xplethMpleth.OUM <- OUwie(tree_xplethMpleth,data.xplethMpleth,model=c("OUM"),root.station=TRUE)
xplethMpleth.BMS <- OUwie(tree_xplethMpleth,data.xplethMpleth,model=c("BMS"),root.station=TRUE)
xplethMpleth.BM1 <- OUwie(tree_xplethMpleth,data.xplethMpleth,model=c("BM1"),root.station=TRUE)

###############################################################
###  xpleth
###############################################################
quartz()
plot(tree_xpleth, cex=.6)     # plot tree with slightly reduced font size
nodelabels(pch=21, bg=nodecol$xpleth)   # plot node lables
tiplabels(pch=21, bg=tipcol$xpleth)   # plot tip labels
title("xpleth")


xpleth.OUMVA <- OUwie(tree_xpleth,data.xpleth,model=c("OUMVA"),root.station=TRUE)
xpleth.OUMV <- OUwie(tree_xpleth,data.xpleth,model=c("OUMV"),root.station=TRUE)
xpleth.OUM <- OUwie(tree_xpleth,data.xpleth,model=c("OUM"),root.station=TRUE)
xpleth.BMS <- OUwie(tree_xpleth,data.xpleth,model=c("BMS"),root.station=TRUE)
xpleth.BM1 <- OUwie(tree_xpleth,data.xpleth,model=c("BM1"),root.station=TRUE)

###############################################################
###  metamorphosis
###############################################################
quartz()
plot(tree_metamorphosis, cex=.6)     # plot tree with slightly reduced font size
nodelabels(pch=21, bg=nodecol$metamorphosis)   # plot node lables
tiplabels(pch=21, bg=tipcol$metamorphosis)   # plot tip labels
title("metamorphosis")


metamorphosis.OUMVA <- OUwie(tree_metamorphosis,data.metamorphosis,model=c("OUMVA"),root.station=TRUE)
metamorphosis.OUMV <- OUwie(tree_metamorphosis,data.metamorphosis,model=c("OUMV"),root.station=TRUE)
metamorphosis.OUM <- OUwie(tree_metamorphosis,data.metamorphosis,model=c("OUM"),root.station=TRUE)
metamorphosis.BMS <- OUwie(tree_metamorphosis,data.metamorphosis,model=c("BMS"),root.station=TRUE)
metamorphosis.BM1 <- OUwie(tree_metamorphosis,data.metamorphosis,model=c("BM1"),root.station=TRUE)

###############################################################
###  Make PDF's of all tree plots
###############################################################

pdf("OUwieplots%03d.pdf", width=7, height=11)

plot(tree_xMpleth, cex=.6)     # plot tree with slightly reduced font size
nodelabels(pch=21, bg=nodecol$xMpleth)   # plot node labels
tiplabels(pch=21, bg=tipcol$xMpleth)   # plot tip labels
title("xMpleth")

plot(tree_xplethMpleth, cex=.6)     # plot tree with slightly reduced font size
nodelabels(pch=21, bg=nodecol$xplethMpleth)   # plot node lables
tiplabels(pch=21, bg=tipcol$xplethMpleth)   # plot tip labels
title("xplethMpleth")

plot(tree_xpleth, cex=.6)     # plot tree with slightly reduced font size
nodelabels(pch=21, bg=nodecol$xpleth)   # plot node lables
tiplabels(pch=21, bg=tipcol$xpleth)   # plot tip labels
title("xpleth")

plot(tree_metamorphosis, cex=.6)     # plot tree with slightly reduced font size
nodelabels(pch=21, bg=nodecol$metamorphosis)   # plot node lables
tiplabels(pch=21, bg=tipcol$metamorphosis)   # plot tip labels
title("metamorphosis")

dev.off()

###############################################################
###  Gather Results to Compare Models
###############################################################

metamorphosis <- list(metamorphosis.BM1, metamorphosis.BMS, metamorphosis.OUM, metamorphosis.OUMV, metamorphosis.OUMVA)
xpleth <- list(xpleth.BM1, xpleth.BMS, xpleth.OUM, xpleth.OUMV, xpleth.OUMVA)
xMpleth <- list(xMpleth.BM1, xMpleth.BMS, xMpleth.OUM, xMpleth.OUMV, xMpleth.OUMVA)
xplethMpleth <- list(xplethMpleth.BM1, xplethMpleth.BMS, xplethMpleth.OUM, xplethMpleth.OUMV, xplethMpleth.OUMVA)

names(metamorphosis) <- names(xpleth) <- names(xMpleth) <- names(xplethMpleth) <- c("BM", "BMS", "OU", "OUS", "OUSA")

aic.meta <- sapply(metamorphosis, function(x) x$AIC)
aic.xpleth <- sapply(xpleth, function(x) x$AIC)
aic.xMpleth <- sapply(xMpleth, function(x) x$AIC)
aic.xplethMpleth <- sapply(xplethMpleth, function(x) x$AIC)

aics <- cbind(aic.meta, aic.xpleth, aic.xMpleth, aic.xplethMpleth)
write.csv(aics, "OUwie_modelfits.csv")

load("OUwie_output.Rdata")
aic.meta <- sapply(metamorphosis, function(x) x$AIC)
aic.xpleth <- sapply(xpleth, function(x) x$AIC)
aic.xMpleth <- sapply(xMpleth, function(x) x$AIC)
aic.xplethMpleth <- sapply(xplethMpleth, function(x) x$AIC)
