require(ape)
require(ouch)

brldat <- read.csv("brlen_data_4-20-12.csv", header=T)
brldat$sp <- paste(brldat$genus, brldat$species, sep="_")
dat <- brldat[c("sp","genomesize","plethodontids","lifehistory")]

tree1 <- read.nexus("av_ultra_fulldataset.nex")
tree2 <- read.nexus("av_ultra_nomiss.nex")
tree3 <- read.nexus("ultra_av_fulldataset.nex")
tree4 <- read.nexus("ultra_av_nomiss.nex")

#Check:
dat$sp [ !dat$sp %in% tree1$tip.label]
tree1$tip.label[!tree1$tip.label %in% dat$sp]

dat$sp [ !dat$sp %in% tree2$tip.label]
tree2$tip.label[!tree2$tip.label %in% dat$sp]

dat$sp [ !dat$sp %in% tree3$tip.label]
tree3$tip.label[!tree3$tip.label %in% dat$sp]

dat$sp [ !dat$sp %in% tree4$tip.label]
tree4$tip.label[!tree4$tip.label %in% dat$sp]
# all OK
#####################################
#
#    TREE 1
#  tree1 <- read.nexus("av_ultra_fulldataset.nex")
#
#####################################

### Tree 1 was hand-edited and 
### data is in tree1.dat.csv
###

# paint(tree, subtree=c(1="M", "69"="D", "7"="P"), branch=c("69"="D", "7"="P", "2"="P"))
# met <- paint(tree, subtree=c("1"="M", "71"="P", "70"="M", "69"="D", "51"="M", "7"="P", "2"="P"), branch=c("69"="D", "70"="M", "71"="M", "51"="M", "7"="P", "2"="P"))
# plot(tree, regimes=met, cex=.6)

NPplethMpleth <- paint(tree, subtree=c("1"="NP", "69"="pleth", "51"="Mpleth", "10"="Mpleth"), branch=c("69"="pleth", "51"="Mpleth", "97"="Mpleth", "10"="Mpleth"))
# met <- paint(tree, subtree=c("1"="NP", "71"="P", "70"="M", "69"="D", "51"="M", "7"="P", "2"="P"), branch=c("69"="D", "70"="M", "71"="M", "51"="M", "7"="P", "2"="P"))

NPplethMpleth <- paint(tree, subtree=c("1"="NP", "68"="pleth", "51"="Mpleth", "10"="Mpleth", "7"="pleth"), branch=c("68"="pleth", "51"="Mpleth", "97"="Mpleth", "10"="Mpleth", "7"="pleth"))
plot(tree, node.names=T, regimes=NPplethMpleth, cex=.6)


dat <- read.csv("tree1.dat.csv")
rownames(dat) <- dat$nodes
tree <- ape2ouch(tree1)

pdf(file="hyp1.pdf", height=15, width=8)
plot(tree, regimes=dat['pleth'], cex=.5, node.names=T)
title("pleth")
plot(tree, regimes=dat['mpd'], cex=.5, node.names=T)
title("mpd")
plot(tree, regimes=dat['metamorphosis'], cex=.5, node.names=T)
title("metamorphosis")
plot(tree, regimes=dat['mpdMpleth'], cex=.5, node.names=T)
title("mpdMpleth")
plot(tree, regimes=dat['xMpleth'], cex=.5, node.names=T)
title("xMpleth")
plot(tree, regimes=dat['MxMpleth'], cex=.5, node.names=T)
title("MxMpleth")
dev.off()

#####################################
#
#    TREE 2
#  tree2 <- read.nexus("av_ultra_nomiss.nex")
#
#####################################
dat <- brldat
tree <- ape2ouch(tree2)
pleth <- as(tree, "data.frame")

dat <- merge(pleth, dat, by.x="labels", by.y="sp", all.x=T)
rownames(dat) <- dat$nodes
oo <- as.numeric(rownames(dat))
dat$lifehistory[ is.na(dat$lifehistory) ] <- "M"
plot(tree, regimes=dat['lifehistory'], node.names=T, cex=.65)

dat$plethodontids[ is.na(dat$plethodontids) ] <- "nonp"
plot(tree, regimes=dat['plethodontids'], node.names=T, cex=.5)
pleth <- paint(tree, subtree=c("1"="nonp", "79"="pleth"), branch=c("1"="nonp", "79"="pleth"))
quartz()
plot(tree, regimes=pleth, node.names=T, cex=.5)

quartz()
mpd <- paint(tree, subtree=c("1"="M", "15"="P", "7"="P", "79"="D", "30"="M", "42"="M", "37"="P"), branch=c("1"="M", "110"="P", "15"="P", "7"="P", "79"="M", "78"="M", "30"="M", "42"="M", "37"="P", "136"="M"))
plot(tree, regimes=mpd, node.names=T, cex=.5)

metamorphosis <- paint(tree, subtree=c("1"="M", "15"="x", "7"="x", "79"="x", "30"="M", "42"="M", "37"="x"), branch=c("1"="M", "110"="x", "15"="x", "7"="x", "79"="M", "78"="M", "30"="M", "42"="M", "37"="x", "136"="M"))
plot(tree, regimes=metamorphosis, node.names=T, cex=.5)
quartz()

mpdMpleth <- paint(tree, subtree=c("1"="M", "15"="P", "7"="P", "79"="D", "30"="Mpleth", "42"="Mpleth", "37"="P"), branch=c("1"="M", "110"="P", "15"="P", "7"="P", "79"="Mpleth", "78"="Mpleth", "30"="Mpleth", "42"="Mpleth", "37"="P", "136"="Mpleth"))
plot(tree, regimes=mpdMpleth, node.names=T, cex=.5)

quartz()
xMpleth <- paint(tree, subtree=c("1"="x", "30"="Mpleth", "42"="Mpleth", "37"="x"), branch=c("1"="x", "79"="Mpleth", "78"="Mpleth", "30"="Mpleth", "42"="Mpleth", "136"="Mpleth", "37"="x"))
plot(tree, regimes=xMpleth, node.names=T, cex=.5)

quartz()
MxMpleth <- paint(tree, subtree=c("1"="M", "15"="x", "7"="x", "79"="x", "30"="Mpleth", "42"="Mpleth", "37"="x"), branch=c("1"="M", "110"="x", "15"="x", "7"="x", "79"="Mpleth", "78"="Mpleth", "30"="Mpleth", "42"="Mpleth", "37"="x", "136"="Mpleth"))
plot(tree, regimes=MxMpleth, node.names=T, cex=.5)

rownames(dat) <- dat$nodes

regimes <- cbind(OU1="global", pleth=as.character(pleth), metamorphosis=as.character(metamorphosis), mpd=as.character(mpd), mpdMpleth=as.character(mpdMpleth), xMpleth=as.character(xMpleth), MxMpleth=as.character(MxMpleth))
regimes <- as.data.frame(apply(regimes, 2, factor))
dat <- cbind(dat, regimes[oo,])

pdf(file="hyp2.pdf", height=15, width=8)
plot(tree, regimes=dat['pleth'], cex=.5, node.names=T)
title("pleth")
plot(tree, regimes=dat['mpd'], cex=.5, node.names=T)
title("mpd")
plot(tree, regimes=dat['metamorphosis'], cex=.5, node.names=T)
title("metamorphosis")
plot(tree, regimes=dat['mpdMpleth'], cex=.5, node.names=T)
title("mpdMpleth")
plot(tree, regimes=dat['xMpleth'], cex=.5, node.names=T)
title("xMpleth")
plot(tree, regimes=dat['MxMpleth'], cex=.5, node.names=T)
title("MxMpleth")
dev.off()

write.csv(dat, file="tree2.dat.csv")
#####################################
#
#    TREE 3
#  tree3 <- read.nexus("ultra_av_fulldataset.nex")
#
#####################################
dat <- brldat
tree <- ape2ouch(tree3)
pleth <- as(tree, "data.frame")

dat <- merge(pleth, dat, by.x="labels", by.y="sp", all.x=T)
rownames(dat) <- dat$nodes
dat$lifehistory[ is.na(dat$lifehistory) ] <- "M"
plot(tree, regimes=dat['lifehistory'], node.names=T, cex=.5)

dat$plethodontids[ is.na(dat$plethodontids) ] <- "nonp"
plot(tree, regimes=dat['plethodontids'], node.names=T, cex=.5)
pleth <- paint(tree, subtree=c("1"="nonp", "89"="pleth"), branch=c("1"="nonp", "89"="pleth"))
quartz()
plot(tree, regimes=pleth, node.names=T, cex=.5)

quartz()
mpd <- paint(tree, subtree=c("1"="M", "18"="P", "8"="P", "89"="D", "34"="M", "48"="M", "42"="P"), branch=c("1"="M", "121"="P", "18"="P", "8"="P", "89"="M", "88"="M", "87"="M", "34"="M", "48"="M", "42"="P", "147"="M"))
plot(tree, regimes=mpd, node.names=T, cex=.5)

metamorphosis <- paint(tree, subtree=c("1"="M", "18"="x", "8"="x", "89"="x", "34"="M", "48"="M", "42"="x"), branch=c("1"="M", "121"="x", "18"="x", "8"="x", "89"="M", "88"="M", "87"="M", "34"="M", "48"="M", "42"="x", "147"="M"))
plot(tree, regimes=metamorphosis, node.names=T, cex=.5)
quartz()

mpdMpleth <- paint(tree, subtree=c("1"="M", "18"="P", "8"="P", "89"="D", "34"="Mpleth", "48"="Mpleth", "42"="P"), branch=c("1"="M", "121"="P", "18"="P", "8"="P", "89"="Mpleth", "88"="Mpleth", "87"="Mpleth", "34"="Mpleth", "48"="Mpleth", "42"="P", "147"="Mpleth"))
plot(tree, regimes=mpdMpleth, node.names=T, cex=.5)

quartz()
xMpleth <- paint(tree, subtree=c("1"="x", "18"="x", "8"="x", "89"="x", "34"="Mpleth", "48"="Mpleth", "42"="x"), branch=c("1"="x", "121"="x", "18"="x", "8"="x", "89"="Mpleth", "88"="Mpleth", "87"="Mpleth", "34"="Mpleth", "48"="Mpleth", "42"="x", "147"="Mpleth"))
plot(tree, regimes=xMpleth, node.names=T, cex=.5)

quartz()
MxMpleth <- paint(tree, subtree=c("1"="M", "18"="x", "8"="x", "89"="x", "34"="Mpleth", "48"="Mpleth", "42"="x"), branch=c("1"="M", "121"="x", "18"="x", "8"="x", "89"="Mpleth", "88"="Mpleth", "87"="Mpleth", "34"="Mpleth", "48"="Mpleth", "42"="x", "147"="Mpleth"))
plot(tree, regimes=MxMpleth, node.names=T, cex=.5)

regimes <- cbind(OU1="global", pleth=as.character(pleth), metamorphosis=as.character(metamorphosis), mpd=as.character(mpd), mpdMpleth=as.character(mpdMpleth), xMpleth=as.character(xMpleth), MxMpleth=as.character(MxMpleth))
regimes <- as.data.frame(apply(regimes, 2, factor))
oo <- as.numeric(rownames(dat))
dat <- cbind(dat, regimes[oo,])


pdf(file="hyp3.pdf", height=15, width=8)
plot(tree, regimes=dat['pleth'], cex=.5, node.names=T)
title("pleth")
plot(tree, regimes=dat['mpd'], cex=.5, node.names=T)
title("mpd")
plot(tree, regimes=dat['metamorphosis'], cex=.5, node.names=T)
title("metamorphosis")
plot(tree, regimes=dat['mpdMpleth'], cex=.5, node.names=T)
title("mpdMpleth")
plot(tree, regimes=dat['xMpleth'], cex=.5, node.names=T)
title("xMpleth")
plot(tree, regimes=dat['MxMpleth'], cex=.5, node.names=T)
title("MxMpleth")
dev.off()

write.csv(dat, file="tree3.dat.csv")

#####################################
#
#    TREE 4
# tree4 <- read.nexus("ultra_av_nomiss.nex")
#
#####################################
dat <- brldat
tree <- ape2ouch(tree4)
pleth <- as(tree, "data.frame")

dat <- merge(pleth, dat, by.x="labels", by.y="sp", all.x=T)
rownames(dat) <- dat$nodes

## find Plethodontid clade
dat$plethodontids[ is.na(dat$plethodontids) ] <- "nonp"
plot(tree, regimes=dat['plethodontids'], node.names=T, cex=.5)
pleth <- paint(tree, subtree=c("1"="nonp", "79"="pleth"), branch=c("1"="nonp", "79"="pleth"))
quartz()
plot(tree, regimes=pleth, node.names=T, cex=.5)

## start life history mapping
dat$lifehistory[ is.na(dat$lifehistory) ] <- "M"
plot(tree, regimes=dat['lifehistory'], node.names=T, cex=.5)

quartz()
mpd <- paint(tree, subtree=c("1"="M", "15"="P", "7"="P", "79"="D", "30"="M", "42"="M", "37"="P"), branch=c("1"="M", "110"="P", "15"="P", "7"="P", "79"="M", "78"="M", "30"="M", "42"="M", "37"="P", "136"="M"))
plot(tree, regimes=mpd, node.names=T, cex=.5)

metamorphosis <- paint(tree, subtree=c("1"="M", "15"="x", "7"="x", "79"="x", "30"="M", "42"="M", "37"="x"), branch=c("1"="M", "110"="x", "15"="x", "7"="x", "79"="M", "78"="M", "30"="M", "42"="M", "37"="x", "136"="M"))
plot(tree, regimes=metamorphosis, node.names=T, cex=.5)
quartz()

mpdMpleth <- paint(tree, subtree=c("1"="M", "15"="P", "7"="P", "79"="D", "30"="Mpleth", "42"="Mpleth", "37"="P"), branch=c("1"="M", "110"="P", "15"="P", "7"="P", "79"="Mpleth", "78"="Mpleth", "30"="Mpleth", "42"="Mpleth", "37"="P", "136"="Mpleth"))
plot(tree, regimes=mpdMpleth, node.names=T, cex=.5)

quartz()
xMpleth <- paint(tree, subtree=c("1"="x", "15"="x", "7"="x", "79"="x", "30"="Mpleth", "42"="Mpleth", "37"="x"), branch=c("1"="x", "110"="x", "15"="x", "7"="x", "79"="Mpleth", "78"="Mpleth", "30"="Mpleth", "42"="Mpleth", "37"="x", "136"="Mpleth"))
plot(tree, regimes=xMpleth, node.names=T, cex=.5)

quartz()
MxMpleth <- paint(tree, subtree=c("1"="M", "15"="x", "7"="x", "79"="x", "30"="Mpleth", "42"="Mpleth", "37"="x"), branch=c("1"="M", "110"="x", "15"="x", "7"="x", "79"="Mpleth", "78"="Mpleth", "30"="Mpleth", "42"="Mpleth", "37"="x", "136"="Mpleth"))
plot(tree, regimes=MxMpleth, node.names=T, cex=.5)

regimes <- cbind(OU1="global", pleth=as.character(pleth), metamorphosis=as.character(metamorphosis), mpd=as.character(mpd), mpdMpleth=as.character(mpdMpleth), xMpleth=as.character(xMpleth), MxMpleth=as.character(MxMpleth))
regimes <- as.data.frame(apply(regimes, 2, factor))
oo <- as.numeric(rownames(dat))
dat <- cbind(dat, regimes[oo,])


pdf(file="hyp4.pdf", height=15, width=8)
plot(tree, regimes=dat['pleth'], cex=.5, node.names=T)
title("pleth")
plot(tree, regimes=dat['mpd'], cex=.5, node.names=T)
title("mpd")
plot(tree, regimes=dat['metamorphosis'], cex=.5, node.names=T)
title("metamorphosis")
plot(tree, regimes=dat['mpdMpleth'], cex=.5, node.names=T)
title("mpdMpleth")
plot(tree, regimes=dat['xMpleth'], cex=.5, node.names=T)
title("xMpleth")
plot(tree, regimes=dat['MxMpleth'], cex=.5, node.names=T)
title("MxMpleth")
dev.off()

write.csv(dat, file="tree4.dat.csv")

