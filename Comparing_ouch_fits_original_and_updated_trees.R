## Analysis of salamander genome size evolution using the older tree and data and the newer tree and data
dat <- read.csv("tree1.dat.csv") ## old trees and regimes
rownames(dat) <- dat$nodes
gsz <- log(dat["genomesize"])
regimes <- dat[c("meta.dd.paed", "meta", "meta.nf.dd.paed","meta.nf")]
tree <- ouchtree(nodes=dat[,"nodes"],ancestors=dat[,"ancestors"], times=dat[,"times"], labels=dat[,"labels"])

bm.ouch <- brown(gsz, tree)
meta.dd.paed.ouch <- hansen(gsz, tree, regimes['meta.dd.paed'], sqrt.alpha=1, sigma=1)
meta.ouch <- hansen(gsz, tree, regimes['meta'], sqrt.alpha=1, sigma=1)
meta.nf.dd.paed.ouch <- hansen(gsz, tree, regimes['meta.nf.dd.paed'], sqrt.alpha=1, sigma=1)
meta.nf.ouch <- hansen(gsz, tree, regimes['meta.nf'], sqrt.alpha=1, sigma=1)

dat2 <- read.csv("ouchtree_ancPlethD.csv")
rownames(dat2) <- dat2$nodes
tab <- read.csv("Data_Supp_Table_9-1-21.csv")
genomesize <- rep(NA, nrow(dat2))
for (i in 1:length(genomesize)) {
  if(dat2$labels[i]=="") genomesize[i] <- NA
  else {
    genomesize[i] <- tab[which(paste(tab[,1],tab[,2],sep="_")==dat2$labels[i]),3]
  }
}
dat2$genomesize <- genomesize
gsz2 <- log(dat2["genomesize"])
regimes2 <- dat2[c("meta.other","meta.dd.paed","metap.np.other","metap.np.dd.paed")]
tree2 <- ouchtree(nodes=dat2[,"nodes"],ancestors=dat2[,"ancestors"], times=dat2[,"times"], labels=dat2[,"labels"])

bm.ouch2 <- brown(gsz2, tree2)
meta.dd.paed.ouch2 <- hansen(gsz2, tree2, regimes2['meta.dd.paed'], sqrt.alpha=1, sigma=1)
meta.ouch2 <- hansen(gsz2, tree2, regimes2['meta.other'], sqrt.alpha=1, sigma=1)
meta.nf.dd.paed.ouch2 <- hansen(gsz2, tree2, regimes2['metap.np.dd.paed'], sqrt.alpha=1, sigma=1)
meta.nf.ouch2 <- hansen(gsz2, tree2, regimes2['metap.np.other'], sqrt.alpha=1, sigma=1)

c(summary(bm.ouch)$aic.c,
  summary(meta.ouch)$aic.c,
  summary(meta.dd.paed.ouch)$aic.c,
  summary(meta.nf.ouch)$aic.c,
  summary(meta.nf.dd.paed.ouch)$aic.c)
 
c(summary(bm.ouch2)$aic.c,
  summary(meta.ouch2)$aic.c,
  summary(meta.dd.paed.ouch2)$aic.c,
  summary(meta.nf.ouch2)$aic.c,
  summary(meta.nf.dd.paed.ouch2)$aic.c)

png(file="~/Desktop/old_ouch_regimes.png", height=8, width=8, units='in', res=600)
plot(tree, regimes=regimes["meta.nf.dd.paed"], cex=0.30)
dev.off()

tree2_for_plotting <- tree2
tree2_for_plotting@nodelabels[118:235] <- paste(tree2@nodelabels[118:235], signif(gsz2[118:235,],3), sep="_")
png(file="~/Desktop/ouch_tree_with_labels.png", height=8, width=8, units='in', res=600)
plot(tree2_for_plotting, regimes=regimes2["metap.np.dd.paed"], cex=0.3)
dev.off()
