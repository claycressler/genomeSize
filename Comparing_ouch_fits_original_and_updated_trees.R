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

dat2 <- read.csv("ouchtree_ancPlethM.csv")
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

plot(tree, regimes=regimes['meta.nf.dd.paed'])

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

## Check to make sure data and regimes are assigned properly
cbind(exp(meta.ouch2@data$genomesize), dat2$genomesize) %>% apply(., 1, function(m) m[1]-m[2]) %>% max(.,na.rm=T)
cbind(meta.ouch2@regimes$meta.other, dat2$meta.other) %>% apply(., 1, function(m) m[1]==m[2]) %>% all





## Prune the ape tree to include only those species found in the original tree, convert to an ouchtree, and rerun the analyses
## Load the tree
tree <- read.nexus("av_ultra_fulldataset.nex")

## Load in the consensus phylogenetic tree 
tree2 <- read.tree("JP_MLtree_118spp_meanbrlens.tre")

to.include <- names(which(sapply(tree2$tip.label, function(t) t%in%tree$tip.label)))
tree3 <- keep.tip(tree2, to.include)
tree4 = ape2ouch(tree3)
## Paint regimes onto the tree
paint(tree4, subtree=c("1"="meta-np","8"="paedomorphosis","21"="paedomorphosis","91"="meta-p","45"="direct development", "48"="paedomorphosis", "36"="meta-p", "87"="direct development"), branch=c("1"="meta-np","36"="meta-p", "8"="paedomorphosis", "21"="paedomorphosis", "121"="paedomorphosis", "48"="paedomorphosis", "87"="direct development")) -> reg
png(filename="~/Downloads/Tree4.png",height=8, width=11, units='in', res=500)
plot(tree4,regimes=reg,node.names=TRUE,cex=0.4)
dev.off()

regimes4 <- data.frame(meta.other=unname(sapply(as.character(reg),function(r) switch(r, "direct development"="other", "paedomorphosis"="other", "meta-np"="meta","meta-p"="meta"))),
                       meta.dd.paed=unname(sapply(as.character(reg),function(r) switch(r, "direct development"="direct development", "paedomorphosis"="paedomorphosis","meta-np"="meta","meta-p"="meta"))),
                       metap.np.other=unname(sapply(as.character(reg),function(r) switch(r, "direct development"="other", "paedomorphosis"="other", "meta-np"="meta-np","meta-p"="meta-p"))),
                       metap.np.dd.paed=as.character(reg))
tab <- read.csv("Data_Supp_Table_9-1-21.csv")
genomesize <- rep(NA, nrow(regimes4))
for (i in 1:length(genomesize)) {
  if(tree4@nodelabels[i]=="") genomesize[i] <- NA
  else {
    genomesize[i] <- log(tab[which(paste(tab[,1],tab[,2],sep="_")==tree4@nodelabels[i]),3])
  }
}
regimes4$genomesize <- genomesize
gsz4 <- regimes4["genomesize"]

bm.ouch4 <- brown(gsz4, tree4)
meta.dd.paed.ouch4 <- hansen(gsz4, tree4, regimes4['meta.dd.paed'], sqrt.alpha=1, sigma=1)
meta.ouch4 <- hansen(gsz4, tree4, regimes4['meta.other'], sqrt.alpha=1, sigma=1)
meta.nf.dd.paed.ouch4 <- hansen(gsz4, tree4, regimes4['metap.np.dd.paed'], sqrt.alpha=1, sigma=1)
meta.nf.ouch4 <- hansen(gsz4, tree4, regimes4['metap.np.other'], sqrt.alpha=1, sigma=1)

c(summary(bm.ouch4)$aic.c,
  summary(meta.ouch4)$aic.c,
  summary(meta.dd.paed.ouch4)$aic.c,
  summary(meta.nf.ouch4)$aic.c,
  summary(meta.nf.dd.paed.ouch4)$aic.c)


