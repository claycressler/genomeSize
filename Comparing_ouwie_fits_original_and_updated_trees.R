## Old tree and data

## Load the tree
tree <- read.nexus("av_ultra_fulldataset.nex")
tree_meta <- tree_meta.nf <- tree_meta.dd.paed <- tree_meta.nf.dd.paed <- tree

## Load the hypotheses and assign to the interior node labels
nodedata <- read.csv('ouwie_nodelabels.csv',stringsAsFactors=FALSE)
tree_meta$node.label <- nodedata$meta
tree_meta.nf$node.label <- nodedata$meta.nf
tree_meta.dd.paed$node.label <- nodedata$meta.dd.paed
tree_meta.nf.dd.paed$node.label <- nodedata$meta.nf.dd.paed

## Load the genome data and tip regimes
dat <- read.csv("tree1.dat.csv")
dat$genomesize <- log(dat$genomesize) ## log-transform the genomesizes
tdat <- dat[!is.na(dat$genomesize),]  ## cut out all of the interior nodes
tdat <- tdat[match(tree$tip.label,tdat$labels),] ## put the tips in the same order as found in the tree

bm.ouwie <- OUwie(tree_meta,tdat[,c("labels","meta","genomesize")],model=c("BM1"),root.station=FALSE,scaleHeight=TRUE,quiet=TRUE)

meta.OUM <- OUwie(tree_meta,tdat[,c("labels","meta","genomesize")],model=c("OUM"),root.station=FALSE, scaleHeight=TRUE, shift.point=1, quiet=TRUE)
meta.BMS <- OUwie(tree_meta,tdat[,c("labels","meta","genomesize")],model=c("BMS"),root.station=FALSE, scaleHeight=TRUE, shift.point=1, quiet=TRUE)
meta.OUMA <- OUwie(tree_meta,tdat[,c("labels","meta","genomesize")],model=c("OUMA"), scaleHeight=TRUE, shift.point=1, quiet=TRUE)
meta.OUMV <- OUwie(tree_meta,tdat[,c("labels","meta","genomesize")],model=c("OUMV"), scaleHeight=TRUE, shift.point=1, quiet=TRUE)
meta.OUMVA <- OUwie(tree_meta,tdat[,c("labels","meta","genomesize")],model=c("OUMVA"), scaleHeight=TRUE, shift.point=1, quiet=TRUE)

meta.nf.OUM <- OUwie(tree_meta.nf,tdat[,c("labels","meta.nf","genomesize")],model=c("OUM"),root.station=FALSE, scaleHeight=TRUE, shift.point=1, quiet=TRUE)
meta.nf.BMS <- OUwie(tree_meta.nf,tdat[,c("labels","meta.nf","genomesize")],model=c("BMS"),root.station=FALSE, scaleHeight=TRUE, shift.point=1, quiet=TRUE)
meta.nf.OUMA <- OUwie(tree_meta.nf,tdat[,c("labels","meta.nf","genomesize")],model=c("OUMA"), scaleHeight=TRUE, shift.point=1, quiet=TRUE)
meta.nf.OUMV <- OUwie(tree_meta.nf,tdat[,c("labels","meta.nf","genomesize")],model=c("OUMV"), scaleHeight=TRUE, shift.point=1, quiet=TRUE)
meta.nf.OUMVA <- OUwie(tree_meta.nf,tdat[,c("labels","meta.nf","genomesize")],model=c("OUMVA"), scaleHeight=TRUE, shift.point=1, quiet=TRUE)

meta.dd.paed.OUM <- OUwie(tree_meta.dd.paed,tdat[,c("labels","meta.dd.paed","genomesize")],model=c("OUM"),root.station=FALSE, scaleHeight=TRUE, shift.point=1, quiet=TRUE)
meta.dd.paed.BMS <- OUwie(tree_meta.dd.paed,tdat[,c("labels","meta.dd.paed","genomesize")],model=c("BMS"),root.station=FALSE, scaleHeight=TRUE, shift.point=1, quiet=TRUE)
meta.dd.paed.OUMA <- OUwie(tree_meta.dd.paed,tdat[,c("labels","meta.dd.paed","genomesize")],model=c("OUMA"), scaleHeight=TRUE, shift.point=1, quiet=TRUE)
meta.dd.paed.OUMV <- OUwie(tree_meta.dd.paed,tdat[,c("labels","meta.dd.paed","genomesize")],model=c("OUMV"), scaleHeight=TRUE, shift.point=1, quiet=TRUE)
meta.dd.paed.OUMVA <- OUwie(tree_meta.dd.paed,tdat[,c("labels","meta.dd.paed","genomesize")],model=c("OUMVA"), scaleHeight=TRUE, shift.point=1, quiet=TRUE)

meta.nf.dd.paed.OUM <- OUwie(tree_meta.nf.dd.paed,tdat[,c("labels","meta.nf.dd.paed","genomesize")],model=c("OUM"),root.station=FALSE, scaleHeight=TRUE, shift.point=1, quiet=TRUE)
meta.nf.dd.paed.BMS <- OUwie(tree_meta.nf.dd.paed,tdat[,c("labels","meta.nf.dd.paed","genomesize")],model=c("BMS"),root.station=FALSE, scaleHeight=TRUE, shift.point=1, quiet=TRUE)
meta.nf.dd.paed.OUMA <- OUwie(tree_meta.nf.dd.paed,tdat[,c("labels","meta.nf.dd.paed","genomesize")],model=c("OUMA"), scaleHeight=TRUE, shift.point=1, quiet=TRUE)
meta.nf.dd.paed.OUMV <- OUwie(tree_meta.nf.dd.paed,tdat[,c("labels","meta.nf.dd.paed","genomesize")],model=c("OUMV"), scaleHeight=TRUE, shift.point=1, quiet=TRUE)
meta.nf.dd.paed.OUMVA <- OUwie(tree_meta.nf.dd.paed,tdat[,c("labels","meta.nf.dd.paed","genomesize")],model=c("OUMVA"), scaleHeight=TRUE, shift.point=1, quiet=TRUE)

results <- list(BM1=bm.ouwie,
                meta.OUM=meta.OUM,
                meta.BMS=meta.BMS,
                meta.OUMA=meta.OUMA,
                meta.OUMV=meta.OUMV,
                meta.OUMVA=meta.OUMVA,
                meta.nf.OUM=meta.nf.OUM,
                meta.nf.BMS=meta.nf.BMS,
                meta.nf.OUMA=meta.nf.OUMA,
                meta.nf.OUMV=meta.nf.OUMV,
                meta.nf.OUMVA=meta.nf.OUMVA,
                meta.dd.paed.OUM=meta.dd.paed.OUM,
                meta.dd.paed.BMS=meta.dd.paed.BMS,
                meta.dd.paed.OUMA=meta.dd.paed.OUMA,
                meta.dd.paed.OUMV=meta.dd.paed.OUMV,
                meta.dd.paed.OUMVA=meta.dd.paed.OUMVA,
                meta.nf.dd.paed.OUM=meta.nf.dd.paed.OUM,
                meta.nf.dd.paed.BMS=meta.nf.dd.paed.BMS,
                meta.nf.dd.paed.OUMA=meta.nf.dd.paed.OUMA,
                meta.nf.dd.paed.OUMV=meta.nf.dd.paed.OUMV,
                meta.nf.dd.paed.OUMVA=meta.nf.dd.paed.OUMVA
)

data.frame(Hypothesis = lapply(results, function(x) x$AICc) %>% unlist %>% names,
           loglik = lapply(results, function(x) x$loglik) %>% unlist %>% signif(3),
           AICc = lapply(results, function(x) x$AICc) %>% unlist %>% signif(3),
           alpha = lapply(results, function(x) unique(x$solution['alpha',]) %>% signif(3)) %>% unlist,
           sigma = lapply(results, function(x) unique(x$solution['sigma.sq',]) %>% sqrt %>% signif(3) %>% paste(., collapse=", ")) %>% unlist,
           theta = lapply(results, function(x) x$theta[,1] %>% signif(3) %>% paste(., collapse=", ")) %>% unlist
) %>% 
  arrange(AICc) -> aics
rownames(aics) <- as.character(1:nrow(aics))
aics

bm.ouwie <- OUwie(tree_meta,tdat[,c("labels","meta","genomesize")],model=c("BM1"),root.station=FALSE,scaleHeight=TRUE,quiet=TRUE)
bm.ouwie$loglik
summary(bm.ouch)$loglik

meta.OUM2 <- OUwie(tree_meta,tdat[,c("labels","meta","genomesize")],model=c("OUM"),root.station=TRUE, scaleHeight=TRUE, shift.point=1, quiet=TRUE)
meta.OUM2$loglik
summary(meta.ouch)$loglik

meta.nf.OUM2 <- OUwie(tree_meta.nf,tdat[,c("labels","meta.nf","genomesize")],model=c("OUM"),root.station=TRUE, scaleHeight=TRUE, shift.point=1, quiet=TRUE)
meta.nf.OUM2$loglik
summary(meta.nf.ouch)$loglik

meta.dd.paed.OUM2 <- OUwie(tree_meta.dd.paed,tdat[,c("labels","meta.dd.paed","genomesize")],model=c("OUM"),root.station=TRUE, scaleHeight=TRUE, shift.point=1, quiet=TRUE)
meta.dd.paed.OUM2$loglik
summary(meta.dd.paed.ouch)$loglik

meta.nf.dd.paed.OUM2 <- OUwie(tree_meta.nf.dd.paed,tdat[,c("labels","meta.nf.dd.paed","genomesize")],model=c("OUM"),root.station=TRUE, scaleHeight=TRUE, shift.point=1, quiet=TRUE)
meta.nf.dd.paed.OUM2$loglik
summary(meta.nf.dd.paed.ouch)$loglik


#######################################################################################

## New tree and data

## Load in the consensus phylogenetic tree 
tree2 <- read.tree("JP_MLtree_118spp_meanbrlens.tre")
## Load in the regime information for internal nodes
nodedata <- readRDS("nodedata_ancPlethM.RDS")
## Load in the genome size and life history data and modify it for OUwie analyses
tab <- read.csv("Data_Supp_Table_9-1-21.csv")
mutate(tab,
       labels=paste(Genus,Species,sep="_"),
       meta.other=unname(sapply(as.character(LifeHistoryStrategy),
                                function(lh)
                                  switch(lh,"D"="other","M"="metamorphosis","Mpleth"="metamorphosis","P"="other"))),
       meta.dd.paed=unname(sapply(as.character(LifeHistoryStrategy),
                                  function(lh) switch(lh,"D"="direct development","M"="metamorphosis","Mpleth"="metamorphosis","P"="paedomorphosis"))),
       metap.np.other=unname(sapply(as.character(LifeHistoryStrategy),
                                    function(lh) switch(lh,"D"="other","M"="meta-np","Mpleth"="meta-p","P"="other"))),
       metap.np.dd.paed=unname(sapply(as.character(LifeHistoryStrategy),
                                      function(lh) switch(lh,"D"="direct development","M"="meta-np","Mpleth"="meta-p","P"="paedomorphosis"))),
       genomesize=log(GenomeSize))[,5:10] -> tipdata
## Put the data in tip.label order
tipdata <- tipdata[match(tree2$tip.label,tipdata$labels),]

tree_meta.other <- tree_meta.dd.paed <- tree_metap.np.other <- tree_metap.np.dd.paed <- tree2
tree_meta.other$node.label <- nodedata$meta.other
tree_meta.dd.paed$node.label <- nodedata$meta.dd.paed
tree_metap.np.other$node.label <- nodedata$metap.np.other
tree_metap.np.dd.paed$node.label <- nodedata$metap.np.dd.paed

tree_for_plotting <- tree_metap.np.dd.paed
tree_for_plotting$tip.label <- paste(tree_metap.np.dd.paed$tip.label, signif(tipdata$genomesize,3),sep="_")
png(filename="~/Desktop/ouwie_tree_with_labels.png", width=8, height=8, units='in', res=600)
plot(tree_for_plotting, cex=0.25)
nodelabels(pch=21, bg=as.numeric(as.factor(tree_for_plotting$node.label)))
dev.off()

bm.ouwie <- OUwie(tree_meta.other,data=tipdata[,c("labels","meta.other","genomesize")],model=c("BM1"),root.station=FALSE,scaleHeight=TRUE,quiet=TRUE)

## Try a different fit method - just the tree with no node labels at all
tree_noreg <- tree2
tree_noreg$node.label <- rep("none",nrow(nodedata))
tipdata$noreg <- "none"
bm.ouwie2 <- OUwie(tree_noreg,data=tipdata[,c("labels","noreg","genomesize")],model=c("BM1"),root.station=FALSE,scaleHeight=TRUE,quiet=TRUE)

## What if you change all the regime names to numbers?
tree_noreg <- tree2
tree_noreg$node.label <- rep("1",nrow(nodedata))
tipdata$noreg <- "1"
bm.ouwie3 <- OUwie(tree_noreg,data=tipdata[,c("labels","noreg","genomesize")],model=c("BM1"),root.station=FALSE,scaleHeight=TRUE,quiet=TRUE)


## Fit all of the model variants (BMS, OUMA, OUMV, OUMVA) to each hypothesis
meta.other.OUM <- OUwie(tree_meta.other,data=tipdata[,c("labels","meta.other","genomesize")],model=c("OUM"),root.station=FALSE, scaleHeight=TRUE, shift.point=1, quiet=TRUE)
meta.other.BMS <- OUwie(tree_meta.other,data=tipdata[,c("labels","meta.other","genomesize")],model=c("BMS"),root.station=FALSE, scaleHeight=TRUE, shift.point=1, quiet=TRUE)
meta.other.OUMA <- OUwie(tree_meta.other,data=tipdata[,c("labels","meta.other","genomesize")],model=c("OUMA"), scaleHeight=TRUE, shift.point=1, quiet=TRUE)
meta.other.OUMV <- OUwie(tree_meta.other,data=tipdata[,c("labels","meta.other","genomesize")],model=c("OUMV"), scaleHeight=TRUE, shift.point=1, quiet=TRUE)
meta.other.OUMVA <- OUwie(tree_meta.other,data=tipdata[,c("labels","meta.other","genomesize")],model=c("OUMVA"), scaleHeight=TRUE, shift.point=1, quiet=TRUE)

meta.dd.paed.OUM <- OUwie(tree_meta.dd.paed,data=tipdata[,c("labels","meta.dd.paed","genomesize")],model=c("OUM"),root.station=FALSE, scaleHeight=TRUE, shift.point=1, quiet=TRUE)
meta.dd.paed.BMS <- OUwie(tree_meta.dd.paed,data=tipdata[,c("labels","meta.dd.paed","genomesize")],model=c("BMS"),root.station=FALSE, scaleHeight=TRUE, shift.point=1, quiet=TRUE)
meta.dd.paed.OUMA <- OUwie(tree_meta.dd.paed,data=tipdata[,c("labels","meta.dd.paed","genomesize")],model=c("OUMA"), scaleHeight=TRUE, shift.point=1, quiet=TRUE)
meta.dd.paed.OUMV <- OUwie(tree_meta.dd.paed,data=tipdata[,c("labels","meta.dd.paed","genomesize")],model=c("OUMV"), scaleHeight=TRUE, shift.point=1, quiet=TRUE)
meta.dd.paed.OUMVA <- OUwie(tree_meta.dd.paed,data=tipdata[,c("labels","meta.dd.paed","genomesize")],model=c("OUMVA"), scaleHeight=TRUE, shift.point=1, quiet=TRUE)

metap.np.other.OUM <- OUwie(tree_metap.np.other,data=tipdata[,c("labels","metap.np.other","genomesize")],model=c("OUM"),root.station=FALSE, scaleHeight=TRUE, shift.point=1, quiet=TRUE)
metap.np.other.BMS <- OUwie(tree_metap.np.other,data=tipdata[,c("labels","metap.np.other","genomesize")],model=c("BMS"),root.station=FALSE, scaleHeight=TRUE, shift.point=1, quiet=TRUE)
metap.np.other.OUMA <- OUwie(tree_metap.np.other,data=tipdata[,c("labels","metap.np.other","genomesize")],model=c("OUMA"), scaleHeight=TRUE, shift.point=1, quiet=TRUE)
metap.np.other.OUMV <- OUwie(tree_metap.np.other,data=tipdata[,c("labels","metap.np.other","genomesize")],model=c("OUMV"), scaleHeight=TRUE, shift.point=1, quiet=TRUE)
metap.np.other.OUMVA <- OUwie(tree_metap.np.other,data=tipdata[,c("labels","metap.np.other","genomesize")],model=c("OUMVA"), scaleHeight=TRUE, shift.point=1, quiet=TRUE)

metap.np.dd.paed.OUM <- OUwie(tree_metap.np.dd.paed,data=tipdata[,c("labels","metap.np.dd.paed","genomesize")],model=c("OUM"),root.station=FALSE, scaleHeight=TRUE, shift.point=1, quiet=TRUE)
metap.np.dd.paed.BMS <- OUwie(tree_metap.np.dd.paed,data=tipdata[,c("labels","metap.np.dd.paed","genomesize")],model=c("BMS"),root.station=FALSE, scaleHeight=TRUE, shift.point=1, quiet=TRUE)
metap.np.dd.paed.OUMA <- OUwie(tree_metap.np.dd.paed,data=tipdata[,c("labels","metap.np.dd.paed","genomesize")],model=c("OUMA"), scaleHeight=TRUE, shift.point=1, quiet=TRUE)
metap.np.dd.paed.OUMV <- OUwie(tree_metap.np.dd.paed,data=tipdata[,c("labels","metap.np.dd.paed","genomesize")],model=c("OUMV"), scaleHeight=TRUE, shift.point=1, quiet=TRUE)
metap.np.dd.paed.OUMVA <- OUwie(tree_metap.np.dd.paed,data=tipdata[,c("labels","metap.np.dd.paed","genomesize")],model=c("OUMVA"), scaleHeight=TRUE, shift.point=1, quiet=TRUE)

results2 <- list(BM1=bm.ouwie,
                meta.other.OUM=meta.other.OUM,
                meta.other.BMS=meta.other.BMS,
                meta.other.OUMA=meta.other.OUMA,
                meta.other.OUMV=meta.other.OUMV,
                meta.other.OUMVA=meta.other.OUMVA,
                metap.np.other.OUM=metap.np.other.OUM,
                metap.np.other.BMS=metap.np.other.BMS,
                metap.np.other.OUMA=metap.np.other.OUMA,
                metap.np.other.OUMV=metap.np.other.OUMV,
                metap.np.other.OUMVA=metap.np.other.OUMVA,
                meta.dd.paed.OUM=meta.dd.paed.OUM,
                meta.dd.paed.BMS=meta.dd.paed.BMS,
                meta.dd.paed.OUMA=meta.dd.paed.OUMA,
                meta.dd.paed.OUMV=meta.dd.paed.OUMV,
                meta.dd.paed.OUMVA=meta.dd.paed.OUMVA,
                metap.np.dd.paed.OUM=metap.np.dd.paed.OUM,
                metap.np.dd.paed.BMS=metap.np.dd.paed.BMS,
                metap.np.dd.paed.OUMA=metap.np.dd.paed.OUMA,
                metap.np.dd.paed.OUMV=metap.np.dd.paed.OUMV,
                metap.np.dd.paed.OUMVA=metap.np.dd.paed.OUMVA
)
data.frame(Hypothesis = lapply(results2, function(x) x$AICc) %>% unlist %>% names,
           loglik = lapply(results2, function(x) x$loglik) %>% unlist %>% signif(3),
           AICc = lapply(results2, function(x) x$AICc) %>% unlist %>% signif(3),
           alpha = lapply(results2, function(x) unique(x$solution['alpha',]) %>% signif(3) %>% paste(., collapse=", ")) %>% unlist,
           sigma = lapply(results2, function(x) unique(x$solution['sigma.sq',]) %>% sqrt %>% signif(3) %>% paste(., collapse=", ")) %>% unlist,
           theta = lapply(results2, function(x) x$theta[,1] %>% signif(3) %>% paste(., collapse=", ")) %>% unlist
) %>% 
  arrange(AICc) -> aics2
rownames(aics2) <- as.character(1:nrow(aics2))
aics2


## Idiot check: are the names, regimes, and genome sizes being assigned properly? 
## Yes they are.
cbind(meta.other.OUM$data[sapply(paste(tab[,1],tab[,2],sep="_"), function(n) which(rownames(meta.other.OUM$data)==n)),1],
      unname(sapply(as.character(tab[,4]), function(n) switch(n, "M"="1","Mpleth"="1","D"="2","P"="2")))) %>%
  apply(., 1, function(x) x[1]==x[2])
cbind(meta.other.OUM$data[sapply(paste(tab[,1],tab[,2],sep="_"), function(n) which(rownames(meta.other.OUM$data)==n)),2],
      log(tab[,3])) %>% apply(., 1, function(x) x[1]==x[2])

cbind(metap.np.other.OUM$data[sapply(paste(tab[,1],tab[,2],sep="_"), function(n) which(rownames(metap.np.other.OUM$data)==n)),1],
      unname(sapply(as.character(tab[,4]), function(n) switch(n, "M"="1","Mpleth"="2","D"="3","P"="3")))) %>%
  apply(., 1, function(x) x[1]==x[2])
cbind(metap.np.other.OUM$data[sapply(paste(tab[,1],tab[,2],sep="_"), function(n) which(rownames(metap.np.other.OUM$data)==n)),2],
      log(tab[,3])) %>% apply(., 1, function(x) x[1]==x[2])

cbind(metap.np.dd.paed.OUM$data[sapply(paste(tab[,1],tab[,2],sep="_"), function(n) which(rownames(metap.np.dd.paed.OUM$data)==n)),1],
      unname(sapply(as.character(tab[,4]), function(n) switch(n, "M"="2","Mpleth"="3","D"="1","P"="4")))) %>%
  apply(., 1, function(x) x[1]==x[2])
cbind(metap.np.dd.paed.OUM$data[sapply(paste(tab[,1],tab[,2],sep="_"), function(n) which(rownames(metap.np.dd.paed.OUM$data)==n)),2],
      log(tab[,3])) %>% apply(., 1, function(x) x[1]==x[2])

cbind(meta.dd.paed.OUM$data[sapply(paste(tab[,1],tab[,2],sep="_"), function(n) which(rownames(meta.dd.paed.OUM$data)==n)),1],
      unname(sapply(as.character(tab[,4]), function(n) switch(n, "M"="2","Mpleth"="2","D"="1","P"="3")))) %>%
  apply(., 1, function(x) x[1]==x[2])
cbind(meta.dd.paed.OUM$data[sapply(paste(tab[,1],tab[,2],sep="_"), function(n) which(rownames(meta.dd.paed.OUM$data)==n)),2],
      log(tab[,3])) %>% apply(., 1, function(x) x[1]==x[2])


## compare against ouch fit of same hypothesis by getting OUwie to assume stationarity of the root
bm.ouwie$loglik
summary(bm.ouch2)$loglik

meta.other.OUM2 <- OUwie(tree_meta.other,data=tipdata[,c("labels","meta.other","genomesize")],model=c("OUM"),root.station=TRUE, scaleHeight=TRUE, shift.point=1, quiet=TRUE)
meta.other.OUM2$loglik
summary(meta.ouch2)$loglik

meta.dd.paed.OUM2 <- OUwie(tree_meta.dd.paed,data=tipdata[,c("labels","meta.dd.paed","genomesize")],model=c("OUM"),root.station=TRUE, scaleHeight=TRUE, shift.point=1, quiet=TRUE)
meta.dd.paed.OUM2$loglik
summary(meta.dd.paed.ouch2)$loglik

metap.np.other.OUM2 <- OUwie(tree_metap.np.other,data=tipdata[,c("labels","metap.np.other","genomesize")],model=c("OUM"),root.station=TRUE, scaleHeight=TRUE, shift.point=1, quiet=TRUE)
metap.np.other.OUM2$loglik
summary(meta.nf.ouch2)$loglik

metap.np.dd.paed.OUM2 <- OUwie(tree_metap.np.dd.paed,data=tipdata[,c("labels","metap.np.dd.paed","genomesize")],model=c("OUM"),root.station=TRUE, scaleHeight=TRUE, shift.point=1, quiet=TRUE)
metap.np.dd.paed.OUM2$loglik
summary(meta.nf.dd.paed.ouch2)$loglik
