library(OUwie)
library(ouch)
library(tidyverse)

## OUwie analysis

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

bm.ouwie <- OUwie(tree_meta.other,data=tipdata[,c("labels","meta.other","genomesize")],model=c("BM1"),root.station=FALSE,scaleHeight=TRUE,quiet=TRUE)
meta.other.OUM <- OUwie(tree_meta.other,data=tipdata[,c("labels","meta.other","genomesize")],model=c("OUM"),root.station=FALSE, scaleHeight=TRUE, shift.point=1, quiet=TRUE)
#meta.other.BMS <- OUwie(tree_meta.other,data=tipdata[,c("labels","meta.other","genomesize")],model=c("BMS"),root.station=FALSE, scaleHeight=TRUE, shift.point=1, quiet=TRUE)
#meta.other.OUMA <- OUwie(tree_meta.other,data=tipdata[,c("labels","meta.other","genomesize")],model=c("OUMA"), scaleHeight=TRUE, shift.point=1, quiet=TRUE)
meta.other.OUMV <- OUwie(tree_meta.other,data=tipdata[,c("labels","meta.other","genomesize")],model=c("OUMV"), scaleHeight=TRUE, shift.point=1, quiet=TRUE)
#meta.other.OUMVA <- OUwie(tree_meta.other,data=tipdata[,c("labels","meta.other","genomesize")],model=c("OUMVA"), scaleHeight=TRUE, shift.point=1, quiet=TRUE)

meta.dd.paed.OUM <- OUwie(tree_meta.dd.paed,data=tipdata[,c("labels","meta.dd.paed","genomesize")],model=c("OUM"),root.station=FALSE, scaleHeight=TRUE, shift.point=1, quiet=TRUE)
#meta.dd.paed.BMS <- OUwie(tree_meta.dd.paed,data=tipdata[,c("labels","meta.dd.paed","genomesize")],model=c("BMS"),root.station=FALSE, scaleHeight=TRUE, shift.point=1, quiet=TRUE)
#meta.dd.paed.OUMA <- OUwie(tree_meta.dd.paed,data=tipdata[,c("labels","meta.dd.paed","genomesize")],model=c("OUMA"), scaleHeight=TRUE, shift.point=1, quiet=TRUE)
meta.dd.paed.OUMV <- OUwie(tree_meta.dd.paed,data=tipdata[,c("labels","meta.dd.paed","genomesize")],model=c("OUMV"), scaleHeight=TRUE, shift.point=1, quiet=TRUE)
#meta.dd.paed.OUMVA <- OUwie(tree_meta.dd.paed,data=tipdata[,c("labels","meta.dd.paed","genomesize")],model=c("OUMVA"), scaleHeight=TRUE, shift.point=1, quiet=TRUE)

metap.np.other.OUM <- OUwie(tree_metap.np.other,data=tipdata[,c("labels","metap.np.other","genomesize")],model=c("OUM"),root.station=FALSE, scaleHeight=TRUE, shift.point=1, quiet=TRUE)
#metap.np.other.BMS <- OUwie(tree_metap.np.other,data=tipdata[,c("labels","metap.np.other","genomesize")],model=c("BMS"),root.station=FALSE, scaleHeight=TRUE, shift.point=1, quiet=TRUE)
#metap.np.other.OUMA <- OUwie(tree_metap.np.other,data=tipdata[,c("labels","metap.np.other","genomesize")],model=c("OUMA"), scaleHeight=TRUE, shift.point=1, quiet=TRUE)
metap.np.other.OUMV <- OUwie(tree_metap.np.other,data=tipdata[,c("labels","metap.np.other","genomesize")],model=c("OUMV"), scaleHeight=TRUE, shift.point=1, quiet=TRUE)
#metap.np.other.OUMVA <- OUwie(tree_metap.np.other,data=tipdata[,c("labels","metap.np.other","genomesize")],model=c("OUMVA"), scaleHeight=TRUE, shift.point=1, quiet=TRUE)

metap.np.dd.paed.OUM <- OUwie(tree_metap.np.dd.paed,data=tipdata[,c("labels","metap.np.dd.paed","genomesize")],model=c("OUM"),root.station=FALSE, scaleHeight=TRUE, shift.point=1, quiet=TRUE)
#metap.np.dd.paed.BMS <- OUwie(tree_metap.np.dd.paed,data=tipdata[,c("labels","metap.np.dd.paed","genomesize")],model=c("BMS"),root.station=FALSE, scaleHeight=TRUE, shift.point=1, quiet=TRUE)
#metap.np.dd.paed.OUMA <- OUwie(tree_metap.np.dd.paed,data=tipdata[,c("labels","metap.np.dd.paed","genomesize")],model=c("OUMA"), scaleHeight=TRUE, shift.point=1, quiet=TRUE)
metap.np.dd.paed.OUMV <- OUwie(tree_metap.np.dd.paed,data=tipdata[,c("labels","metap.np.dd.paed","genomesize")],model=c("OUMV"), scaleHeight=TRUE, shift.point=1, quiet=TRUE)
#metap.np.dd.paed.OUMVA <- OUwie(tree_metap.np.dd.paed,data=tipdata[,c("labels","metap.np.dd.paed","genomesize")],model=c("OUMVA"), scaleHeight=TRUE, shift.point=1, quiet=TRUE)

results2 <- list(BM1=bm.ouwie,
                 meta.other.OUM=meta.other.OUM,
                 #meta.other.BMS=meta.other.BMS,
                 #meta.other.OUMA=meta.other.OUMA,
                 meta.other.OUMV=meta.other.OUMV,
                 #meta.other.OUMVA=meta.other.OUMVA,
                 metap.np.other.OUM=metap.np.other.OUM,
                 #metap.np.other.BMS=metap.np.other.BMS,
                 #metap.np.other.OUMA=metap.np.other.OUMA,
                 metap.np.other.OUMV=metap.np.other.OUMV,
                 #metap.np.other.OUMVA=metap.np.other.OUMVA,
                 meta.dd.paed.OUM=meta.dd.paed.OUM,
                 #meta.dd.paed.BMS=meta.dd.paed.BMS,
                 #meta.dd.paed.OUMA=meta.dd.paed.OUMA,
                 meta.dd.paed.OUMV=meta.dd.paed.OUMV,
                 #meta.dd.paed.OUMVA=meta.dd.paed.OUMVA,
                 metap.np.dd.paed.OUM=metap.np.dd.paed.OUM,
                 #metap.np.dd.paed.BMS=metap.np.dd.paed.BMS,
                 #metap.np.dd.paed.OUMA=metap.np.dd.paed.OUMA,
                 metap.np.dd.paed.OUMV=metap.np.dd.paed.OUMV
                 #metap.np.dd.paed.OUMVA=metap.np.dd.paed.OUMVA
)

## Note that, for some reason, there are two theta values for the BM1 model
data.frame(Hypothesis = lapply(results2, function(x) x$AICc) %>% unlist %>% names,
           loglik = lapply(results2, function(x) x$loglik) %>% unlist %>% signif(3),
           AICc = lapply(results2, function(x) x$AICc) %>% unlist %>% signif(3),
           alpha = lapply(results2, function(x) unique(x$solution['alpha',]) %>% signif(3) %>% paste(., collapse=", ")) %>% unlist,
           sigma = lapply(results2, function(x) unique(x$solution['sigma.sq',]) %>% sqrt %>% signif(3) %>% paste(., collapse=", ")) %>% unlist,
           theta = lapply(results2, function(x) x$theta[,1] %>% signif(3) %>% paste(., collapse=", ")) %>% unlist
) %>% 
  arrange(AICc) -> aics.ouwie


## ouch analysis
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
regimes2$meta.dd.paed <- as.factor(regimes2$meta.dd.paed)
regimes2$metap.np.dd.paed <- as.factor(regimes2$metap.np.dd.paed)
regimes2$meta.other <- as.factor(regimes2$meta.other)
regimes2$metap.np.other <- as.factor(regimes2$metap.np.other)


bm.ouch2 <- brown(gsz2, tree2)
meta.dd.paed.ouch2 <- hansen(gsz2, tree2, regimes2['meta.dd.paed'], sqrt.alpha=1, sigma=1)
meta.other.ouch2 <- hansen(gsz2, tree2, regimes2['meta.other'], sqrt.alpha=1, sigma=1)
metap.np.dd.paed.ouch2 <- hansen(gsz2, tree2, regimes2['metap.np.dd.paed'], sqrt.alpha=1, sigma=1)
metap.np.other.ouch2 <- hansen(gsz2, tree2, regimes2['metap.np.other'], sqrt.alpha=1, sigma=1)

data.frame(Hypothesis = c("BM1","meta.other.OUM","meta.dd.paed.OUM","metap.np.other.OUM","metap.np.dd.paed.OUM"),
           loglik = c(bm.ouch2@loglik, meta.other.ouch2@loglik, meta.dd.paed.ouch2@loglik, metap.np.other.ouch2@loglik, metap.np.dd.paed.ouch2@loglik) %>% signif(3),
           AICc = c(summary(bm.ouch2)$aic.c, summary(meta.other.ouch2)$aic.c, summary(meta.dd.paed.ouch2)$aic.c, summary(metap.np.other.ouch2)$aic.c, summary(metap.np.dd.paed.ouch2)$aic.c) %>% signif(3),
           alpha = c(NA,meta.other.ouch2@sqrt.alpha^2, meta.dd.paed.ouch2@sqrt.alpha^2, metap.np.other.ouch2@sqrt.alpha^2, metap.np.dd.paed.ouch2@sqrt.alpha^2) %>% signif(3),
           sigma = c(bm.ouch2@sigma,meta.other.ouch2@sigma, meta.dd.paed.ouch2@sigma, metap.np.other.ouch2@sigma, metap.np.dd.paed.ouch2@sigma) %>% signif(3),
           theta = c(signif(bm.ouch2@theta$genomesize,3), paste(signif(meta.other.ouch2@theta$genomesize,3), collapse=", "), paste(signif(meta.dd.paed.ouch2@theta$genomesize,3), collapse=", "), paste(signif(metap.np.other.ouch2@theta$genomesize,3), collapse=", "), paste(signif(metap.np.dd.paed.ouch2@theta$genomesize,3), collapse=", "))
) %>% arrange(AICc) -> aics.ouch

## Comparing the OUM fits between OUwie and ouch reveal that they are very similar. The difference just come down to the stationarity assumption.
## The major difference is with the BM1 fits.
aics.ouwie
aics.ouch
