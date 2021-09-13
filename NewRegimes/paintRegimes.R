require(ape)
require(ouch)

dat <- read.csv("Data_Supp_Table_corrected_all_regimes.csv")
ouwie.dat <- read.csv("ouwie_nodelabels.csv")

tree4a <- read.nexus("Mplethodontid_Mnon-plethodontid_DD_Paed_ancestral_plethodontid_directdevelops.tre")
tree <- ape2ouch(tree4a)

plot(tree, text_opts = list(cex=.3, offset=0.1), margin=c(.3), lwd=3, node.names=T)

treea <- read.tree("JP_MLtree_118spp_meanbrlens.tre")
plot(treea, cex=.4, no.margin=T)
tree <- ape2ouch(treea)

treedat <- as(tree, "data.frame")
treedat <- merge(treedat, dat, by="labels", all=T)
rownames(treedat) <- treedat$nodes

regimes <- treedat[c("meta.nf.dd.paed", "meta.nf", "meta.dd.paed", "meta")]
regimes[is.na(regimes)] <- "M"

regimes2 <- regimes
regimes2[c("meta", "meta.dd.paed")][regimes2[c("meta", "meta.dd.paed")]=="M"] <- "metamorphosis"  #2
regimes2[regimes2=="D"] <- "direct development"  #1
regimes2[regimes2=="M"] <- "meta-np"   # mf   #3
regimes2[regimes2=="Mpleth"] <- "meta-p"  #mnf  #4
regimes2[regimes2=="P"] <- "paedomorphosis"  #5
regimes2[regimes2=="x"] <- "other"  #6

# Again, Mplethodontid = Meta_nf
# And Mnon-plethodontid = Meta_f

reg <- palette <- c("#000000", "chartreuse",  "#E69F00","#009E73", "#0072B2", "#CC79A7")
names(reg) <- c("direct development", "metamorphosis", "meta-np", "meta-p", "paedomorphosis", "other" )

mypal <- function(x) {return(reg[names(reg)])}

##### Painting new regimes 

# r4 = "metap.np.dd.paed"
  r4 <- paint(tree, subtree=c("1"="meta-np", 
  											"9"="paedomorphosis",
  											"22"="paedomorphosis",
  											"114"="meta-np",
  											"112"="direct development",
  											"37"="meta-p",
  											"111"="meta-p",
  											"49"="paedomorphosis",
  											"109"="direct development"
 											), 
  							branch=c("1"="meta-np",
  							                "9"="paedomorphosis",
  											"22"="paedomorphosis",
  											"114"="meta-np",
   											"144"="paedomorphosis",
  											"112"="direct development",
 											"37"="meta-p",
  											"111"="meta-p",
  											"49"="paedomorphosis",
  											"109"="direct development"
 											))
# r3 = "metap.np.other"
  r3 <- paint(tree, subtree=c("1"="meta-np", 
  											"9"="other",
  											"22"="other",
  											"114"="meta-np",
  											"112"="other",
  											"37"="meta-p",
  											"111"="meta-p",
  											"49"="other",
  											"109"="other"
 											), 
  							branch=c("1"="meta-np",
  							                "9"="other",
  											"22"="other",
  											"114"="meta-np",
   											"144"="other",
  											"112"="other",
 											"37"="meta-p",
  											"111"="meta-p",
  											"49"="other",
  											"109"="other"
 											))

# r3a = "meta.dd.paed"
  r3a <- paint(tree, subtree=c("1"="metamorphosis", 
  											"9"="paedomorphosis",
  											"22"="paedomorphosis",
  											"114"="metamorphosis",
  											"112"="direct development",
  											"37"="metamorphosis",
  											"111"="metamorphosis",
  											"49"="paedomorphosis",
  											"109"="direct development"
 											), 
  							branch=c("1"= "metamorphosis",
  							                "9"="paedomorphosis",
  											"22"="paedomorphosis",
  											"114"="metamorphosis",
   											"144"="paedomorphosis",
  											"112"="direct development",
 											"37"="metamorphosis",
  											"111"="metamorphosis",
  											"49"="paedomorphosis",
  											"109"="direct development"
 											))
# r2 = "meta.other"
  r2 <- paint(tree, subtree=c("1"="metamorphosis", 
  											"9"="other",
  											"22"="other",
  											"114"="metamorphosis",
  											"112"="other",
  											"37"="metamorphosis",
  											"111"="metamorphosis",
  											"49"="other",
  											"109"="other"
 											), 
  							branch=c("1"= "metamorphosis",
  							                "9"="other",
  											"22"="other",
  											"114"="metamorphosis",
   											"144"="other",
  											"112"="other",
 											"37"="metamorphosis",
  											"111"="metamorphosis",
  											"49"="other",
  											"109"="other"
 											))

pdf(file="Regimes_ancPlethD.pdf")
plot(tree, text_opts = list(cex=.3, offset=0.1), margin=c(.3), regimes= r4, palette=mypal, lwd=3, main="metap.np.dd.paed")
plot(tree, text_opts = list(cex=.3, offset=0.1), margin=c(.3), regimes= r3, palette=mypal, lwd=3, main="metap.np.other")
plot(tree, text_opts = list(cex=.3, offset=0.1), margin=c(.3), regimes= r3a, palette=mypal, lwd=3, main="meta.dd.paed")
plot(tree, text_opts = list(cex=.3, offset=0.1), margin=c(.3), regimes= r2, palette=mypal, lwd=3, main="meta.other")
dev.off()

regimes <- data.frame("meta.other"=r2, "meta.dd.paed"=r3a, "metap.np.other"=r3, "metap.np.dd.paed"=r4)
rownames(regimes) <- names(r4)
regimes$nodes <- names(r4)
treefile <- as(tree,"data.frame")
treefile <- merge(treefile, regimes, by="nodes")
oo <- order(as.numeric(treefile$nodes))
write.csv(treefile[oo,], file="ouchtree_ancPlethD.csv", row.names=F)

############  as above but assuming ancestor to Plethodons (node 112) = M or Mp

# r4 = "metap.np.dd.paed"
  r4 <- paint(tree, subtree=c("1"="meta-np", 
  											"9"="paedomorphosis",
  											"22"="paedomorphosis",
  											"114"="meta-np",
  											"46"="direct development",
   											"37"="meta-p",
  											"111"="meta-p",
  											"49"="paedomorphosis",
  											"109"="direct development"
 											), 
  							branch=c("1"="meta-np",
  							                "9"="paedomorphosis",
  											"22"="paedomorphosis",
  											"114"="meta-np",
   											"144"="paedomorphosis",
  											"46"="direct development",
  											"112"="meta-p",
 											"37"="meta-p",
  											"111"="meta-p",
  											"49"="paedomorphosis",
  											"109"="direct development"
 											))
# r3 = "metap.np.other"
  r3 <- paint(tree, subtree=c("1"="meta-np", 
  											"9"="other",
  											"22"="other",
  											"114"="meta-np",
  											"46"="other",
  											"37"="meta-p",
  											"111"="meta-p",
  											"49"="other",
  											"109"="other"
 											), 
  							branch=c("1"="meta-np",
  										    "9"="other",
  											"22"="other",
  											"114"="meta-np",
   											"144"="other",
  											"46"="other",
  											"112"="meta-p",
 											"37"="meta-p",
  											"111"="meta-p",
  											"49"="other",
  											"109"="other"
 											))

# r3a = "meta.dd.paed"
  r3a <- paint(tree, subtree=c("1"="metamorphosis", 
  											"9"="paedomorphosis",
  											"22"="paedomorphosis",
  											"114"="metamorphosis",
  											"46"="direct development",
  											"37"="metamorphosis",
  											"111"="metamorphosis",
  											"49"="paedomorphosis",
  											"109"="direct development"
 											), 
  							branch=c("1"="metamorphosis",
  							                "9"="paedomorphosis",
  											"22"="paedomorphosis",
  											"114"="metamorphosis",
   											"144"="paedomorphosis",
  											"46"="direct development",
  											"112"="metamorphosis",
 											"37"="metamorphosis",
  											"111"="metamorphosis",
  											"49"="paedomorphosis",
  											"109"="direct development"
 											))
# r2 = "meta.other"
  r2 <- paint(tree, subtree=c("1"="metamorphosis", 
  											"9"="other",
  											"22"="other",
  											"114"="metamorphosis",
  											"46"="other",
  											"37"="metamorphosis",
  											"111"="metamorphosis",
  											"49"="other",
  											"109"="other"
 											), 
  							branch=c("1"="metamorphosis",
  										    "9"="other",
  											"22"="other",
  											"114"="metamorphosis",
   											"144"="other",
  											"46"="other",
 											"37"="metamorphosis",
  											"111"="metamorphosis",
  											"49"="other",
  											"109"="other"
 											))

pdf(file="Regimes_ancPlethM.pdf")
plot(tree, text_opts = list(cex=.3, offset=0.1), margin=c(.3), regimes= r4, palette=mypal, lwd=3, main="metap.np.dd.paed")
plot(tree, text_opts = list(cex=.3, offset=0.1), margin=c(.3), regimes= r3, palette=mypal, lwd=3, main="metap.np.other")
plot(tree, text_opts = list(cex=.3, offset=0.1), margin=c(.3), regimes= r3a, palette=mypal, lwd=3, main="meta.dd.paed")
plot(tree, text_opts = list(cex=.3, offset=0.1), margin=c(.3), regimes= r2, palette=mypal, lwd=3, main="meta.other")
dev.off()

regimes <- data.frame("meta.other"=r2, "meta.dd.paed"=r3a, "metap.np.other"=r3, "metap.np.dd.paed"=r4)
rownames(regimes) <- names(r4)
regimes$nodes <- names(r4)
treefile <- as(tree,"data.frame")
treefile <- merge(treefile, regimes, by="nodes")
oo <- order(as.numeric(treefile$nodes))
write.csv(treefile[oo,], file="ouchtree_ancPlethM.csv", row.names=F)





