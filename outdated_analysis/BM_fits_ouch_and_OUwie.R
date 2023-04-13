library(ouch, quietly=TRUE)
library(OUwie, quietly=TRUE)
data(tworegime)

ouch.tree <- ape2ouch(tree)
tr <- as(ouch.tree, 'data.frame')
tr$nodes <- as.character(tr$nodes)
tr$labels <- as.character(tr$labels)
tr$ancestors <- as.character(tr$ancestors)
rownames(tr) <- tr$nodes
ouch.trait <- rep(NA, ouch.tree@nnodes)
## match tips of the ouchtree with tips of the OUwie tree
ouch.trait[sapply(as.character(trait$Genus_species), function(n) which(ouch.tree@nodelabels==n))] <- trait$X
names(ouch.trait) <- ouch.tree@nodes

ouwiefit <- OUwie(tree, trait, model="BM1", scaleHeight=TRUE, root.station=TRUE, quiet=TRUE)
ouchfit1 <- brown(ouch.trait, tree=ouch.tree)

fits <- data.frame(
  alpha=c(ouwiefit$solution["alpha",1], summary(ouchfit1)$alpha),
  sigma.sq=c(ouwiefit$solution["sigma.sq",1], summary(ouchfit1)$sigma.sq),
  loglik=c(ouwiefit$loglik,summary(ouchfit1)$loglik)
)
rownames(fits) <- c("OUwie", "ouch")
fits



