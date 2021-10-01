## Load in the consensus phylogenetic tree 
tree2 <- read.tree("JP_MLtree_118spp_meanbrlens.tre")

## Read in the internal node regime paintings (need to match the nodes)
ouchtree.data <- read.csv('ouchtree_ancPlethD.csv',stringsAsFactors=FALSE)
## convert to an ouchtree object to make this a bit easier
with(ouchtree.data, ouchtree(nodes=nodes, ancestors=ancestors, times=times, labels=labels)) -> ouch.tree
## for each internal node 1-117 in the ouchtree, 
## - determine the names of tips descended from that node, 
## - then use findMRCA to find the node number in the phylo tree that corresponds to that node in the ouchtree
sapply(1:117, function(n) findMRCA(tree2, ouch.tree@nodelabels[seq(118,235)[which((lapply(ouch.tree@lineages, function(l) n%in%l) %>% unlist)[118:235])]])) -> reorder.nodes
nodedata <- data.frame(node=119:235,
                       meta.other=ouchtree.data$meta.other[sapply(119:235, function(n) which(reorder.nodes==n))],
                       meta.dd.paed=ouchtree.data$meta.dd.paed[sapply(119:235, function(n) which(reorder.nodes==n))],
                       metap.np.other=ouchtree.data$metap.np.other[sapply(119:235, function(n) which(reorder.nodes==n))],
                       metap.np.dd.paed=ouchtree.data$metap.np.dd.paed[sapply(119:235, function(n) which(reorder.nodes==n))])
saveRDS(nodedata, file="nodedata_ancPlethD.RDS")

## Read in the internal node regime paintings (need to match the nodes)
ouchtree.data <- read.csv('ouchtree_ancPlethM.csv',stringsAsFactors=FALSE)
## convert to an ouchtree object to make this a bit easier
with(ouchtree.data, ouchtree(nodes=nodes, ancestors=ancestors, times=times, labels=labels)) -> ouch.tree
## for each internal node 1-117, determine the names of tips descended from that node, then use findMRCA to find the node number in the phylo tree that corresponds to that node in the ouchtree
sapply(1:117, function(n) findMRCA(tree2, ouch.tree@nodelabels[seq(118,235)[which((lapply(ouch.tree@lineages, function(l) n%in%l) %>% unlist)[118:235])]])) -> reorder.nodes
nodedata2 <- data.frame(node=119:235,
                       meta.other=ouchtree.data$meta.other[sapply(119:235, function(n) which(reorder.nodes==n))],
                       meta.dd.paed=ouchtree.data$meta.dd.paed[sapply(119:235, function(n) which(reorder.nodes==n))],
                       metap.np.other=ouchtree.data$metap.np.other[sapply(119:235, function(n) which(reorder.nodes==n))],
                       metap.np.dd.paed=ouchtree.data$metap.np.dd.paed[sapply(119:235, function(n) which(reorder.nodes==n))])
saveRDS(nodedata2, file="nodedata_ancPlethM.RDS")
