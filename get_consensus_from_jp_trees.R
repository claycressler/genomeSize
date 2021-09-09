library(ape)
library(phytools)

#file = "VertLife_output.nex" #original 106 spp tree
file = "VertLife_output_updated_122_taxa_dataset.nex" #122 spp tree
t1 <- read.nexus(file)
species <- t1[[1]]$tip.label  #get spp list

t2 <- read.tree("jetz_pyron_data_doi_10.5061_dryad.cc3n6j5__v1/amph_shl_new.tre") #consensus tree
spp_wo_moldata <- setdiff(species, t2$tip.label) #spp w/o mol data are inconsistent
species <- setdiff(species, spp_wo_moldata) #get just spp w molecular data

t1.2 <-lapply(t1,keep.tip,tip=species) #all 1000 trees have only spp w mol data
dist.topo(t1.2[[1]], t1.2[[21]]) #distances are now 0 b/c all trees are from ml est
dist.topo(t1[[1]], t1[[21]]) #distances are >0 b/c some species w/o data assigned randomly

for (i in 1:100){
  for (j in 1:100){
    if (dist.topo(t1.2[[i]], t1.2[[j]]) > 0){
      print(i)
      print(j)
    }
  }
}  #all good - if using t1 (w B. robustus) there are differences

#### mean brlen
t1c<-consensus.edges(t1.2)
plot(t1c, cex = .4)
dist.topo(t1.2[[21]], t1c) #ok
#write.tree(t1c, file = "JP_MLtree_105spp_meanbrlens.tre")
write.tree(t1c, file = "JP_MLtree_122spp_meanbrlens.tre")

#### comparing to old tree
oldtree <- read.nexus("av_ultra_fulldataset.nex")
dist.topo(drop.tip(oldtree,"Batrachoseps_robustus"), t1c)
setdiff(oldtree$tip.label, species) #tip labels don't quite match
setdiff(species,oldtree$tip.label)
is.rooted(oldtree)
plot(oldtree, cex = .4)

#double checked - new names are correct
oldtree$tip.label[64] <- "Pseudoeurycea_lineola"
oldtree$tip.label[5] <- "Paradactylodon_mustersi"
oldtree$tip.label[81] <- "Batrachoseps_bramei"
oldtree<-drop.tip(oldtree, "Batrachoseps_robustus")
dist.topo(oldtree, t1c)