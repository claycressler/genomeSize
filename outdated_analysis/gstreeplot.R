require(ggtree)
require(treeio)
require(ggplot2)
require(dplyr)
require(cowplot)
require(aplot)

## Load in the consensus phylogenetic tree 
tree2 <- read.newick("JP_MLtree_118spp_meanbrlens.tre")
tree2$tip.label <- gsub("_", " ", tree2$tip.label)
tree <- as.treedata(tree2)  ## convert to ggtree "treedata" format 

## Load in the regime information for internal nodes
nodedata <- readRDS("nodedata_ancPlethM.RDS")
## Change regime names
nodedata <- lapply(nodedata, function(x) gsub("meta-p", "abrupt", x))
nodedata <- as.data.frame(lapply(nodedata, function(x) gsub("meta-np", "gradual", x)))
names(nodedata)[4:5] <- c("abrupt.gradual.other", "abrupt.gradual.dd.paed")

## Parse tipdata: Load in the genome size and life history data and parse into regimes
tab <- read.csv("Data_Supp_Table_9-1-21.csv")
tipdata <- mutate(tab,
       labels=paste(Genus,Species,sep=" "),
       meta.other=unname(sapply(as.character(LifeHistoryStrategy),
                                function(lh)
                                  switch(lh,"D"="other","M"="metamorphosis","Mpleth"="metamorphosis","P"="other"))),
										       meta.dd.paed=unname(sapply(as.character(LifeHistoryStrategy),
                                  function(lh) switch(lh,"D"="direct  development","M"="metamorphosis","Mpleth"="metamorphosis","P"="paedomorphosis"))),
       abrupt.gradual.other=unname(sapply(as.character(LifeHistoryStrategy),
                                    function(lh) switch(lh,"D"="other","M"="gradual","Mpleth"="abrupt","P"="other"))),
       abrupt.gradual.dd.paed=unname(sapply(as.character(LifeHistoryStrategy),
                                      function(lh) switch(lh,"D"="direct development","M"="gradual","Mpleth"="abrupt","P"="paedomorphosis"))),
       genomesize=log(GenomeSize))[-4] 
## Put the data in tip.label order
tipdata <- tipdata[match(tree2$tip.label,tipdata$labels),]
tipdata <- cbind(node=1:118, tipdata)   # add node numbers (ape order)

regimes <- rbind( tipdata,cbind(nodedata[1], Genus="", Species="", GenomeSize="", labels="", nodedata[-1], genomesize=""))  # all nodes+tips
tipdata <- tipdata[c(5,1:4,6:10)]   # reorder columns putting labels first 

## Colors used to paint life history regimes - a named vector of colors, named by life history
# "other", "metamorphosis", "abrupt", "gradual", "direct development", "paedomorphosis"
cols <- c("#66cdaa", "#ffa500", "#00ff00", "#1e90ff", "#0000ff", "#ff1493")
names(cols) <- c("other", "metamorphosis", "abrupt", "gradual", "direct development", "paedomorphosis")

# custom margins around each subplot
mar <- c(1, 1, 25, 1)

# the trees, one per regime painting
a <- ggtree(tree, aes(color=as.factor(regimes$meta.other)), size=1) +   # associate line color with regimes
     theme(plot.title = element_text(hjust = 0.5, vjust=-3, size=10),     # customize title side, justification, no legend, margins
       		legend.position="none", 
			plot.margin= margin(mar)
     		) + 
    labs(title = "meta-other")+          # character string for title
 	scale_colour_manual(values = cols[names(cols) %in% c("metamorphosis", "other")])+         # associating color with regime category
#	annotate("text", -4, 121, hjust=0, size=4, label="A") +
 	annotate("point", 1, 0,  shape=22, size=2, color=cols["metamorphosis"], fill=cols["metamorphosis"] ) +   # custom legend on bottom as annotation
	annotate("text", 13, 0, hjust=0, size=2.5, label="metamorphosis") +
 	annotate("point", 1, -2,  shape=22, size=2, color=cols["other"], fill=cols["other"] ) +
	annotate("text", 13, -2, hjust=0, size=2.5, label="other") +
	coord_cartesian(ylim=c(0,118), clip = "off")
	 
b <- ggtree(tree, aes(color=as.factor(regimes$meta.dd.paed)), size=1) + 
      theme(plot.title = element_text(hjust = 0.5, vjust=-3, size=10),
       		legend.position="none", 
			plot.margin= margin(mar)
    		) + 
    labs(title = "meta-dd-paed")+
    guides(color = guide_legend(override.aes=list(size=1.8))) + # customize legend line size
 	scale_colour_manual(values = cols[names(cols) %in% c("metamorphosis", "direct development", "paedomorphosis")])+
#	annotate("text", -4, 121, hjust=0, size=4, label="B") +
 	annotate("point", 1, 0,  shape=22, size=2, color=cols["metamorphosis"], fill=cols["metamorphosis"] ) +
	annotate("text", 13, 0, hjust=0, size=2.5, label="metamorphosis") +
 	annotate("point", 1, -2,  shape=22, size=2, color=cols["direct development"], fill=cols["direct development"] ) +
	annotate("text", 13, -2, hjust=0, size=2.5, label="direct development") +
 	annotate("point", 1, -4,  shape=22, size=2, color=cols["paedomorphosis"], fill=cols["paedomorphosis"] ) +
	annotate("text", 13, -4, hjust=0, size=2.5, label="paedomorphosis") +
	coord_cartesian(ylim=c(0,118), clip = "off")
	 
c <- ggtree(tree, aes(color=as.factor(regimes$abrupt.gradual.other)), size=1) + 
      theme(plot.title = element_text(hjust = 0.5, vjust=-3, size=9),
       		legend.position="none", 
			plot.margin= margin(mar)
     		) + 
    labs(title = expression("meta"[ab]*"-meta"[gr]*"-other"))+
    guides(color = guide_legend(override.aes=list(size=2))) + # customize legend line size
 	scale_colour_manual(values = cols[names(cols) %in% c("abrupt", "gradual", "other")])+
#	annotate("text", -4, 121, hjust=0, size=4, label="C") +
 	annotate("point", 1, 0,  shape=22, size=2, color=cols["abrupt"], fill=cols["abrupt"] ) +
	annotate("text", 13, 0, hjust=0, size=2.5, label="metamorphosis abrupt") +
 	annotate("point", 1, -2,  shape=22, size=2, color=cols["gradual"], fill=cols["gradual"] ) +
	annotate("text", 13, -2, hjust=0, size=2.5, label="metamorphosis gradual") +
 	annotate("point", 1, -4,  shape=22, size=2, color=cols["other"], fill=cols["other"] ) +
	annotate("text", 13, -4, hjust=0, size=2.5, label="other") +
	coord_cartesian(ylim=c(0,118), clip = "off")

d <- ggtree(tree, aes(color=as.factor(regimes$abrupt.gradual.dd.paed)), size=1) + 
       theme(plot.title = element_text(hjust = 0.5, vjust=-3, size=9),
       		legend.position="none", 
			plot.margin= margin(mar)
     		) + 
     labs(title = expression("meta"[ab]*"-meta"[gr]*"-dd-paed"))+
     guides(color = guide_legend(override.aes=list(size=2))) + # customize legend line size
 	 scale_colour_manual(values = cols[names(cols) %in% c("abrupt", "gradual", "direct development", "paedomorphosis")]) +
 	 geom_tiplab(size=1.5, color="black", offset=3, fontface=3)+
 	 xlim(0, 340) +
#	annotate("text", 0, 121, hjust=0, size=4, label="D") +
 	annotate("point", 1, 0,  shape=22, size=2, color=cols["abrupt"], fill=cols["abrupt"] ) +
	annotate("text", 13, 0, hjust=0, size=2.5, label="metamorphosis abrupt") +
 	annotate("point", 1, -2,  shape=22, size=2, color=cols["gradual"], fill=cols["gradual"] ) +
	annotate("text", 13, -2, hjust=0, size=2.5, label="metamorphosis gradual") +
 	annotate("point", 1, -4,  shape=22, size=2, color=cols["direct development"], fill=cols["direct development"] ) +
	annotate("text", 13, -4, hjust=0, size=2.5, label="direct development") +
 	annotate("point", 1, -6,  shape=22, size=2, color=cols["paedomorphosis"], fill=cols["paedomorphosis"] ) +
	annotate("text", 13, -6, hjust=0, size=2.5, label="paedomorphosis") +
	coord_cartesian(ylim=c(0,118), clip = "off")

bars <- ggplot(tipdata, aes(x=label,y=exp(genomesize))) + 
     geom_col(col="white",fill="blue") + coord_flip() + 
     theme(plot.margin= margin(mar + c(2,0,-5.5,0)),
     	plot.title = element_text(hjust = 0.5, vjust=-2, size=10),
      axis.line.y = element_blank(),
         axis.text.y = element_blank(),
         axis.ticks.y = element_blank(),
         axis.title.x= element_blank(),
         axis.title.y= element_blank(),
         panel.background = element_blank()
         ) +
#  	 annotate("text", 118,-4, hjust=0, size=4, label="E") +
     labs(title="Genome Size (pg)") 
	
pdf("figures/lifehistoryregimes.pdf")
multiplot(a,b,c,d,bars, ncol=5, widths=c(1,1,1,1.7,1)) 
dev.off()


# plot_grid(a,b,c,d,bars, labels="AUTO", align="h", nrow=1, rel_widths=c(1,1,1,1.5,1)) #cowplot function - not quite aligned

# d <- d %<+% tipdata  # attach dataframe to treedata
# trees <- list(meta.other=tree2,meta.dd.paed=tree2,abrupt.gradual.other=tree2,abrupt.gradual.dd.paed=tree2)		# list of trees
# class(trees) <- "treedataList"
# ggtree(trees) + facet_wrap(~.id)    # plot list of trees with facet_wrap
# plot_list(a, b, c, d, ncol=4, tag_levels="A")   #aplot function

# #get_taxa_name(d)   # get plotting order of taxa in tree useful for reordering data to match

## Check with ape plotting function
#tree_abrupt.gradual.dd.paed <- tree2     
# tree_abrupt.gradual.dd.paed$node.label <- nodedata$abrupt.gradual.dd.paed
# tree_for_plotting <- tree_abrupt.gradual.dd.paed
# tree_for_plotting$tip.label <- paste(tree_abrupt.gradual.dd.paed$tip.label, signif(tipdata$genomesize,3),sep="_")
# png(filename="ouwie_tree_with_labels.png", width=8, height=8, units='in', res=600)
# plot(tree_for_plotting, cex=0.25)
# nodelabels(pch=21, bg=as.numeric(as.factor(tree_for_plotting$node.label)))
# dev.off()

