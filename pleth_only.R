require(ouch)
dat <- read.csv("tree1.dat.pleth.csv")
rownames(dat) <- dat$nodes
dat$times <- dat$times - min(dat$times)

tree <- with(dat, ouchtree(nodes=nodes, ancestors=ancestors, times=times/max(times)))

regimes <- dat[c("mpd", "metamorphosis")]
regimes$OU1 <- "global"
gsz <- log(dat['genomesize'])
res.bm <- brown(gsz,  tree)
res.ou <- apply(regimes, 2, function(x) hansen(gsz, tree, factor(x), sqrt.alpha=1, sigma=1))

aic1.dat <- as.data.frame(t(sapply(c(BM=res.bm, res.ou),function(x)unlist(summary(x)[c('aic','aic.c','sic', 'dof')]))))

#                   aic     aic.c       sic dof
#BM            2.348372  2.506267  7.087268   2
#mpd           6.368039  7.189957 18.215278   5
#metamorphosis 4.444229  4.984770 13.922020   4
#OU1           9.697033 10.017033 16.805377   3
