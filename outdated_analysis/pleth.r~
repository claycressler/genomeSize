require(ape)
require(ouch)

dat1 <- read.csv("tree1.dat.csv")
dat2 <- read.csv("tree2.dat.csv")
dat3 <- read.csv("tree3.dat.csv")
dat4 <- read.csv("tree4.dat.csv")

tree1 <- read.nexus("av_ultra_fulldataset.nex")
tree2 <- read.nexus("av_ultra_nomiss.nex")
tree3 <- read.nexus("ultra_av_fulldataset.nex")
tree4 <- read.nexus("ultra_av_nomiss.nex")

####
#### Analysis Tree 1
####
tree <- ape2ouch(tree1)
dat <- dat1
rownames(dat) <- dat$nodes
regimes <- dat[c("metamorphosis", "mpd", "mpdMpleth", "mpdND", "MxMpleth", "xMpleth", "xpleth", "xplethMpleth", "xplethMplethPpleth", "OU1")]
gsz <- log(dat['genomesize'])
res.bm <- brown(gsz,  tree)
res.ou <- apply(regimes, 2, function(x) hansen(gsz, tree, factor(x), sqrt.alpha=1, sigma=1))

aic1.dat <- as.data.frame(t(sapply(c(BM=res.bm, res.ou),function(x)unlist(summary(x)[c('aic','aic.c','sic', 'dof')]))))

tree1 <- tree
regimes1 <- regimes
gsz1 <- gsz
res.bm1 <- res.bm
res.ou1 <- res.ou

####
#### Analysis Tree 2
####
tree <- ape2ouch(tree2)
dat <- dat2
rownames(dat) <- dat$nodes
regimes <- dat[c("mpdMpleth", "xMpleth", "MxMpleth", "mpd", "metamorphosis")]
gsz <- log(dat['genomesize'])
res.bm <- brown(gsz,  tree)
res.ou <- apply(regimes, 2, function(x) hansen(gsz, tree, factor(x), sqrt.alpha=1, sigma=1))

aic2.dat <- as.data.frame(t(sapply(c(BM=res.bm, res.ou),function(x)unlist(summary(x)[c('aic','aic.c','sic', 'dof')]))))

tree2 <- tree
regimes2 <- regimes
gsz2 <- gsz
res.bm2 <- res.bm
res.ou2 <- res.ou

####
#### Analysis Tree 3
####
tree <- ape2ouch(tree3)
dat <- dat3
rownames(dat) <- dat$nodes
regimes <- dat[c("mpdMpleth", "xMpleth", "MxMpleth", "mpd", "metamorphosis")]
gsz <- log(dat['genomesize'])
res.bm <- brown(gsz,  tree)
res.ou <- apply(regimes, 2, function(x) hansen(gsz, tree, factor(x), sqrt.alpha=1, sigma=1))

aic3.dat <- as.data.frame(t(sapply(c(BM=res.bm, res.ou),function(x)unlist(summary(x)[c('aic','aic.c','sic', 'dof')]))))

tree3 <- tree
regimes3 <- regimes
gsz3 <- gsz
res.bm3 <- res.bm
res.ou3 <- res.ou


####
#### Analysis Tree 4
####
tree <- ape2ouch(tree4)
dat <- dat4
rownames(dat) <- dat$nodes
regimes <- dat[c("mpdMpleth", "xMpleth", "MxMpleth", "mpd", "metamorphosis")]
gsz <- log(dat['genomesize'])
res.bm <- brown(gsz,  tree)
res.ou <- apply(regimes, 2, function(x) hansen(gsz, tree, factor(x), sqrt.alpha=1, sigma=1))

aic4.dat <- as.data.frame(t(sapply(c(BM=res.bm, res.ou),function(x)unlist(summary(x)[c('aic','aic.c','sic', 'dof')]))))

tree4 <- tree
regimes4 <- regimes
gsz4 <- gsz
res.bm4 <- res.bm
res.ou4 <- res.ou

print(aic1.dat, digits=4)
print(aic2.dat, digits=4)
print(aic3.dat, digits=4)
print(aic4.dat, digits=4)

aic1 <- cbind(aic1.dat, tree="av_ultra_fulldataset.nex")
aic2 <- cbind(aic2.dat, tree="av_ultra_nomiss.nex")
aic3 <- cbind(aic3.dat, tree="ultra_av_fulldataset.nex")
aic4 <- cbind(aic4.dat, tree="ultra_av_nomiss.nex")
aics <- rbind(aic1, aic2, aic3, aic4)
write.csv(aics, "modelfits.csv")

######
#
#   Plot genome size
#
######

quartz()
plot(gsz1[[1]], regimes1$mpdMpleth)
quartz()
quartz()
hist(gsz1[[1]])
quartz()
qqnorm(gsz1[[1]])

######
#
#   Get Parameters, IC output
#
######

return.param <- function( results.ou, results.bm ) {
  param.names <- c('sigma.sq.matrix', 'alpha.matrix', 'theta.genomesize.D', 'theta.genomesize.M', 'theta.genomesize.P', 'theta.genomesize.pleth', 'theta.genomesize.Mpleth', 'theta.genomesize.Ppleth', 'theta.genomesize.x', 'theta.genomesize.global')
  ou <- as.data.frame(t(sapply(results.ou, function(x) unlist(coef(x))[param.names])))
  bm <- as.data.frame(t(unlist(coef(results.bm))[param.names]))
  rownames(bm) <- "bm"
  names(ou) <- names(bm) <- param.names
  out <- rbind( bm, ou)
  out$model <- rownames(out)
  rownames(out) <- NULL
  return( out[,c(8, 1:7)] )
}

param1 <- cbind(return.param(res.ou1, res.bm1), treename="av_ultra_fulldataset.nex", tree="tree1")
#param2 <- cbind(return.param(res.ou2, res.bm2), treename="av_ultra_nomiss.nex", tree="tree2")
#param3 <- cbind(return.param(res.ou3, res.bm3), treename="ultra_av_fulldataset.nex", tree="tree3")
#param4 <- cbind(return.param(res.ou4, res.bm4), treename="ultra_av_nomiss.nex", tree="tree4")
#parameters <- rbind(param1, param2, param3, param4)
parameters <- param1
parameters$model <- c("bm", names(res.ou1))

write.csv(parameters, "parameter.estimates.csv")

best.params <- parameters[parameters$model=="mpdMpleth" & parameters$tree!="tree3" | parameters$model=="xMpleth" & parameters$tree=="tree3",]
best.params$parameter <- "estimate"
names(best.params) <- sub(".matrix", "", names(best.params))
names(best.params) <- sub(".sq", ".squared", names(best.params))
names(best.params) <- sub("theta", "optima", names(best.params))
best.params <- best.params[,c(10, 11, 1:9)]   # get columns in right order


######
#
#   Parameter Estimate Bootstraps
#
######

# Tree 1 aic: mpdMpleth, sic: BM
# Tree 2 aic: mpdMpleth, sic: BM
# Tree 3 aic: xMpleth, sic: BM
# Tree 4 aic: mpdMpleth, sic: BM

boot1 <- bootstrap(res.ou1[["mpdMpleth"]], nboot=2000)
boot2 <- bootstrap(res.ou2[["mpdMpleth"]], nboot=2000)
boot3 <- bootstrap(res.ou3[["xMpleth"]], nboot=2000)
boot4 <- bootstrap(res.ou4[["mpdMpleth"]], nboot=2000)

boot1.bm <- bootstrap(res.bm1, nboot=2000)
boot2.bm <- bootstrap(res.bm2, nboot=2000)
boot3.bm <- bootstrap(res.bm3, nboot=2000)
boot4.bm <- bootstrap(res.bm4, nboot=2000)

boot1.ci <- sapply(boot1[1:6], function(x) quantile(x, c(0.025, 0.975)))
boot2.ci <- sapply(boot2[1:6], function(x) quantile(x, c(0.025, 0.975)))
boot3.ci <- sapply(boot3[1:4], function(x) quantile(x, c(0.025, 0.975)))
boot4.ci <- sapply(boot4[1:6], function(x) quantile(x, c(0.025, 0.975)))

boot1.bm.ci <- quantile(boot1.bm$sigma.squared, c(0.025, 0.975))
boot2.bm.ci <- quantile(boot2.bm$sigma.squared, c(0.025, 0.975))
boot3.bm.ci <- quantile(boot3.bm$sigma.squared, c(0.025, 0.975))
boot4.bm.ci <- quantile(boot4.bm$sigma.squared, c(0.025, 0.975))

make.ci.data <- function (bootx, treex) {
	return( cbind(tree=treex, parameter=rownames(bootx), as.data.frame(bootx)))
}
boot1.ci <- make.ci.data( boot1.ci, "tree1")
boot2.ci <- make.ci.data( boot2.ci, "tree2")
boot3.ci <- make.ci.data( boot3.ci, "tree3")
boot4.ci <- make.ci.data( boot4.ci, "tree4")

parameters.boot.ci <- merge( rbind(boot1.ci, boot2.ci, boot4.ci), boot3.ci, all=T )
parameters.boot.ci <- parameters.boot.ci[c(1:4, 7, 8, 5, 6), ]
parameters.boot.ci$model <- rep(c("mpdMpleth", "mpdMpleth", "xMpleth", "mpdMpleth"), each=2)

param.out <- merge(best.params, parameters.boot.ci, all=T)  # merge the estimates with 95% CI's
param.out <- param.out[ c(3,1,2,6,4,5,9,7,8,12,10,11),]    # get the rows in the right order

write.csv(param.out, "parameters.csv")


save( boot1, boot2, boot3, boot4, boot1.bm, boot2.bm, boot3.bm, boot4.bm, file="Rdata/parameterboots.rda")
save( boot1.ci, boot2.ci, boot3.ci, boot4.ci, boot1.bm.ci, boot2.bm.ci, boot3.bm.ci, boot4.bm.ci, file="Rdata/parameterboots.ci.rda")




######
#
#   Model Selection Bootstraps  for tree 1
#
######

tree1.mpdMpleth.sims <- simulate(res.ou1[["mpdMpleth"]], nsim=2000)
nnames <- rownames(regimes1)

bm.sim.fit <- lapply(tree1.mpdMpleth.sims, function(y) brown(y, tree1) )

mpdMpleth.fit <- lapply(tree1.mpdMpleth.sims, function(y) hansen(y, tree1, sqrt.alpha=1, sigma=1, regimes=regimes1["mpdMpleth"]))
save(mpdMpleth.fit, file="Rdata/mpdMpleth.sims.fit.rda")
rm(mpdMpleth.fit, bm.sim.fit)

xpleth.fit <- lapply(tree1.mpdMpleth.sims, function(y) hansen(y, tree1, sqrt.alpha=1, sigma=1, regimes=regimes1["xpleth"]))
save(xpleth.fit, file="Rdata/xpleth.sims.fit.rda")
rm(xpleth.fit)

xMpleth.fit <- lapply(tree1.mpdMpleth.sims, function(y) hansen(y, tree1, sqrt.alpha=1, sigma=1, regimes=regimes1["xMpleth"]))
save(xMpleth.fit, file="Rdata/xMpleth.sims.fit.rda")
rm(xMpleth.fit)

xplethMpleth.fit <- lapply(tree1.mpdMpleth.sims, function(y) hansen(y, tree1, sqrt.alpha=1, sigma=1, regimes=regimes1["xplethMpleth"]))
save(xplethMpleth.fit, file="Rdata/xplethMpleth.sims.fit.rda")
rm(xplethMpleth.fit)

xplethMplethPpleth.fit <- lapply(tree1.mpdMpleth.sims, function(y) hansen(y, tree1, sqrt.alpha=1, sigma=1, regimes=regimes1["xplethMplethPpleth"]))
save(xplethMplethPpleth.fit, file="Rdata/xplethMplethPpleth.sims.fit.rda")
rm(xplethMplethPpleth.fit)

MxMpleth.fit <- lapply(tree1.mpdMpleth.sims, function(y) hansen(y, tree1, sqrt.alpha=1, sigma=1, regimes=regimes1["MxMpleth"]))
save(MxMpleth.fit, file="Rdata/MxMpleth.sims.fit.rda")
rm(MxMpleth.fit)

mpd.fit <- lapply(tree1.mpdMpleth.sims, function(y) hansen(y, tree1, sqrt.alpha=1, sigma=1, regimes=regimes1["mpd"]))
save(mpd.fit, file="Rdata/mpd.sims.fit.rda")
rm(mpd.fit)

metamorphosis.fit <- lapply(tree1.mpdMpleth.sims, function(y) hansen(y, tree1, sqrt.alpha=1, sigma=1, regimes=regimes1["metamorphosis"]))
save(metamorphosis.fit, file="Rdata/metamorphosis.sims.fit.rda")
rm(metamorphsis.fit)

regimes <- dat[c("metamorphosis", "mpd", "mpdMpleth", "xpleth", "xplethMpleth", "xplethMplethPpleth", "xMpleth", "OU1")]


regimes <- dat[c("mpdMpleth", "xMpleth", "MxMpleth", "mpd", "metamorphosis")]

aic <- vector(mode="list")
aicc <- vector(mode="list")
sic <- vector(mode="list")

load("Rdata/bm.sim.fit.rda")
fits <- bm.sim.fit
aic[[1]] <- sapply(fits, function(x) unlist(summary(x)['aic']))
aicc[[1]] <- sapply(fits, function(x) unlist(summary(x)['aic.c']))
sic[[1]] <- sapply(fits, function(x) unlist(summary(x)['sic']))
rm(bm.sim.fit)

load("Rdata/mpdMpleth.sims.fit.rda")
fits <- mpdMpleth.fit
aic[[2]] <- sapply(fits, function(x) unlist(summary(x)['aic']))
aicc[[2]] <- sapply(fits, function(x) unlist(summary(x)['aic.c']))
sic[[2]] <- sapply(fits, function(x) unlist(summary(x)['sic']))
rm(mpdMpleth.fit)

load("Rdata/xMpleth.sims.fit.rda")
fits <- xMpleth.fit
aic[[3]] <- sapply(fits, function(x) unlist(summary(x)['aic']))
aicc[[3]] <- sapply(fits, function(x) unlist(summary(x)['aic.c']))
sic[[3]] <- sapply(fits, function(x) unlist(summary(x)['sic']))
rm(xMpleth.fit)

load("Rdata/MxMpleth.sims.fit.rda")
fits <- MxMpleth.fit
aic[[4]] <- sapply(fits, function(x) unlist(summary(x)['aic']))
aicc[[4]] <- sapply(fits, function(x) unlist(summary(x)['aic.c']))
sic[[4]] <- sapply(fits, function(x) unlist(summary(x)['sic']))
rm(MxMpleth.fit)

load("Rdata/mpd.sims.fit.rda")
fits <- mpd.fit
aic[[5]] <- sapply(fits, function(x) unlist(summary(x)['aic']))
aicc[[5]] <- sapply(fits, function(x) unlist(summary(x)['aic.c']))
sic[[5]] <- sapply(fits, function(x) unlist(summary(x)['sic']))
rm(mpd.fit)

load("Rdata/metamorphosis.sims.fit.rda")
fits <- metamorphosis.fit
aic[[6]] <- sapply(fits, function(x) unlist(summary(x)['aic']))
aicc[[6]] <- sapply(fits, function(x) unlist(summary(x)['aic.c']))
sic[[6]] <- sapply(fits, function(x) unlist(summary(x)['sic']))
rm(metamorphosis.fit)

names(aic) <- names(aicc) <- names(sic) <- c('bm', 'mpdMpleth', 'xMpleth', 'MxMpleth', 'mpd', 'metamorphosis')
aic <- as.data.frame(aic)
aicc <- as.data.frame(aicc)
sic <- as.data.frame(sic)

min.aic <- t(apply( aic, 1, function(x) x <= 2+min(x)))
min.aicc <- t(apply( aicc, 1, function(x) x <= 2+min(x)))
min.sic <- t(apply( sic, 1, function(x) x <= 2+min(x)))

modsel.aic.win2 <- colSums(min.aic)/2000
modsel.aicc.win2 <- colSums(min.aicc)/2000
modsel.sic.win2 <- colSums(min.sic)/2000

## number of times that each model is the best model

min.aic <- t(apply( aic, 1, function(x) x == min(x)))
min.aicc <- t(apply( aicc, 1, function(x) x == min(x)))
min.sic <- t(apply( sic, 1, function(x) x == min(x)))

modsel.aic <- colSums(min.aic)/2000
modsel.aicc <- colSums(min.aicc)/2000
modsel.sic <- colSums(min.sic)/2000

#> modsel.aic
#           bm     mpdMpleth       xMpleth      MxMpleth           mpd metamorphosis
#       0.0740        0.6995        0.0845        0.0845        0.1040        0.0380
#> modsel.aicc
#           bm     mpdMpleth       xMpleth      MxMpleth           mpd metamorphosis
#       0.0870        0.6665        0.0885        0.0885        0.1125        0.0455
#> modsel.sic
#           bm     mpdMpleth       xMpleth      MxMpleth           mpd metamorphosis
#       0.4960        0.2725        0.0655        0.0655        0.0960        0.0700

## number of times the best model is more than 2 better than the next best

min.aic <- t(apply( aic, 1, function(x) { oo <- order(x); x + 2 <= x[oo[2]] }))
min.aicc <- t(apply( aicc, 1, function(x) { oo <- order(x); x + 2 <= x[oo[2]] }))
min.sic <- t(apply( sic, 1, function(x) { oo <- order(x); x + 2 <= x[oo[2]] }))

modsel.aic2 <- colSums(min.aic)/2000
modsel.aicc2 <- colSums(min.aicc)/2000
modsel.sic2 <- colSums(min.sic)/2000

#> modsel.aic2
#           bm     mpdMpleth       xMpleth      MxMpleth           mpd metamorphosis
#       0.0365        0.4825        0.0000        0.0000        0.0000        0.0000
#> modsel.aicc2
#           bm     mpdMpleth       xMpleth      MxMpleth           mpd metamorphosis
#       0.0435        0.4505        0.0000        0.0000        0.0150        0.0010
#> modsel.sic2
#           bm     mpdMpleth       xMpleth      MxMpleth           mpd metamorphosis
#       0.3900        0.1725        0.0000        0.0000        0.0340        0.0220

save( aic, aicc, sic, min.aic, min.aicc, min.sic, modsel.aic, modsel.aicc, modsel.sic, file="Rdata/modelfits.summaries.rda")


## To dos
# 1. test additional trees
# 2. get model parameters
# 3. Run boostraps on model fits and parameters
# 4. Explore metamorphosis as ancestral for plethodontids

#lapply(res.ou, summary)
write.csv(aic.dat, file="Rdata/pleth_aic.csv")
write.csv(pleth, file="Rdata/pleth_data.csv")
save(res.bm, res.ou, file="Rdata/pleth_bmresults.rda")
save(tree.ape, file="Rdata/pleth_apetree.rda")
save(bm.sim.fit, file="Rdata/bm.sim.fit.rda")
save(tree1.mpdMpleth.sims, file="Rdata/tree1.mpdMpleth.sims.rda")
