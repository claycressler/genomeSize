require(ape)
require(ouch)

dat1 <- read.csv("tree1.dat.csv")
tree1 <- read.nexus("av_ultra_fulldataset.nex")

####
#### Analysis Tree 1
####
tree <- ape2ouch(tree1)
dat <- dat1
rownames(dat) <- dat$nodes
regimes <- dat[c("metamorphosis", "mpd", "mpdMpleth", "MxMpleth", "xMpleth", "xpleth", "xplethMpleth", "xplethMplethPpleth", "OU1")]
gsz <- log(dat['genomesize'])
res.bm <- brown(gsz,  tree)
res.ou <- apply(regimes, 2, function(x) hansen(gsz, tree, factor(x), sqrt.alpha=1, sigma=1))

aic1.dat <- as.data.frame(t(sapply(c(BM=res.bm, res.ou),function(x)unlist(summary(x)[c('aic','aic.c','sic', 'dof')]))))

tree1 <- tree
regimes1 <- regimes
gsz1 <- gsz
res.bm1 <- res.bm
res.ou1 <- res.ou

print(aic1.dat, digits=4)
aic1 <- cbind(aic1.dat, tree="av_ultra_fulldataset.nex")
write.csv(aic1, "modelfits_tree1.csv")



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

parameters <- cbind(return.param(res.ou1, res.bm1), treename="av_ultra_fulldataset.nex", tree="tree1")
parameters$model <- c("bm", names(res.ou1))
write.csv(parameters, "parameter_estimates_tree1.csv")

best.params <- parameters[parameters$model=="mpdMpleth" & parameters$tree!="tree3",]
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

# Tree 1 aic & aicc: mpdMpleth, sic: BM

boot1 <- bootstrap(res.ou1[["mpdMpleth"]], nboot=2000)
boot1.bm <- bootstrap(res.bm1, nboot=2000)

boot1.ci <- sapply(boot1[1:6], function(x) quantile(x, c(0.025, 0.975)))
boot1.bm.ci <- quantile(boot1.bm$sigma.squared, c(0.025, 0.975))

make.ci.data <- function (bootx, treex) {
	return( cbind(tree=treex, parameter=rownames(bootx), as.data.frame(bootx)))
}
boot1.ci <- make.ci.data( boot1.ci, "tree1")

parameters.boot.ci <- boot1.ci[ ,c(1:4, 7, 8, 5, 6) ]
parameters.boot.ci$model <- c("mpdMpleth", "mpdMpleth")

param.out <- merge(best.params, parameters.boot.ci, all=T)  # merge the estimates with 95% CI's
#param.out <- param.out[ c(3,1,2,6,4,5,9,7,8,12,10,11),]    # get the rows in the right order

write.csv(param.out, "parameters_tree1.csv")


save( boot1, boot1.bm, file="Rdata/parameterboots_tree1.rda")
save( boot1.ci, boot1.bm.ci, file="Rdata/parameterboots_ci_tree1.rda")


######
#
#   Model Selection Bootstraps  for tree 1
#
######

tree1.mpdMpleth.sims <- simulate(res.ou1[["mpdMpleth"]], nsim=2000)
save(tree1.mpdMpleth.sims, file="Rdata/tree1.mpdMpleth.sims")

bm.fit <- lapply(tree1.mpdMpleth.sims, function(y) brown(y, tree1) )
save(bm.fit, file="Rdata/bm.sims.fit.rda")
rm(bm.fit)

metamorphosis.fit <- lapply(tree1.mpdMpleth.sims, function(y) hansen(y, tree1, sqrt.alpha=1, sigma=1, regimes=regimes1["metamorphosis"]))
save(metamorphosis.fit, file="Rdata/metamorphosis.sims.fit.rda")
rm(metamorphosis.fit)

mpd.fit <- lapply(tree1.mpdMpleth.sims, function(y) hansen(y, tree1, sqrt.alpha=1, sigma=1, regimes=regimes1["mpd"]))
save(mpd.fit, file="Rdata/mpd.sims.fit.rda")
rm(mpd.fit)

mpdMpleth.fit <- lapply(tree1.mpdMpleth.sims, function(y) hansen(y, tree1, sqrt.alpha=1, sigma=1, regimes=regimes1["mpdMpleth"]))
save(mpdMpleth.fit, file="Rdata/mpdMpleth.sims.fit.rda")
rm(mpdMpleth.fit, bm.sim.fit)

MxMpleth.fit <- lapply(tree1.mpdMpleth.sims, function(y) hansen(y, tree1, sqrt.alpha=1, sigma=1, regimes=regimes1["MxMpleth"]))
save(MxMpleth.fit, file="Rdata/MxMpleth.sims.fit.rda")
rm(MxMpleth.fit)

xMpleth.fit <- lapply(tree1.mpdMpleth.sims, function(y) hansen(y, tree1, sqrt.alpha=1, sigma=1, regimes=regimes1["xMpleth"]))
save(xMpleth.fit, file="Rdata/xMpleth.sims.fit.rda")
rm(xMpleth.fit)

xpleth.fit <- lapply(tree1.mpdMpleth.sims, function(y) hansen(y, tree1, sqrt.alpha=1, sigma=1, regimes=regimes1["xpleth"]))
save(xpleth.fit, file="Rdata/xpleth.sims.fit.rda")
rm(xpleth.fit)

xplethMpleth.fit <- lapply(tree1.mpdMpleth.sims, function(y) hansen(y, tree1, sqrt.alpha=1, sigma=1, regimes=regimes1["xplethMpleth"]))
save(xplethMpleth.fit, file="Rdata/xplethMpleth.sims.fit.rda")
rm(xplethMpleth.fit)

xplethMplethPpleth.fit <- lapply(tree1.mpdMpleth.sims, function(y) hansen(y, tree1, sqrt.alpha=1, sigma=1, regimes=regimes1["xplethMplethPpleth"]))
save(xplethMplethPpleth.fit, file="Rdata/xplethMplethPpleth.sims.fit.rda")
rm(xplethMplethPpleth.fit)

ou1.fit <- lapply(tree1.mpdMpleth.sims, function(y) hansen(y, tree1, sqrt.alpha=1, sigma=1, regimes=regimes1["OU1"]))
save(ou1.fit, file="Rdata/ou1.sims.fit.rda")
rm(ou1.fit)


######
#
#   Assemble AIC, AICC and SIC 
#
######


aic <- vector(mode="list", length=10)
aicc <- vector(mode="list", length=10)
sic <- vector(mode="list", length=10)
names(aic) <- names(aicc) <- names(sic) <- c('bm', 'metamorphosis', 'mpd', 'mpdMpleth', 'MxMpleth', 'xMpleth', 'xpleth', 'xplethMpleth', 'xplethMplethPpleth', "ou1" )

load("Rdata/bm.sims.fit.rda")
fits <- bm.fit
aic[['bm']] <- sapply(fits, function(x) unlist(summary(x)['aic']))
aicc[['bm']] <- sapply(fits, function(x) unlist(summary(x)['aic.c']))
sic[['bm']] <- sapply(fits, function(x) unlist(summary(x)['sic']))
rm(bm.fit)

load("Rdata/metamorphosis.sims.fit.rda")
fits <- metamorphosis.fit
aic[['metamorphosis']] <- sapply(fits, function(x) unlist(summary(x)['aic']))
aicc[['metamorphosis']] <- sapply(fits, function(x) unlist(summary(x)['aic.c']))
sic[['metamorphosis']] <- sapply(fits, function(x) unlist(summary(x)['sic']))
rm(metamorphosis.fit)

load("Rdata/mpd.sims.fit.rda")
fits <- mpd.fit
aic[['mpd']] <- sapply(fits, function(x) unlist(summary(x)['aic']))
aicc[['mpd']] <- sapply(fits, function(x) unlist(summary(x)['aic.c']))
sic[['mpd']] <- sapply(fits, function(x) unlist(summary(x)['sic']))
rm(mpd.fit)

load("Rdata/mpdMpleth.sims.fit.rda")
fits <- mpdMpleth.fit
aic[['mpdMpleth']] <- sapply(fits, function(x) unlist(summary(x)['aic']))
aicc[['mpdMpleth']] <- sapply(fits, function(x) unlist(summary(x)['aic.c']))
sic[['mpdMpleth']] <- sapply(fits, function(x) unlist(summary(x)['sic']))
rm(mpdMpleth.fit)

load("Rdata/MxMpleth.sims.fit.rda")
fits <- MxMpleth.fit
aic[['MxMpleth']] <- sapply(fits, function(x) unlist(summary(x)['aic']))
aicc[['MxMpleth']] <- sapply(fits, function(x) unlist(summary(x)['aic.c']))
sic[['MxMpleth']] <- sapply(fits, function(x) unlist(summary(x)['sic']))
rm(MxMpleth.fit)

load("Rdata/xMpleth.sims.fit.rda")
fits <- xMpleth.fit
aic[['xMpleth']] <- sapply(fits, function(x) unlist(summary(x)['aic']))
aicc[['xMpleth']] <- sapply(fits, function(x) unlist(summary(x)['aic.c']))
sic[['xMpleth']] <- sapply(fits, function(x) unlist(summary(x)['sic']))
rm(xMpleth.fit)

load("Rdata/xpleth.sims.fit.rda")
fits <- xpleth.fit
aic[['xpleth']] <- sapply(fits, function(x) unlist(summary(x)['aic']))
aicc[['xpleth']] <- sapply(fits, function(x) unlist(summary(x)['aic.c']))
sic[['xpleth']] <- sapply(fits, function(x) unlist(summary(x)['sic']))
rm(xpleth.fit)

load("Rdata/xplethMpleth.sims.fit.rda")
fits <- xplethMpleth.fit
aic[['xplethMpleth']] <- sapply(fits, function(x) unlist(summary(x)['aic']))
aicc[['xplethMpleth']] <- sapply(fits, function(x) unlist(summary(x)['aic.c']))
sic[['xplethMpleth']] <- sapply(fits, function(x) unlist(summary(x)['sic']))
rm(xplethMpleth.fit)

load("Rdata/xplethMplethPpleth.sims.fit.rda")
fits <- xplethMplethPpleth.fit
aic[['xplethMplethPpleth']] <- sapply(fits, function(x) unlist(summary(x)['aic']))
aicc[['xplethMplethPpleth']] <- sapply(fits, function(x) unlist(summary(x)['aic.c']))
sic[['xplethMplethPpleth']] <- sapply(fits, function(x) unlist(summary(x)['sic']))
rm(xplethMplethPpleth.fit)

load("Rdata/ou1.sims.fit.rda")
fits <- ou1.fit
aic[['ou1']] <- sapply(fits, function(x) unlist(summary(x)['aic']))
aicc[['ou1']] <- sapply(fits, function(x) unlist(summary(x)['aic.c']))
sic[['ou1']] <- sapply(fits, function(x) unlist(summary(x)['sic']))
rm(ou1.fit)


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
#bm		metamorphosis	mpd		mpdMpleth	MxMpleth	xMpleth	xpleth	xplethMpleth	xplethMplethPpleth	ou1 
#0.0585	0.0335			0.1040	0.6090		0.0380		0.0625	0.0060	0.0625			0.0260				0.0000 
#> modsel.aicc
#bm		metamorphosis	mpd		mpdMpleth	MxMpleth	xMpleth	xpleth	xplethMpleth	xplethMplethPpleth	ou1 
#0.0735	0.0405			0.1130	0.5680		0.0410		0.0705	0.0065	0.0640			0.0230				0.0000 
#> modsel.sic
#bm      metamorphosis	mpd		mpdMpleth	MxMpleth	xMpleth	xpleth	xplethMpleth	xplethMplethPpleth	ou1 
#0.4005	0.0730			0.0915	0.2205		0.0095		0.1440	0.0070	0.0490			0.0050				0.0000 

## number of times the best model is more than 2 better than the next best

min.aic <- t(apply( aic, 1, function(x) { oo <- order(x); x + 2 <= x[oo[2]] }))
min.aicc <- t(apply( aicc, 1, function(x) { oo <- order(x); x + 2 <= x[oo[2]] }))
min.sic <- t(apply( sic, 1, function(x) { oo <- order(x); x + 2 <= x[oo[2]] }))

modsel.aic2 <- colSums(min.aic)/2000
modsel.aicc2 <- colSums(min.aicc)/2000
modsel.sic2 <- colSums(min.sic)/2000

#> modsel.aic2
#bm      metamorphosis	mpd		mpdMpleth	MxMpleth	xMpleth	xpleth	xplethMpleth	xplethMplethPpleth	ou1 
#0.0245	0.0000			0.0000	0.3690		0.0000		0.0000	0.0000	0.0000			0.0080				0.0000 
#> modsel.aicc2
#bm      metamorphosis	mpd		mpdMpleth	MxMpleth	xMpleth	xpleth	xplethMpleth	xplethMplethPpleth	ou1 
#0.0295	0.0010			0.0180	0.3415		0.0010		0.0015	0.0010	0.0080			0.0070				0.0000 
#> modsel.sic2
#bm      metamorphosis	mpd		mpdMpleth	MxMpleth	xMpleth	xpleth	xplethMpleth	xplethMplethPpleth	ou1 
#0.2745	0.0180			0.0385	0.1280		0.0010		0.0525	0.0015	0.0150			0.0015				0.0000 

save( aic, aicc, sic, min.aic, min.aicc, min.sic, modsel.aic, modsel.aicc, modsel.sic, file="Rdata/modelfits.summaries.rda")
