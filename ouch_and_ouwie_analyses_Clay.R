require(ape)
require(ouch)
require(OUwie)

##############################################################
## OUwie setup
##############################################################
tree <- read.nexus("av_ultra_fulldataset.nex")
tree_metamorphosis <- tree_mpd <- tree_mpdMpleth <- tree
## Load the hypotheses
nodedata <- read.csv('ouwie_nodelabels.csv')
nodedata <- nodedata[,c('node','metamorphosis','mpd','mpdMpleth')]
## metamorphosis has two regimes: M, x
## mpd has three regimes: D, M, P
## mpdMpleth has four regimes: D, M, Mpleth, P
## Node Colors for plots
nodecol <- nodedata
nodecol[-1] <- sapply(nodecol[-1], as.character)
nodecol[nodecol=="x"] <- "brown"
nodecol[nodecol=="M"] <- "red"
nodecol[nodecol=="D"] <- "black"
nodecol[nodecol=="P"] <- "dodgerblue1"
nodecol[nodecol=="Mpleth"] <- "purple"
## code Regime names as numbers for the OUwie algorithm
nn <- nodedata
nn[-1] <- sapply(nn[-1], as.character)
nn[nn=="x"] <- 1
nn[nn=="M"] <- 2
nn[nn=="D"] <- 3
nn[nn=="P"] <- 4
nn[nn=="Mpleth"] <- 5
## Add node labels to trees
tree_metamorphosis$node.label <- nn$metamorphosis
tree_mpd$node.label <- nn$mpd
tree_mpdMpleth$node.label <- nn$mpdMpleth
## Tip data
dat <- read.csv("tree1.dat.csv")
dat$genomesize <- log(dat$genomesize)
tdat <- dat[!is.na(dat$genomesize),]  # tip data
ndat <- dat[is.na(dat$genomesize),]   # node data
## reorder tdat to match order of tip labels in phylo tree
oo <- sapply(tree$tip.label, function(x) grep(x, tdat$labels))
tipdat <- tdat[oo,]
tipdat <- tipdat[c("labels", "metamorphosis", "mpd", "mpdMpleth", "genomesize")]
## Tip Colors for plots
tipcol <- tipdat
tipcol[,1:(ncol(tipdat)-1)] <- sapply(tipcol[,1:(ncol(tipdat)-1)], as.character)
tipcol[tipcol=="x"] <- "brown"
tipcol[tipcol=="M"] <- "red"
tipcol[tipcol=="D"] <- "black"
tipcol[tipcol=="P"] <- "dodgerblue1"
tipcol[tipcol=="Mpleth"] <- "purple"
## code Regime names as numbers for the OUwie algorithm
tt <- tipdat
tt[,1:(ncol(tipdat)-1)] <- sapply(tt[,1:(ncol(tipdat)-1)], as.character)
tt[tt=="x"] <- 1
tt[tt=="M"] <- 2
tt[tt=="D"] <- 3
tt[tt=="P"] <- 4
tt[tt=="Mpleth"] <- 5
### input dataframes - one for each hypothesis
data.metamorphosis <- tt[c("labels", "metamorphosis", "genomesize")]
data.mpd <- tt[c("labels", "mpd", "genomesize")]
data.mpdMpleth <- tt[c("labels", "mpdMpleth", "genomesize")]

## Single-rate BM
mpdMpleth.BM1 <- OUwie(tree_mpdMpleth,data.mpdMpleth,model=c("BM1"),root.station=TRUE)
## Multiple-rate BM
mpdMpleth.BMS <- OUwie(tree_mpdMpleth,data.mpdMpleth,model=c("BMS"),root.station=TRUE)
mpd.BMS <- OUwie(tree_mpd,data.mpd,model=c("BMS"),root.station=TRUE)
meta.BMS <- OUwie(tree_metamorphosis,data.metamorphosis,model=c("BMS"),root.station=TRUE)

##############################################################
## OUCH setup
##############################################################
dat1 <- read.csv("tree1.dat.csv")
tree1 <- read.nexus("av_ultra_fulldataset.nex")
tree <- ape2ouch(tree1)
dat <- dat1
rownames(dat) <- dat$nodes
regimes <- dat[c("metamorphosis", "mpd", "mpdMpleth","OU1","fivereg")]
gsz <- log(dat['genomesize'])

## Fit hansen model to the data
bm.ouch <- brown(gsz, tree)
mpdMpleth.ouch <- hansen(gsz, tree, regimes['mpdMpleth'], sqrt.alpha=1, sigma=1)
mpd.ouch <- hansen(gsz, tree, regimes['mpd'], sqrt.alpha=1, sigma=1)
meta.ouch <- hansen(gsz, tree, regimes['metamorphosis'], sqrt.alpha=1, sigma=1)
ou1.ouch <- hansen(gsz, tree, regimes['OU1'], sqrt.alpha=1, sigma=1)


sims = simulate(mpdMpleth.ouch, nsim=200)
mpdMpleth.boots <- sapply(sims, function(x) summary(hansen(x, tree, regimes['mpdMpleth'], sqrt.alpha=0.1, sigma=0.1))$aic.c)
mpd.boots <- sapply(sims, function(x) summary(hansen(x, tree, regimes['mpd'], sqrt.alpha=0.1, sigma=0.1))$aic.c)
meta.boots <- sapply(sims, function(x) summary(hansen(x, tree, regimes['metamorphosis'], sqrt.alpha=0.1, sigma=0.1))$aic.c)
ou1.boots <- sapply(sims, function(x) summary(hansen(x, tree, regimes['OU1'], sqrt.alpha=0.1, sigma=0.1))$aic.c)
bm.boots <- sapply(sims, function(x) summary(brown(x, tree))$aic.c)

ouwieord <- as.numeric(sapply(data.mpdMpleth$labels, function(x) which(tree@nodelabels==x)))

ouwie.BM1.boots <- vector(length=length(sims))
ouwie.BMS.boots <- vector(length=length(sims))
for (i in 1:length(sims)) {
    this.dat <- data.mpdMpleth
    this.dat$genomesize <- sims[[i]][ouwieord,1]
    ouwie.BMS.boots[i] <- OUwie(tree_mpdMpleth,this.dat,model=c("BMS"),root.station=TRUE)$AICc
    ouwie.BM1.boots[i] <- OUwie(tree_mpdMpleth,this.dat,model=c("BM1"),root.station=TRUE)$AICc
}

bms.mpd.boots <- vector(length=length(sims))
bms.meta.boots <- vector(length=length(sims))
for (i in 1:length(sims)) {
    this.dat <- data.mpd
    this.dat$genomesize <- sims[[i]][ouwieord,1]
    bms.mpd.boots[i] <- OUwie(tree_mpd,this.dat,model=c("BMS"),root.station=TRUE)$AICc
    this.dat <- data.metamorphosis
    this.dat$genomesize <- sims[[i]][ouwieord,1]
    bms.meta.boots[i] <- OUwie(tree_metamorphosis,this.dat,model=c("BMS"),root.station=TRUE)$AICc
}

boots <- data.frame(mpdMpleth=mpdMpleth.boots, mpd=mpd.boots,
                    meta=meta.boots, ou1=ou1.boots,
                    bm=bm.boots, bms.mpdMpleth=ouwie.BMS.boots,
                    bms.mpd=bms.mpd.boots, bms.meta=bms.meta.boots)
save(boots, file='Model_selection_bootstraps.rda')
table(apply(boots, 1, which.min))/200

## Boettiger et al. 2012 phylogenetic Monte Carlo approach
mpdMpleth <- hansen(gsz, tree, regimes['mpdMpleth'], sqrt.alpha=0.1, sigma=0.1)
mpd <- hansen(gsz, tree, regimes['mpd'], sqrt.alpha=0.1, sigma=0.1)
meta <- hansen(gsz, tree, regimes['metamorphosis'], sqrt.alpha=0.1, sigma=0.1)
ou1 <- hansen(gsz, tree, regimes['OU1'], sqrt.alpha=0.1, sigma=0.1)
bm <- brown(gsz, tree)

## test statistic
lr <- -2*(summary(bm)$loglik-summary(mpdMpleth)$loglik)

## to estimate the distribution of the test statistic -2*(log L0-log L1), generate 1000 datasets under the simpler model
simdata <- simulate(bm, nsim=1000)
bm.fits <- lapply(simdata, function(x) brown(x, tree))
mpdMpleth.fits <- lapply(simdata, function(x) hansen(x, tree, regimes['mpdMpleth'], sqrt.alpha=0.1, sigma=0.1))
lr.dist <- sapply(1:length(bm.fits), function(x) -2*(summary(bm.fits[[x]])$loglik-summary(mpdMpleth.fits[[x]])$loglik))
## choose a value of the test statistic so that 95% of the simulated values fall below it - reject model 0 (BM) if the observed value of the test statistic is greater than this value
lr > sort(lr.dist)[950] ## true
## put another way, the proportion of the simulated values larger than the test statistic provides an approximation to the p-value for the test, the probability that a difference at least as large would be seen under model 0 (BM)
sum(lr < lr.dist)/length(lr.dist) ## "p-value" of 0.042

## repeat the same process, simulating under the more complex model to get a monte carlo estimate of power
simdata.2 <- simulate(mpdMpleth, nsim=1000)
bm.fits.2 <- lapply(simdata.2, function(x) brown(x, tree))
mpdMpleth.fits.2 <- lapply(simdata.2, function(x) hansen(x, tree, regimes['mpdMpleth'], sqrt.alpha=0.1, sigma=0.1))
lr.dist.2 <- sapply(1:length(bm.fits.2), function(x) -2*(summary(bm.fits.2[[x]])$loglik-summary(mpdMpleth.fits.2[[x]])$loglik))
## the amount of this distribution to the right of the value chosen above approximates the probability of rejecting model 0 when the data are produced by model 1 (e.g., power)
sum(lr.dist.2 > sort(lr.dist)[950])/length(lr.dist.2) ## 0.86

## plot this to show it another way
par(mar=c(5,5,1,1), oma=rep(0,4))
d1 <- density(lr.dist)
d2 <- density(lr.dist.2)
plot.new()
plot.window(xlim=range(c(d1$x,d2$x)), ylim=range(c(d1$y,d2$y)))
axis(1);axis(2);box('plot')
mtext(side=1, line=3, expression(delta))
mtext(side=2, line=3, 'Density')
lines(d1$x, d1$y, col=gray(0.7))
with(d1, polygon(x=c(min(x), x[x<max(x)], rev(x[x<max(x)]), min(x)), y=c(0, 0*x[x<max(x)], rev(y[x<max(x)]), 0), col=gray(0.7, alpha=0.5)))
lines(d2$x, d2$y, col=gray(0.2))
with(d2, polygon(x=c(min(x), x[x<max(x)], rev(x[x<max(x)]), min(x)), y=c(0, 0*x[x<max(x)], rev(y[x<max(x)]), 0), col=gray(0.2, alpha=0.5)))
abline(v=lr, lwd=2)
legend(x='topright', c('BM','mpdMpleth'), fill=c(gray(0.7),gray(0.2)), bty='n')


###############################################################
###  mpdMpleth
###############################################################
## Confirm that the two paintings are identical
plot(tree_mpdMpleth, cex=.6)     # plot tree with slightly reduced font size
nodelabels(pch=21, bg=nodecol$mpdMpleth)   # plot node lables
tiplabels(pch=21, bg=tipcol$mpdMpleth)   # plot tip labels
title("mpdMpleth")
quartz()
plot(tree, regimes=regimes['mpdMpleth'])

## And check that the tip assignation is the same
ouch.tip.reg <- dat1[tree@term,c('labels','mpdMpleth')]
sort.order <- as.numeric(sapply(as.character(ouch.tip.reg[,1]), function(x) which(tipdat[,'labels']==x)))
ouwie.tip.reg <- tipdat[sort.order,c('labels','mpdMpleth')]
cbind(ouch.tip.reg, ouwie.tip.reg[,2])
ouch.tip.reg[,2]==ouwie.tip.reg[,2]




## Update the estimate, using subplex for the fit
mpdMpleth1.up <- update(mpdMpleth1, method='subplex')
## You can see that the update caused some very minor changes in a,
## s^2, and the optima.
unlist(summary(mpdMpleth1)[4:6])-unlist(summary(mpdMpleth1.up)[4:6])

fivereg1 <- hansen(gsz, tree, regimes['fivereg'], sqrt.alpha=1, sigma=1)

## calculate eta, phi, snr
dtheta=diff(summary(mpdMpleth1)$optima$genomesize[c('D','M')])
phi <- 2*summary(mpdMpleth1)$alpha*dtheta/summary(mpdMpleth1)$sigma.squared
eta <- summary(mpdMpleth1)$alpha
snr <- sqrt(eta)*phi
## so we are definitely in a region where you would expect BM to be frequently chosen, and yet AIC favors a more complicated model

## modify very slightly by repainting some of the internal branches
## leading to metamorphic plethodontids
regimes2 <- regimes
regimes2[c(45,46,68),'mpdMpleth'] <- as.factor('D')
mpdMpleth2 <- hansen(gsz, tree, regimes2['mpdMpleth'], sqrt.alpha=1, sigma=1)
## No real change here

## Drop hemidactylium scudatum from the analysis
sub <- c(90:96,98:195)
subtree <- subsample.tree(tree, sub)
subreg <- regimes2
subreg <- subreg[-which(rownames(regimes2)=='97'),]
subgsz <- data.frame(genomesize=gsz[-which(rownames(gsz)=='97'),])
rownames(subgsz) <- rownames(gsz)[-which(rownames(gsz)=='97')]
subgsz <- gsz[-which(rownames(gsz)=='97'),]
mpdMpleth.sub <- hansen(subgsz, subtree, subreg['mpdMpleth'], sqrt.alpha=1, sigma=1)
## this actually makes the estimate of the Mpleth optimum even smaller; it also really drops the estimate of the D optimum

## look at the genome size relative to the amount of time each tip species has been in its final regime state
gsz.against.time <- array(NA, dim=c(tree@nterm,4))
for (i in 1:tree@nterm) {
    line <- tree@lineages[[tree@term[i]]] ## get the lineage for this tip
    line.reg <- regimes[sapply(as.character(line), function(x) which(rownames(regimes)==x)),'mpdMpleth'] ## get the regimes for this lineage
    ## find the first node in the lineage not in the same regime as the tip and gets its time
    time.in.reg <- tree@times[tree@term[i]]-tree@times[line[which(line.reg!=line.reg[1])[1]]]
    if (is.na(time.in.reg)) time.in.reg <- tree@times[tree@term[i]]
    this.gsz <- gsz[rownames(gsz)==as.character(tree@term[i]),]
    this.sp <- tree@nodelabels[tree@term[i]]
    this.reg <- regimes[which(rownames(regimes)==as.character(tree@term[i])),'mpdMpleth']
    gsz.against.time[i,] <- c(this.sp, as.character(this.reg), round(time.in.reg,3), round(this.gsz,3))
}

## Plot time in regime against gsz for each of the
plot(as.numeric(gsz.against.time[,3]), as.numeric(gsz.against.time[,4]), pch=21, bg=as.numeric(as.factor(gsz.against.time[,2])))
legend(x='topleft', levels(as.factor(gsz.against.time[,2])), fill=1:4)
## interesting. there is definitely an increasing trend - the longer in the regime, the larger the genome. That may be part of the problem for the mPleth - genomes should hypothetically be smaller the longer in the regime, not larger

## although, actually, that isn't really the right thing to be looking at, because it is really about the regime that it transitions from.

## What if you use different starting parameter guesses?
mpdMpleth2 <- hansen(gsz, tree, regimes['mpdMpleth'], sqrt.alpha=0.1, sigma=0.1)
## They converge to almost exactly the same place.
unlist(summary(mpdMpleth1)[4:6])-unlist(summary(mpdMpleth2)[4:6])

## Bootstrap the fit
mpdMpleth.bs <- bootstrap(mpdMpleth1, nboot=1000)
mpdMpleth.ci <- apply(mpdMpleth.bs, 2, function(x) quantile(x, c(0.025, 0.975)))
## The MLEs of alpha and the optima Mpleth and P are not within the
## 95% confidence intervals estimated by bootstrapping.
sapply(names(unlist(summary(mpdMpleth1)[4:6])), function(x) findInterval(unlist(summary(mpdMpleth1)[4:6])[x], mpdMpleth.ci[,x])==1)
## This is very strange. The MLE parameter estimates are used to
## simulate phenotypic datasets, which are then fit using hansen(), so
## this implies that the generating parameter set is not its own MLE!
## That doesn't make any sense.

## Use one of the bootstrapped parameter estimates to fit the true
## genomesize data
mpdMpleth.bs.fit <- hansen(data=gsz, tree=tree, regimes=regimes['mpdMpleth'], sqrt.alpha=sqrt(mpdMpleth.bs$alpha[which.min(mpdMpleth.bs$aic.c)]), sigma=sqrt(mpdMpleth.bs$sigma.squared[which.min(mpdMpleth.bs$aic.c)]))
## There aren't really any differences here.
unlist(summary(mpdMpleth1)[4:6])-unlist(summary(mpdMpleth.bs.fit)[4:6])

## Take one of the bootstrapped parameter estimates and simulate a dataset. Fit hansen() to that simulated dataset to get the MLE parameter estimates. Then bootstrap that fitted hansentree object.
simdata <- simulate(mpdMpleth.bs.fit, nsim=1, seed=10001)
simfit <- hansen(data=simdata$rep.1, tree=tree, regimes=regimes['mpdMpleth'], sqrt.alpha=sqrt(as.numeric(summary(mpdMpleth.bs.fit)$alpha)), sigma=sqrt(as.numeric(summary(mpdMpleth.bs.fit)$sigma.squared)))
boot.simfit <- bootstrap(simfit, nboot=1000)
boot.simfit.ci <- apply(boot.simfit, 2, quantiles(x, c(0.025, 0.975)))
## In this case, the data-generating parameter set is contained within
## the 95% confidence interval of the bootstrapped parameter sets.
sapply(names(unlist(summary(simfit)[4:6])), function(x) findInterval(unlist(summary(simfit)[4:6])[x], boot.simfit.ci[,x])==1)

## So what does this mean? The MLE parameter estimates for the true
## data, when used to simulate new data, are not good estimates of the
## simulated data. However, if you use one of the bootstrap MLEs to
## generate a new simulated dataset, fit it, and then bootstrap that,
## the data-generating parameter set is contained in the 95%
## confidence interval. This suggests to me that there is something
## strange about the original data.

## Either that, or there is something strange about this particular
## regime painting. What if we fit one of the other hypotheses, and
## then repeat this bootstrapping analysis?

###############################################################
###  mpd
###############################################################
## This is a much simpler hypothesis

## Confirm that the two paintings are identical
plot(tree_mpd, cex=.6)     # plot tree with slightly reduced font size
nodelabels(pch=21, bg=nodecol$mpd)   # plot node lables
tiplabels(pch=21, bg=tipcol$mpd)   # plot tip labels
title("mpd")
quartz()
plot(tree, regimes=regimes['mpd'])

## And check that the tip assignation is the same
ouch.tip.reg <- dat1[tree@term,c('labels','mpd')]
sort.order <- as.numeric(sapply(as.character(ouch.tip.reg[,1]), function(x) which(tipdat[,'labels']==x)))
ouwie.tip.reg <- tipdat[sort.order,c('labels','mpd')]
cbind(ouch.tip.reg, ouwie.tip.reg[,2])
ouch.tip.reg[,2]==ouwie.tip.reg[,2]

## Fit hansen model to the data
mpd <- hansen(gsz, tree, regimes['mpd'], sqrt.alpha=1, sigma=1)
## Bootstrap the fit
mpd.bs <- bootstrap(mpd, nboot=1000)
mpd.ci <- apply(mpd.bs, 2, function(x) quantile(x, c(0.025, 0.975)))
## For this hypothesis, the MLEs of alpha and the optima D and P are
## not within the 95% confidence intervals estimated by bootstrapping.
sapply(names(unlist(summary(mpd)[4:6])), function(x) findInterval(unlist(summary(mpd)[4:6])[x], mpd.ci[,x])==1)
## This reinforces the possibility that there is something strange
## about the dataset itself. In particular, for both the mpdMpleth and
## mpd hypotheses, the alpha estimate and the selective optima
## estimate for paedomorphs are not within the 95% confidence interval
## estimated using datasets generated by these parameter values.

###############################################################
###  metamorphosis
###############################################################
## Confirm that the two paintings are identical
plot(tree_metamorphosis, cex=.6)     # plot tree with slightly reduced font size
nodelabels(pch=21, bg=nodecol$metamorphosis)   # plot node lables
tiplabels(pch=21, bg=tipcol$metamorphosis)   # plot tip labels
title("metamorphosis")
quartz()
plot(tree, regimes=regimes['metamorphosis'])

## And check that the tip assignation is the same
ouch.tip.reg <- dat1[tree@term,c('labels','metamorphosis')]
sort.order <- as.numeric(sapply(as.character(ouch.tip.reg[,1]), function(x) which(tipdat[,'labels']==x)))
ouwie.tip.reg <- tipdat[sort.order,c('labels','metamorphosis')]
cbind(ouch.tip.reg, ouwie.tip.reg[,2])
ouch.tip.reg[,2]==ouwie.tip.reg[,2]

## Fit hansen model to the data
meta <- hansen(gsz, tree, regimes['metamorphosis'], sqrt.alpha=1, sigma=1)
## Bootstrap the fit
meta.bs <- bootstrap(meta, nboot=1000)
meta.ci <- apply(meta.bs, 2, function(x) quantile(x, c(0.025, 0.975)))
## For this hypothesis, the MLEs of alpha and the optima x are not in the 95% confidence interval
sapply(names(unlist(summary(meta)[4:6])), function(x) findInterval(unlist(summary(meta)[4:6])[x], meta.ci[,x])==1)

## This is beginning to suggest something. In particular, the two
## regimes in the 'mpd' hypotheses that are combined to create the 'x'
## regime in the 'metamorphosis' hypothesis are 'P' and 'D', which
## were exactly the two regimes whose MLE estimate were not in the 95%
## confidence interval. However, splitting out the Mpleth group from
## the M group somehow allows 'D' to be estimated, which is not
## consistent. Regardless, putting all of this together definitely
## seems to suggest that there is something about the data itself that
## makes the fits behave in a particularly strange way. I wonder if
## that is sort of what OUwie was indicating, with the alpha estimate
## always running away towards the boundary - it was trying to tell
## you that the data do not have a nice MLE.

## Here is what I think: the MLEs from fitting the data are so extreme
## that the data produced from simulating at the MLEs on the tree
## produce data that is not extreme. This is suggested by the CIs -
## the most extreme parameters (alpha, the P optima and the Mpleth
## optima) all fall outside the 95% confidence intervals.

## Bootstrapping
