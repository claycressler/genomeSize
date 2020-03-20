## From Jeremy Beaulieu, Mar. 19, 2020
#Necessary packages:
library(OUwie)
library(geiger)

## Load a simulated dataset from OUwie:
data(tworegime)

## Run OUwie:
m1 <- OUwie(tree, trait, model="OU1")

## Now geiger
data.sort <- data.frame(vals1=trait[,3], val2=trait[,3], row.names=trait[,1])
data.sort <- data.sort[tree$tip.label,]
dat <- data.sort
m2 <- fitContinuous(tree, dat["vals1"], model="OU")

## Should be near identical (some differences due to optimization routine used)

## Compare with ouch
library(ouch)
outree <- ape2ouch(tree, scale=FALSE)
