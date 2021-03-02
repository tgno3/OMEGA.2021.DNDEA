#########################################################################################################################
### Project  : Dynamic Network DEA with distributable carryovers
### Script   : DNDEA on GIT.R
### Contents : Non-cooperative model with dynamic carryovers
#########################################################################################################################

#########################################################################################################################
### Setting up environment
#########################################################################################################################

# Load library
pkgs <- c("DJL")
sapply(pkgs, require, character.only = T)

# Load data & parameters
load("RI.14.18.Rdata")


#########################################################################################################################
### Motivational example
#########################################################################################################################

# Toy data
df.toy <- array(c(10, 10, 10, 10, 10,  5, 
                  10, 10, 10, 10, 10, 15), c(2,3,2),
                dimnames = list(c("A", "B"), c("X", "Z", "Y"), c("T1", "T2")))

# DMU B underestimated in Period 1 while DMU A underestimated in Period 2
dm.dynamic.network(df.toy[,1,], NULL, df.toy[,2,], NULL, df.toy[,3,])$eff.s2


#########################################################################################################################
### Analysis
#########################################################################################################################

# DNDEA
rts <- "vrs"
res.dndea.ns <- dm.dynamic.network(df.3d[,id.x.s1,], df.3d[,id.y.s1,], df.3d[,id.z,], 
                                   df.3d[,id.x.s2,], df.3d[,id.y.s2,], rts = rts)

res.dndea.ss <- dm.dynamic.network(df.3d[,id.x.s1,], df.3d[,id.y.s1,], df.3d[,id.z,], 
                                   df.3d[,id.x.s2,], df.3d[,id.y.s2,], alpha = c(0.4, 0.4, 0.2), rts = rts)

res.dndea.ds <- dm.dynamic.network(df.3d[,id.x.s1,], df.3d[,id.y.s1,], df.3d[,id.z,], 
                                   df.3d[,id.x.s2,], df.3d[,id.y.s2,], alpha = "free", max.cp = 2, rts = rts)

data.frame(No.Split = res.dndea.ns$eff.sys, Static.Split = res.dndea.ss$eff.sys, Dynamic.Split = res.dndea.ds$eff.sys)

data.frame(No.split      = res.dndea.ns$eff.s2[,,3],
           Static.split  = res.dndea.ss$eff.s2[,,3],
           Dynamic.split = res.dndea.ds$eff.s2[,,3],
           DS.win        = ifelse(res.dndea.ss$eff.s2[,,3] <= res.dndea.ds$eff.s2[,,3], 1, 0))

data.frame(No.split      = apply(res.dndea.ns$eff.s2[,,1:5], 1 , "mean"),
           Static.split  = apply(res.dndea.ss$eff.s2[,,1:5], 1 , "mean"),
           Dynamic.split = apply(res.dndea.ds$eff.s2[,,1:5], 1 , "mean"),
           DS.win        = ifelse(apply(res.dndea.ss$eff.s2[,,1:5], 1 , "mean") <= apply(res.dndea.ds$eff.s2[,,1:5], 1 , "mean"), 1, 0),
           Diff          = apply(res.dndea.ds$eff.s2[,,1:5], 1 , "mean") - apply(res.dndea.ns$eff.s2[,,1:5], 1 , "mean"))

data.frame(Alpha = res.dndea.ds$alpha[,,1:3],
           Beta  = res.dndea.ds$beta[,,1:3],
           Gamma = res.dndea.ds$gamma[,,1:3])[,c(1, 4, 7, 2, 5, 8, 3, 6, 9)]
