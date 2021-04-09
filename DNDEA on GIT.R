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
load("DNDEA.14.18.Rdata")


#########################################################################################################################
### Analysis
#########################################################################################################################

# Table 2. Illustrative Example
res.eff.ns.t1 <- dm.dea(df.toy[,2,1], df.toy[,3,1], orientation = "i")$eff
res.eff.ns.t2 <- dm.dea(df.toy[,2,2], df.toy[,3,2], orientation = "i")$eff
res.eff.ss.t1 <- dm.dea(df.toy[,2,1]*0.7, df.toy[,3,1], orientation = "i")$eff
res.eff.ss.t2 <- dm.dea(df.toy[,2,1]*0.3 + df.toy[,2,2], df.toy[,3,2], orientation = "i")$eff
res.eff.ds.t1 <- dm.dea(df.toy[,2,1] + c(0, -5), df.toy[,3,1], orientation = "i")$eff
res.eff.ds.t2 <- dm.dea(df.toy[,2,2] + c(0, +5), df.toy[,3,2], orientation = "i")$eff

table.2 <- data.frame(df.toy[,,1],
                      E_no.t1 = sprintf("%.2f", res.eff.ns.t1,3),
                      E_ss.t1 = sprintf("%.2f", res.eff.ss.t1,3),
                      E_ds.t1 = sprintf("%.2f", res.eff.ds.t1,3),
                      df.toy[,,2],
                      E_no.t2 = sprintf("%.2f", res.eff.ns.t2,3),
                      E_ss.t2 = sprintf("%.2f", res.eff.ss.t2,3),
                      E_ds.t2 = sprintf("%.2f", res.eff.ds.t2,3))


# Table 3. Descriptive Statistics
table.3 <- data.frame(Min  = apply(df.3d, 2, "min"),
                      Med  = apply(df.3d, 2, "median"),
                      Mean = round(apply(df.3d, 2, "mean")),
                      Max  = apply(df.3d, 2, "max"),
                      Std  = round(apply(df.3d, 2, "sd")))


# DNDEA
res.dndea.ns <- dm.dynamic.network(df.3d[,id.x.s1,], df.3d[,id.y.s1,], df.3d[,id.z,], 
                                   df.3d[,id.x.s2,], df.3d[,id.y.s2,], rts)

res.dndea.ss <- dm.dynamic.network(df.3d[,id.x.s1,], df.3d[,id.y.s1,], df.3d[,id.z,], 
                                   df.3d[,id.x.s2,], df.3d[,id.y.s2,], rts, alpha = c(0.4, 0.4, 0.2))

res.dndea.ds <- dm.dynamic.network(df.3d[,id.x.s1,], df.3d[,id.y.s1,], df.3d[,id.z,], 
                                     df.3d[,id.x.s2,], df.3d[,id.y.s2,], rts, alpha = "free", max.cp = 2)

res.dndea.ds.b <- dm.dynamic.network(df.3d[,id.x.s1,], df.3d[,id.y.s1,], df.3d[,id.z,], 
                                   df.3d[,id.x.s2,], df.3d[,id.y.s2,], rts, alpha = "free", max.cp = 2, LB = c(0.2,0.2,0.1))

# Comparison of composite indicies
res.ns.c.index <- (res.dndea.ns$eff.s1 + res.dndea.ns$eff.s2) * 0.5
res.ss.c.index <- (res.dndea.ss$eff.s1 + res.dndea.ss$eff.s2) * 0.5
res.ds.c.index <- (res.dndea.ds$eff.s1 + res.dndea.ds$eff.s2) * 0.5
data.frame(No.Split      = apply(res.ns.c.index[,,1:5], 1, "mean"),
           Static.Split  = apply(res.ss.c.index[,,1:5], 1, "mean"),
           Dynamic.Split = apply(res.ds.c.index[,,1:5], 1, "mean"))

data.frame(No.split      = res.dndea.ns$eff.s2[,,3],
           Static.split  = res.dndea.ss$eff.s2[,,3],
           Dynamic.split = res.dndea.ds$eff.s2[,,3],
           Diff          = round(res.dndea.ds$eff.s2[,,3] - res.dndea.ns$eff.s2[,,3], 2),
           DS.win        = ifelse(res.dndea.ss$eff.s2[,,3] <= res.dndea.ds$eff.s2[,,3], 1, 0))

data.frame(No.split      = apply(res.dndea.ns$eff.s2[,,1:5], 1, "mean"),
           Static.split  = apply(res.dndea.ss$eff.s2[,,1:5], 1, "mean"),
           Dynamic.split = apply(res.dndea.ds$eff.s2[,,1:5], 1, "mean"),
           DS.win        = ifelse(apply(res.dndea.ss$eff.s2[,,1:5], 1, "mean") <= apply(res.dndea.ds$eff.s2[,,1:5], 1, "mean"), 1, 0),
           Diff          = round(apply(res.dndea.ds$eff.s2[,,1:5], 1, "mean") - apply(res.dndea.ns$eff.s2[,,1:5], 1, "mean"), 2))

data.frame(Alpha = res.dndea.ds$alpha[,,1:3],
           Beta  = res.dndea.ds$beta[,,1:3],
           Gamma = res.dndea.ds$gamma[,,1:3])[,c(1, 4, 7, 2, 5, 8, 3, 6, 9)]

