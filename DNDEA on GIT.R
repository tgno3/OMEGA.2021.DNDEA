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
t.w <- c(0.1,0.2,0.4,0.2,0.1)
res.dndea.ns <- dm.dynamic.network(df.3d[,id.x.s1,], df.3d[,id.y.s1,], df.3d[,id.z,], 
                                   df.3d[,id.x.s2,], df.3d[,id.y.s2,], rts, t.w = t.w)

res.dndea.ss <- dm.dynamic.network(df.3d[,id.x.s1,], df.3d[,id.y.s1,], df.3d[,id.z,], 
                                   df.3d[,id.x.s2,], df.3d[,id.y.s2,], rts, alpha = c(0.4, 0.4, 0.2), t.w = t.w)

res.dndea.ds <- dm.dynamic.network(df.3d[,id.x.s1,], df.3d[,id.y.s1,], df.3d[,id.z,],
                                   df.3d[,id.x.s2,], df.3d[,id.y.s2,], rts, alpha = "free", max.cp = 2, t.w = t.w)

res.dndea.bs <- dm.dynamic.network(df.3d[,id.x.s1,], df.3d[,id.y.s1,], df.3d[,id.z,], 
                                   df.3d[,id.x.s2,], df.3d[,id.y.s2,], rts, alpha = "free", max.cp = 2, LB = c(0.2,0.2,0.0), t.w = t.w)

# Comparison of composite indicies
s.period <- 1:5
res.ns.c.index <- (res.dndea.ns$eff.s1 + res.dndea.ns$eff.s2) * 0.5
res.ss.c.index <- (res.dndea.ss$eff.s1 + res.dndea.ss$eff.s2) * 0.5
res.ds.c.index <- (res.dndea.ds$eff.s1 + res.dndea.ds$eff.s2) * 0.5
res.bs.c.index <- (res.dndea.bs$eff.s1 + res.dndea.bs$eff.s2) * 0.5
# Composite Index
data.frame(No.Split      = apply(res.ns.c.index[,,s.period], 1, "mean"),
           Static.Split  = apply(res.ss.c.index[,,s.period], 1, "mean"),
           Dynamic.Split = apply(res.ds.c.index[,,s.period], 1, "mean"),
           Bounded.Split = apply(res.bs.c.index[,,s.period], 1, "mean"),
           Diff.NS.BS    = round(apply(res.bs.c.index[,,s.period], 1, "mean") - apply(res.ns.c.index[,,s.period], 1, "mean"), 2),
           Diff.SS.BS    = round(apply(res.bs.c.index[,,s.period], 1, "mean") - apply(res.ss.c.index[,,s.period], 1, "mean"), 2),
           Diff.DS.BS    = round(apply(res.bs.c.index[,,s.period], 1, "mean") - apply(res.ds.c.index[,,s.period], 1, "mean"), 2))

# Stage 2 Eff
data.frame(No.split      = apply(res.dndea.ns$eff.s2[,,s.period], 1, "mean"),
           Static.split  = apply(res.dndea.ss$eff.s2[,,s.period], 1, "mean"),
           Dynamic.split = apply(res.dndea.ds$eff.s2[,,s.period], 1, "mean"),
           Bounded.split = apply(res.dndea.bs$eff.s2[,,s.period], 1, "mean"),
           NS.vs.DS      = ifelse(apply(res.dndea.ns$eff.s2[,,s.period], 1, "mean") <= apply(res.dndea.ds$eff.s2[,,s.period], 1, "mean"), 1, 0),
           Diff.NS.BS    = round(apply(res.dndea.bs$eff.s2[,,s.period], 1, "mean") - apply(res.dndea.ns$eff.s2[,,s.period], 1, "mean"), 2),
           Diff.SS.BS    = round(apply(res.dndea.bs$eff.s2[,,s.period], 1, "mean") - apply(res.dndea.ss$eff.s2[,,s.period], 1, "mean"), 2),
           Diff.DS.BS    = round(apply(res.dndea.bs$eff.s2[,,s.period], 1, "mean") - apply(res.dndea.ds$eff.s2[,,s.period], 1, "mean"), 2))

s.period <- 3
data.frame(No.split      = res.dndea.ns$eff.s2[,,s.period],
           Static.split  = res.dndea.ss$eff.s2[,,s.period],
           Dynamic.split = res.dndea.ds$eff.s2[,,s.period],
           Bounded.split = res.dndea.bs$eff.s2[,,s.period],
           NS.vs.BS      = ifelse(res.dndea.ns$eff.s2[,,s.period] <= res.dndea.bs$eff.s2[,,s.period], 1, 0),
           Diff.NS.BS    = round(res.dndea.bs$eff.s2[,,s.period] - res.dndea.ns$eff.s2[,,s.period], 2),
           Diff.SS.BS    = round(res.dndea.bs$eff.s2[,,s.period] - res.dndea.ss$eff.s2[,,s.period], 2),
           Diff.DS.BS    = round(res.dndea.bs$eff.s2[,,s.period] - res.dndea.ds$eff.s2[,,s.period], 2))

data.frame(Alpha = res.dndea.bs$alpha[,,s.period],
           Beta  = res.dndea.bs$beta[,,s.period],
           Gamma = res.dndea.bs$gamma[,,s.period])

