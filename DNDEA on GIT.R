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

# Table 2. Illustrative example
res.eff.ns.t1 <- dm.dea(df.toy[,2,1], df.toy[,3,1], orientation = "i")$eff
res.eff.ns.t2 <- dm.dea(df.toy[,2,2], df.toy[,3,2], orientation = "i")$eff
res.eff.ss.t1 <- dm.dea(df.toy[,2,1] * 0.7, df.toy[,3,1], orientation = "i")$eff
res.eff.ss.t2 <- dm.dea(df.toy[,2,1] * 0.3 + df.toy[,2,2], df.toy[,3,2], orientation = "i")$eff
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


# Table 4. Descriptive statistics
table.4 <- data.frame(Min  = apply(df.3d, 2, "min"),
                      Med  = apply(df.3d, 2, "median"),
                      Mean = round(apply(df.3d, 2, "mean")),
                      Max  = apply(df.3d, 2, "max"),
                      Std  = round(apply(df.3d, 2, "sd")))


# Table 5. Comparative results
res.dndea.ns <- dm.dynamic.network(df.3d[,id.x.s1,], df.3d[,id.y.s1,], df.3d[,id.z,], 
                                   df.3d[,id.x.s2,], df.3d[,id.y.s2,], rts)

res.dndea.ss <- dm.dynamic.network(df.3d[,id.x.s1,], df.3d[,id.y.s1,], df.3d[,id.z,], 
                                   df.3d[,id.x.s2,], df.3d[,id.y.s2,], rts, alpha = c(0.5, 0.3, 0.2))

res.dndea.fs <- dm.dynamic.network(df.3d[,id.x.s1,], df.3d[,id.y.s1,], df.3d[,id.z,],
                                   df.3d[,id.x.s2,], df.3d[,id.y.s2,], rts, alpha = "free", max.cp = 2)

res.dndea.bs <- dm.dynamic.network(df.3d[,id.x.s1,], df.3d[,id.y.s1,], df.3d[,id.z,], 
                                   df.3d[,id.x.s2,], df.3d[,id.y.s2,], rts, alpha = "free", max.cp = 2, LB = c(0.2, 0.2, 0.1))

table.5 <- data.frame(Name          = rownames(df.3d),
                      S1.overall    = apply(res.dndea.ns$eff.s1, 1, "mean"),
                      S1.2016       = res.dndea.ns$eff.s1[,,3],
                      S2.overall.NS = apply(res.dndea.ns$eff.s2, 1, "mean"),
                      S2.2016.NS    = res.dndea.ns$eff.s2[,,3],
                      S2.overall.SS = apply(res.dndea.ss$eff.s2, 1, "mean"),
                      S2.2016.SS    = res.dndea.ss$eff.s2[,,3],
                      S2.overall.FS = apply(res.dndea.fs$eff.s2, 1, "mean"),
                      S2.2016.FS    = res.dndea.fs$eff.s2[,,3],
                      S2.overall.BS = apply(res.dndea.bs$eff.s2, 1, "mean"),
                      S2.2016.BS    = res.dndea.bs$eff.s2[,,3])


# Discussions
# KIST and KRISS: BS.2016 < SS.2016 but BS.All > SS.All
id.ex <- c(14, 20)
data.frame(Eff.2016.SS = round(res.dndea.ss$eff.s2[id.ex,,3], 4),
           Eff.2016.BS = round(res.dndea.bs$eff.s2[id.ex,,3], 4),
           Eff.All.SS  = round(apply(res.dndea.ss$eff.s2[id.ex,,], 1, "mean"), 4),
           Eff.All.BS  = round(apply(res.dndea.bs$eff.s2[id.ex,,], 1, "mean"), 4))


# KICT: FS.2016 > BS.2016 > NS.2016 > SS.2016
id.ex      <- c(8)
data.frame(Eff.2016.NS = round(res.dndea.ss$eff.s2[id.ex,,3], 4),
           Eff.2016.FS = round(res.dndea.fs$eff.s2[id.ex,,3], 4),
           Eff.2016.BS = round(res.dndea.bs$eff.s2[id.ex,,3], 4))
res.fs.pro <- data.frame(Alpha = res.dndea.fs$alpha[id.ex,,1:3],
                         Beta  = res.dndea.fs$beta [id.ex,,1:3],
                         Gamma = res.dndea.fs$gamma[id.ex,,1:3])
res.bs.pro <- data.frame(Alpha = res.dndea.bs$alpha[id.ex,,1:3],
                         Beta  = res.dndea.bs$beta [id.ex,,1:3],
                         Gamma = res.dndea.bs$gamma[id.ex,,1:3])
data.frame(Patent.actual = df.3d[id.ex, id.z, 3],
           Patent.acc.ss = df.3d[id.ex, id.z, 1:3] %*% c(0.2,0.3,0.5),
           Patent.acc.fs = df.3d[id.ex, id.z, 1:3] %*% c(res.fs.pro[1, 3], res.fs.pro[2, 2], res.fs.pro[3, 1]),
           Patent.acc.bs = df.3d[id.ex, id.z, 1:3] %*% c(res.bs.pro[1, 3], res.bs.pro[2, 2], res.bs.pro[3, 1]))


data.frame(Eff.2016.NS = round(table.5$S2.2016.NS, 4),
           Eff.2016.SS = round(table.5$S2.2016.SS, 4),
           Win         = round(table.5$S2.2016.SS - table.5$S2.2016.NS, 4))

data.frame(Eff.2016.NS = round(table.5$S2.2016.NS, 4),
           Eff.2016.FS = round(table.5$S2.2016.FS, 4),
           Win         = round(table.5$S2.2016.FS - table.5$S2.2016.NS, 4))

data.frame(Eff.2016.NS = round(table.5$S2.2016.NS, 4),
           Eff.2016.BS = round(table.5$S2.2016.BS, 4),
           Win         = round(table.5$S2.2016.BS - table.5$S2.2016.NS, 4))

data.frame(Eff.2016.FS = round(table.5$S2.2016.FS, 4),
           Eff.2016.BS = round(table.5$S2.2016.BS, 4),
           Win         = round(table.5$S2.2016.FS - table.5$S2.2016.BS, 4))




# Further investigation
p.period <- 3
data.frame(No.split      = res.dndea.ns$eff.s2[,,p.period],
           Static.split  = res.dndea.ss$eff.s2[,,p.period],
           Free.Split    = res.dndea.fs$eff.s2[,,p.period],
           Bounded.split = res.dndea.bs$eff.s2[,,p.period],
           NS.vs.BS      = ifelse(res.dndea.ns$eff.s2[,,p.period] <= res.dndea.bs$eff.s2[,,p.period], 1, 0),
           Diff.NS.BS    = round(res.dndea.bs$eff.s2[,,p.period] - res.dndea.ns$eff.s2[,,p.period], 2),
           Diff.SS.BS    = round(res.dndea.bs$eff.s2[,,p.period] - res.dndea.ss$eff.s2[,,p.period], 2),
           Diff.FS.BS    = round(res.dndea.bs$eff.s2[,,p.period] - res.dndea.fs$eff.s2[,,p.period], 2))

data.frame(Alpha = res.dndea.fs$alpha[,,1:p.period],
           Beta  = res.dndea.fs$beta[,,1:p.period],
           Gamma = res.dndea.fs$gamma[,,1:p.period])

data.frame(Alpha = res.dndea.bs$alpha[,,1:p.period],
           Beta  = res.dndea.bs$beta[,,1:p.period],
           Gamma = res.dndea.bs$gamma[,,1:p.period])

data.frame(Alpha = res.dndea.ss$alpha[,,1:p.period],
           Beta  = res.dndea.ss$beta[,,1:p.period],
           Gamma = res.dndea.ss$gamma[,,1:p.period])


###############
# Discussions #
###############

# DMU 19: (0,0.44,0) FS & (0.15, 0.36, 0.2) BS & (0.2, 0.3, 0.5) SS > NS - How time lag consideration changes the results

# DMU 16, 18: No differences - zero weight attached to z