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

# Figure TBD
id.ex <- c(8, 14, 19, 20)
df.gg <- data.frame(Name     = rep(rownames(df.3d)[id.ex], each = 4),
                    Model    = rep(c("NS", "SS", "FS", "BS"), each = 4),
                    Eff.2014 = res.dndea.ns$eff.s2[id.ex,,1],
                    Eff.2015 = res.dndea.ns$eff.s2[id.ex,,2],
                    Eff.2016 = res.dndea.ns$eff.s2[id.ex,,3],
                    Eff.2017 = res.dndea.ns$eff.s2[id.ex,,4],
                    Eff.2018 = res.dndea.ns$eff.s2[id.ex,,5])

# KRICT
id.ex <- c(19)
beta  <- round(c(res.dndea.ss$beta [id.ex,,2], res.dndea.fs$beta [id.ex,,2], res.dndea.bs$beta [id.ex,,2]), 2)
gamma <- round(c(res.dndea.ss$gamma[id.ex,,1], res.dndea.fs$gamma[id.ex,,1], res.dndea.bs$gamma[id.ex,,1]), 2)
beta + gamma


# KIST and KRISS: BS.2016 < SS.2016 but BS.All > SS.All
id.ex <- c(14, 20)
data.frame(Name = table.5$Name[id.ex],
           Eff.2016.SS = round(res.dndea.ss$eff.s2[id.ex,,3], 4),
           Eff.2016.BS = round(res.dndea.bs$eff.s2[id.ex,,3], 4),
           Eff.All.SS  = round(apply(res.dndea.ss$eff.s2[id.ex,,], 1, "mean"), 4),
           Eff.All.BS  = round(apply(res.dndea.bs$eff.s2[id.ex,,], 1, "mean"), 4))


# KICT: FS.2016 > BS.2016 > NS.2016 > SS.2016
id.ex <- c(8)
res.fs.pro <- data.frame(Alpha = res.dndea.fs$alpha[id.ex,,3],
                         Beta  = res.dndea.fs$beta [id.ex,,2],
                         Gamma = res.dndea.fs$gamma[id.ex,,1])
res.bs.pro <- data.frame(Alpha = res.dndea.bs$alpha[id.ex,,3],
                         Beta  = res.dndea.bs$beta [id.ex,,2],
                         Gamma = res.dndea.bs$gamma[id.ex,,1])
data.frame(Name = table.5$Name[id.ex],
           Patent.actual = df.3d[id.ex, id.z, 3],
           Patent.effective.ss = sum(df.3d[id.ex, id.z, 1:3] * c(0.2,0.3,0.5)),
           Patent.effective.fs = sum(df.3d[id.ex, id.z, 1:3] * res.fs.pro),
           Patent.effective.bs = sum(df.3d[id.ex, id.z, 1:3] * res.bs.pro))


# KIT and KRIBB: no change
id.ex <- c(2, 15, 16, 18)
data.frame(Name = table.5$Name[id.ex],
           W.NS = res.dndea.ns$w.s1[id.ex,,3],
           W.SS = res.dndea.ss$w.s1[id.ex,,3],
           W.FS = res.dndea.fs$w.s1[id.ex,,3],
           W.BS = res.dndea.bs$w.s1[id.ex,,3])


# Further investigation
p.period <- 3
data.frame(Name          = rownames(df.3d),
           No.split      = res.dndea.ns$eff.s2[,,p.period],
           Fixed.split   = res.dndea.ss$eff.s2[,,p.period],
           Free.Split    = res.dndea.fs$eff.s2[,,p.period],
           Bounded.split = res.dndea.bs$eff.s2[,,p.period],
           NS.vs.BS      = ifelse(res.dndea.ns$eff.s2[,,p.period] <= res.dndea.bs$eff.s2[,,p.period], 1, 0),
           Diff.NS.SS    = round(res.dndea.ss$eff.s2[,,p.period] - res.dndea.ns$eff.s2[,,p.period], 4),
           Diff.NS.BS    = round(res.dndea.bs$eff.s2[,,p.period] - res.dndea.ns$eff.s2[,,p.period], 4),
           Diff.SS.BS    = round(res.dndea.bs$eff.s2[,,p.period] - res.dndea.ss$eff.s2[,,p.period], 4),
           Diff.FS.BS    = round(res.dndea.bs$eff.s2[,,p.period] - res.dndea.fs$eff.s2[,,p.period], 4))

