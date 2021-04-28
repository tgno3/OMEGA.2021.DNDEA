#########################################################################################################################
### Project  : Dynamic Network DEA
### Script   : DNDEA on GIT.R
### Contents : Measuring dynamic efficiency with variable time lag effects
#########################################################################################################################

#########################################################################################################################
### Setting up environment
#########################################################################################################################

# Load library
pkgs <- c("DJL", "ggplot2", "showtext")
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
# No-split scenario
res.dndea.ns <- dm.dynamic.network(df.3d[,id.x.s1,], df.3d[,id.y.s1,], df.3d[,id.z,], 
                                   df.3d[,id.x.s2,], df.3d[,id.y.s2,], rts)
# Fixed-split scenario
res.dndea.ss <- dm.dynamic.network(df.3d[,id.x.s1,], df.3d[,id.y.s1,], df.3d[,id.z,], 
                                   df.3d[,id.x.s2,], df.3d[,id.y.s2,], rts, alpha = c(0.5, 0.3, 0.2))
# Free-split scenario
res.dndea.fs <- dm.dynamic.network(df.3d[,id.x.s1,], df.3d[,id.y.s1,], df.3d[,id.z,],
                                   df.3d[,id.x.s2,], df.3d[,id.y.s2,], rts, alpha = "free", max.cp = 2)
# Bounded-split scenario
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


# KAERI, KISTI, KIT and KRIBB: no change
id.ex <- c(2, 15, 16, 18)
data.frame(Name = table.5$Name[id.ex],
           W.NS = res.dndea.ns$w.s1[id.ex,,3],
           W.SS = res.dndea.ss$w.s1[id.ex,,3],
           W.FS = res.dndea.fs$w.s1[id.ex,,3],
           W.BS = res.dndea.bs$w.s1[id.ex,,3])


# Figure 4. Efficiency changes
id.ex     <- c(3, 5, 8, 19, 20, 22)
id.p      <- 4
{
  id.dmu    <- id.ex[id.p]
  scenario  <- c("No-split", "Fixed-split", "Free-split", "Bounded-split")
  df.p.per  <- df.3d[,id.z,] / rowSums(df.3d[,id.z,]) * 100
  df.r.per  <- df.3d[,id.y.s2,] / rowSums(df.3d[,id.y.s2,]) * 100
  df.eff.s2 <- array(c(res.dndea.ns$eff.s2, res.dndea.ss$eff.s2, res.dndea.fs$eff.s2, res.dndea.bs$eff.s2), 
                     c(nrow(df.3d), 5, 4), dimnames = list(rownames(df.3d), 2014:2018, scenario))
  df.y1.ex  <- data.frame(Year       = rep(2014:2018, 4),
                          Scenario   = rep(scenario, each = 5),
                          Efficiency = c((df.eff.s2[id.dmu,,])))
  df.y2.ex  <- data.frame(Year       = rep(2014:2018, 2),
                          Type       = rep(c("Patent", "Royalty"), each = 5),
                          Proportion = c(df.p.per[id.dmu,], df.r.per[id.dmu,]))
  df.scale  <- data.frame(DMU = id.ex,
                          MUL = c(120 ,100, 75, 75, 100, 100),
                          MIS = c(36, 40, 45, 45, 30, 50))
  font_add('Cambria', 'Cambria.ttf'); showtext_auto()
  ggplot(df.y2.ex) +
    geom_col (df.y2.ex, mapping = aes(x = Year, y = Proportion, fill = Type), position = 'dodge', width = 0.5) +
    geom_line(df.y1.ex, mapping = aes(x = Year, y = Efficiency * df.scale[id.p, 2] - df.scale[id.p, 3], 
                                      color = Scenario), size = 1.6) +
    scale_y_continuous(limits   = c(0, 30), breaks = seq(0, 30, 10), 
                       labels   = paste0(seq(0, 30, 10), "%"), guide = guide_axis(position = 'right'), 
                       sec.axis = sec_axis(~./df.scale[id.p, 2] + df.scale[id.p, 3]/df.scale[id.p, 2], 
                                           guide = guide_axis(position = 'left'))) + 
    scale_fill_manual (values = c("Patent" = "#e2f0d9", "Royalty" = "#fbe5d6")) +
    scale_color_manual(values = c("Bounded-split" = "#0070c0", "Fixed-split" = "#7c7c7c", "Free-split" = "#f9312c", "No-split" = "#0c0b0a")) +
    ggtitle(rownames(df.3d)[id.dmu]) +
    theme(legend.position = "bottom", 
          legend.background = element_rect(fill = "white"),
          legend.key = element_rect(fill = "white"),
          legend.title = element_blank(),
          legend.text = element_text(family = "Cambria"),
          axis.ticks = element_blank(),
          axis.text = element_text(family = "Cambria", size = 20),
          axis.title.y.left = element_blank(),
          axis.title.y.right = element_blank(),
          axis.title.x = element_blank(),
          plot.title = element_text(family = "Cambria",face = "bold", hjust = 0.5, size = 30),
          panel.background = element_rect(fill = "white", colour = "white"),
          panel.grid.major.y = element_line(color = "grey"))
  
}

