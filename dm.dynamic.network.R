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
### DNDEA - max.cp of 2
#########################################################################################################################

dm.dynamic.network <- function(xdata.s1, ydata.s1 = NULL, zdata, xdata.s2 = NULL, ydata.s2, 
                               rts = "crs", orientation = "i", type = "nc", leader = "1st", 
                               alpha = 1, max.cp = 0, t.w = NULL, o = NULL, exhaust = T){
  
  # Initial checks
  if(is.na(match(rts, c("crs", "vrs", "irs", "drs"))))          stop('rts must be "crs", "vrs", "irs", or "drs".')
  if(is.na(match(orientation, c("i", "o"))))                    stop('orientation must be either "i" or "o".')
  if(is.na(match(type, c("co", "nc"))))                         stop('type must be "co" or "nc".')
  if(is.na(match(leader, c("1st", "2nd"))))                     stop('leader must be "1st" or "2nd".')
  if(!is.null(t.w) && dim(xdata.s1)[3] != length(t.w))          stop('t.w must have a length of the period.')
  if(!is.null(o) && !all(o <= nrow(xdata.s1)))                  stop('o must be element(s) of n.')
  if(!isTRUE(alpha == "free") && sum(alpha) != 1)               stop('alpha must sum to one.')

  # Parameters
  xdata.s1 <- if(length(dim(xdata.s1)) != 3) array(xdata.s1, c(dim(xdata.s1)[1], 1, dim(xdata.s1)[2])) else as.array(xdata.s1)
  ydata.s2 <- if(length(dim(ydata.s2)) != 3) array(ydata.s2, c(dim(ydata.s2)[1], 1, dim(ydata.s2)[2])) else as.array(ydata.s2)
  zdata    <- if(length(dim(zdata   )) != 3) array(zdata,    c(dim(zdata)[1],    1, dim(zdata)[2]))    else as.array(zdata)
  if(!is.null(ydata.s1)){ydata.s1 <- if(length(dim(ydata.s1)) != 3) array(ydata.s1, c(dim(ydata.s1)[1], 1, dim(ydata.s1)[2])) else as.array(ydata.s1)}
  if(!is.null(xdata.s2)){xdata.s2 <- if(length(dim(xdata.s2)) != 3) array(xdata.s2, c(dim(xdata.s2)[1], 1, dim(xdata.s2)[2])) else as.array(xdata.s2)}
  n        <- nrow(xdata.s1)
  m.s1     <- ncol(xdata.s1)
  s.s1     <- ifelse(is.null(ydata.s1), 0, ncol(ydata.s1))
  p        <- ncol(zdata)
  m.s2     <- ifelse(is.null(xdata.s2), 0, ncol(xdata.s2))
  s.s2     <- ncol(ydata.s2)
  o        <- if(is.null(o))  c(1:n) else as.vector(o)
  #alpha    <- if(max.cp == 0) c(1)   else alpha
  max.cp   <- ifelse(is.numeric(alpha), length(alpha) - 1, max.cp)
  t        <- dim(xdata.s1)[3]
  t.w      <- if(is.null(t.w)) rep(1, t)
  
  # Data frames
  res.eff.sys <- array(NA, c(n,    1, t))
  res.eff.s1  <- array(NA, c(n,    1, t))
  res.eff.s2  <- array(NA, c(n,    1, t))
  res.v.s1    <- array(NA, c(n, m.s1, t)) 
  res.u.s1    <- array(NA, c(n, s.s1, t)) 
  res.w.s1    <- array(NA, c(n,    1, t)) 
  res.v.s2    <- array(NA, c(n, m.s2, t)) 
  res.u.s2    <- array(NA, c(n, s.s2, t)) 
  res.w.s2    <- array(NA, c(n,    1, t))
  res.p       <- array(NA, c(n,    p, t)) 
  res.a       <- array(NA, c(n,    p, t)) 
  res.b       <- array(NA, c(n,    p, t))
  res.c       <- array(NA, c(n,    p, t)) 
  
  # Pre-emptive phase for non-cooperative model
  if(type == "nc"){
    for(q in 1:t){res.eff.s1[,,q] <- dm.dea(xdata.s1[,,q], cbind(ydata.s1[,,q], zdata[,,q]), rts, orientation)$eff}  
  }
  
  # Main phase
  # Indices for coding convenience
  no.dv.t <- (m.s1 + s.s1 + p + p + p + p + 1 + m.s2 + s.s2 + 1) 
  id.v.s1 <- 1:m.s1 + rep(no.dv.t * 0:(t-1), each = m.s1)
  id.u.s1 <- if(is.null(ydata.s1)) 0 else id.v.s1 + m.s1
  id.p    <- (m.s1 + s.s1 + 1):(m.s1 + s.s1 + p) + rep(no.dv.t * 0:(t-1), each = p)
  id.a    <- (m.s1 + s.s1 + p + 1):(m.s1 + s.s1 + p + p) + rep(no.dv.t * 0:(t-1), each = p)
  id.b    <- (m.s1 + s.s1 + p + p + 1):(m.s1 + s.s1 + p + p + p) + rep(no.dv.t * 0:(t-1), each = p)
  id.c    <- (m.s1 + s.s1 + p + p + p + 1):(m.s1 + s.s1 + p + p + p + p) + rep(no.dv.t * 0:(t-1), each = p)
  id.w.s1 <-  m.s1 + s.s1 + p + p + p + p + 1 + rep(no.dv.t * 0:(t-1), each = p)
  id.v.s2 <- if(is.null(xdata.s2)) 0 else (m.s1 + s.s1 + p + p + p + p + 2):(m.s1 + s.s1 + p + p + p + p + 1 + m.s2) + rep(no.dv.t * 0:(t-1), each = m.s2) 
  id.u.s2 <- (m.s1 + s.s1 + p + p + p + p + 1 + m.s2 + 1):(m.s1 + s.s1 + p + p + p + p + 1 + m.s2 + s.s2) + rep(no.dv.t * 0:(t-1), each = s.s2)
  id.w.s2 <-  m.s1 + s.s1 + p + p + p + p + 1 + m.s2 + s.s2 + 1 + rep(no.dv.t * 0:(t-1), each = p)
  id.z    <- z.ag <- list()
  for(i in 1:t){
    if(i == 1){
      id.z[[i]] <- c(id.a[i])
      z.ag[[i]] <- matrix(zdata[,,1], n)
    }else if(i == 2){
      id.z[[i]] <- c(id.a[i], id.b[i - 1])
      z.ag[[i]] <- cbind(zdata[,,i], zdata[,,i - 1])
    }else{
      id.z[[i]] <- c(id.a[i], id.b[i - 1], id.c[i - 2])
      z.ag[[i]] <- cbind(zdata[,,i], zdata[,,i - 1], zdata[,,i - 2])
    }
  }

  # LP  
  for(k in o){
    # Declare LP
    lp.s2 <- make.lp(0, no.dv.t * t) # v (+ u) + p + a + b + c + w (+ v) + u + w
    
    # Labelling
    temp <- c()
    for(i in 1:t){temp <- c(temp, 
                            paste0(i, "v1", 1:m.s1), 
                            if(!is.null(ydata.s1)) paste0(i, "u1", 1:s.s1), 
                            paste0(i, "p", 1:p), paste0(i, "a", 1:p), paste0(i, "b", 1:p), paste0(i, "c", 1:p), 
                            paste0(i, "w1"), 
                            if(!is.null(xdata.s2)) paste0(i, "v2", 1:m.s2), 
                            paste0(i, "u2", 1:s.s2), paste0(i, "w2"))}
    dimnames(lp.s2)[[2]] <- temp
    
    # Set objective
    if(orientation == "i"){
      if(type == "co"){
        set.objfn(lp.s2, c(-ydata.s2[k,,] * t.w, rep(t.w, 2)), indices = c(id.u.s2, id.w.s1, id.w.s2))
      }else{
        set.objfn(lp.s2, c(-ydata.s2[k,,] * t.w, t.w), indices = c(id.u.s2, id.w.s2))
      }
    }else{
      # Obj for output orientation: TBA
    }
    
    # RTS
    if(rts == "crs") add.constraint(lp.s2, rep(1, 2*t), indices = c(id.w.s1, id.w.s2), "=", 0)
    #if(rts == "irs") add.constraint(lp.s2, c(rep(1, n)), indices = c(1:n), ">=", 1)
    #if(rts == "drs") add.constraint(lp.s2, c(rep(1, n)), indices = c(1:n), "<=", 1)
    
    # Constraints
    for(i in 1:t){
      if(orientation == "i"){
        #############
        ## Stage 1 ##
        #############
        
        # # v1x1 = 1 for DMU o
        # add.constraint(lp.s2, xdata.s1[k,,i],
        #                indices = id.v.s1[(i*m.s1 - (m.s1 - 1)):(i*m.s1)], "=", 1)
        # 
        # # (u1y1 +) pz - w = eff.s1 for DMU o
        # if(type == "nc"){
        #   if(is.null(ydata.s1)){
        #     add.constraint(lp.s2, c(zdata[k,,i], -1),
        #                    indices = c(id.p[(i*p - (p - 1)):(i*p)],
        #                                id.w.s1[i]), "=", res.eff.s1[k,,i])
        #   }else{
        #     add.constraint(lp.s2, c(ydata.s1[k,,i], zdata[k,,i], -1),
        #                    indices = c(id.u.s1[(i*m.s1 - (m.s1 - 1)):(i*m.s1)],
        #                                id.p[(i*p - (p - 1)):(i*p)],
        #                                id.w.s1[i]), "=", res.eff.s1[k,,i])
        #   }
        # }

        # (u1y1 +) pz - w - eff.s1 * v1x1 = 0 for DMU o
        if(type == "nc"){
          if(is.null(ydata.s1)){
            add.constraint(lp.s2, c(-xdata.s1[k,,i] * res.eff.s1[k,,i], zdata[k,,i], -1),
                           indices = c(id.v.s1[(i*m.s1 - (m.s1 - 1)):(i*m.s1)],
                                       id.p[(i*p - (p - 1)):(i*p)],
                                       id.w.s1[i]), "=", 0)        
          }else{
            add.constraint(lp.s2, c(-xdata.s1[k,,i] * res.eff.s1[k,,i], ydata.s1[k,,i], zdata[k,,i], -1),
                           indices = c(id.v.s1[(i*m.s1 - (m.s1 - 1)):(i*m.s1)],
                                       id.u.s1[(i*m.s1 - (m.s1 - 1)):(i*m.s1)], 
                                       id.p[(i*p - (p - 1)):(i*p)],
                                       id.w.s1[i]), "=", 0)            
          }
        }
        
        # -v1x1 (+ u1y1) + pz - w <= 0 for all
        for(d in o){
          if(is.null(ydata.s1)){
            add.constraint(lp.s2, c(-xdata.s1[d,,i], zdata[d,,i], -1), 
                           indices = c(id.v.s1[(i*m.s1 - (m.s1 - 1)):(i*m.s1)], 
                                       id.p[(i*p - (p - 1)):(i*p)],
                                       id.w.s1[i]), "<=", 0)
          }else{
            add.constraint(lp.s2, c(-xdata.s1[d,,i], ydata.s1[d,,i], zdata[d,,i], -1), 
                           indices = c(id.v.s1[(i*m.s1 - (m.s1 - 1)):(i*m.s1)], 
                                       id.u.s1[(i*m.s1 - (m.s1 - 1)):(i*m.s1)], 
                                       id.p[(i*p - (p - 1)):(i*p)],
                                       id.w.s1[i]), "<=", 0)  
          }
        }
        
        #############
        ## Stage 2 ##
        #############
        # abcz (+ v2x2) = 1 for DMU o when there exist external input
        if(type == "nc" & !(is.null(xdata.s2))){
          add.constraint(lp.s2, c(z.ag[[i]][k,], xdata.s2[k,,i]), 
                         indices = c(id.z[[i]],
                                     id.v.s2[(i*m.s2 - (m.s2 - 1)):(i*m.s2)]), "=", 1)  
        }
        
        # (-v2x2 +) + u2y2 - abcz - w <= 0 for all
        for(d in o){
          if(is.null(xdata.s2)){
            add.constraint(lp.s2, c(ydata.s2[d,,i], -z.ag[[i]][d, ], -1), 
                           indices = c(id.u.s2[(i*s.s2 - (s.s2 - 1)):(i*s.s2)], 
                                       id.z[[i]],
                                       id.w.s2[i]), "<=", 0)
          }else{
            add.constraint(lp.s2, c(-xdata.s2[d,,i], ydata.s2[d,,i], -z.ag[[i]][d, ], -1), 
                           indices = c(id.v.s2[(i*m.s2 - (m.s2 - 1)):(i*m.s2)], 
                                       id.u.s2[(i*s.s2 - (s.s2 - 1)):(i*s.s2)], 
                                       id.z[[i]],
                                       id.w.s2[i]), "<=", 0)  
          }
        }
        
        # Treatments for carry-over
        if(isTRUE(alpha == "free")){
          for(j in 1:p){
            # p = a + b + c for all periods
            add.constraint(lp.s2, c(1, -1, -1, -1), 
                                       indices = c(id.p[(i-1)*p + j], id.a[(i-1)*p + j], id.b[(i-1)*p + j], id.c[(i-1)*p + j]), "=", 0)
            
            # b & c == 0 when no carry-over
            if(max.cp == 0){add.constraint(lp.s2, rep(1, length(c(id.b, id.c))), indices = c(id.b, id.c), "=", 0)}
            
            # c = 0 for the period second to the last when exhaust is T
            if(isTRUE(exhaust) & i == (t - 1)){add.constraint(lp.s2, c(1), indices = c(id.c[(i-1)*p + j]), "=", 0)}
            
            # p = a for the last period when exhaust is T
            if(isTRUE(exhaust) & i ==  t     ){add.constraint(lp.s2, c(1, -1), indices = c(id.p[(i-1)*p + j], id.a[(i-1)*p + j]), "=", 0)}  
          }
        }else if(is.numeric(alpha)){
          for(j in 1:p){
            # p = a + b + c for all periods
            add.constraint(lp.s2, c(1, -1, -1, -1), 
                           indices = c(id.p[(i-1)*p + j], id.a[(i-1)*p + j], id.b[(i-1)*p + j], id.c[(i-1)*p + j]), "=", 0)
            for(a in 1:(max.cp + 1)){
              if(a == 1     ){
                if(isTRUE(exhaust) & i == t){
                  # p = a for the last period when exhaust is T
                  add.constraint(lp.s2, c(1, -1), indices = c(id.p[(i-1)*p + j], id.a[(i-1)*p + j]), "=", 0)
                }else{
                  # p * alpha[1] = a
                  add.constraint(lp.s2, c(alpha[a], -1), indices = c(id.p[(i-1)*p + j], id.a[(i-1)*p + j]), "=", 0)  
                }
              }
              else if(a == 2){
                if(isTRUE(exhaust) & i == (t - 1) & max.cp == 2){
                  # p * (alpha[2] + alpha[3]) = b when exhaust is T and max.cp is 2
                  add.constraint(lp.s2, c(alpha[a] + alpha[a + 1], -1), indices = c(id.p[(i-1)*p + j], id.b[(i-1)*p + j]), "=", 0)
                }else if(i < t){
                  # p * alpha[2] = b
                  add.constraint(lp.s2, c(alpha[a], -1), indices = c(id.p[(i-1)*p + j], id.b[(i-1)*p + j]), "=", 0)
                }
              }
              else if(a == 3){
                if(i < (t - 1)){
                  # p * alpha[3] = c
                  add.constraint(lp.s2, c(alpha[a], -1), indices = c(id.p[(i-1)*p + j], id.c[(i-1)*p + j]), "=", 0)
                }  
              }
            }      
          }
        }

      }else{ 
        # Constraints for output-orientation: TBA
      }
    }
    
    # Bounds
    if(rts == "vrs"){
      temp <- rep(0, no.dv.t * t)
      temp[c(id.w.s1, id.w.s2)] <- -Inf
      set.bounds(lp.s2, lower = temp)  
    }
    
    # Solve
    solve.lpExtPtr(lp.s2)
    
    # Plan B
    if(solve.lpExtPtr(lp.s2) == 3){
      for(i in 1:t){
        # v1x1 = 1 for DMU o
        add.constraint(lp.s2, xdata.s1[k,,i],
                       indices = id.v.s1[(i*m.s1 - (m.s1 - 1)):(i*m.s1)], "=", 1)
      }
      
      # Re-solve
      solve.lpExtPtr(lp.s2)
    }
    
    # Get results
    res.all            <- get.variables(lp.s2)
    if(type == "co"){
      res.temp.pz      <- colSums(matrix(res.all[id.p], ncol = t) * zdata[k,,])
      res.temp.uy      <- colSums(matrix(res.all[id.u.s2], ncol = t) * ydata.s2[k,,])
      res.temp.w1      <- res.all[id.w.s1]
      res.temp.w2      <- res.all[id.w.s2]
      res.eff.s1[k,,]  <- res.temp.pz - res.temp.w1
      res.eff.s2[k,,]  <- (res.temp.uy - res.temp.w2) / res.temp.pz
      res.eff.sys[k,,] <- res.temp.uy - res.temp.w1 - res.temp.w2
    }else{
      res.temp.pzw     <- colSums(matrix(res.all[id.u.s2], ncol = t) * ydata.s2[k,,]) - res.all[id.w.s2]
      res.eff.s2[k,,]  <- if(is.null(xdata.s2) & alpha[1] == 1) res.temp.pzw/res.eff.s1[k,,] else res.temp.pzw
      res.eff.sys[k,,] <- res.eff.s1[k,,] * res.eff.s2[k,,]
    }
    res.v.s1[k,,]      <- res.all[id.v.s1]
    res.u.s1[k,,]      <- if(is.null(ydata.s1)) NULL else res.all[id.u.s1]
    res.w.s1[k,,]      <- res.all[id.w.s1]
    res.v.s2[k,,]      <- if(is.null(xdata.s2)) NULL else res.all[id.v.s2]
    res.u.s2[k,,]      <- res.all[id.u.s2]
    res.w.s2[k,,]      <- res.all[id.w.s2]
    res.p[k,,]         <- res.all[id.p]
    res.a[k,,]         <- res.all[id.a]
    res.b[k,,]         <- res.all[id.b]
    res.c[k,,]         <- res.all[id.c]
  }
  results <- list(eff.sys = res.eff.sys, eff.s1 = res.eff.s1, eff.s2 = res.eff.s2, 
                  v.s1 = res.v.s1, u.s1 = res.u.s1, w.s1 = res.w.s1,
                  v.s2 = res.v.s2, u.s2 = res.u.s2, w.s2 = res.w.s2, 
                  p = res.p, a = res.a, b = res.b, c = res.c, alpha = res.a/res.p, beta = res.b/res.p, gamma = res.c/res.p)
  return(results)
}


#########################################################################################################################
### Analysis
#########################################################################################################################

# Descriptive Stats
table.1 <- data.frame(Min  = unlist(aggregate(df.2d[, c(id.x, id.y) + 3], by = list(df.2d$Location), "min")),
                      Med  = unlist(aggregate(df.2d[, c(id.x, id.y) + 3], by = list(df.2d$Location), "median")),
                      Mean = unlist(aggregate(df.2d[, c(id.x, id.y) + 3], by = list(df.2d$Location), "mean")),
                      Max  = unlist(aggregate(df.2d[, c(id.x, id.y) + 3], by = list(df.2d$Location), "max")),
                      Std  = unlist(aggregate(df.2d[, c(id.x, id.y) + 3], by = list(df.2d$Location), "sd")))

print(cbind(Loc  = rep(levels(loc), 3), format(round(table.1[-c(1:3),]), big.mark = ",")))

# Toy model
id.toy  <- c(6, 16, 18, 19, 20)
res.toy <- dm.dynamic.network(df.3d[id.toy,id.x.s1,1:3], df.3d[id.toy,id.y.s1,1:3], df.3d[id.toy,id.z,1:3, drop = F], 
                              df.3d[id.toy,id.x.s2,1:3], df.3d[id.toy,id.y.s2,1:3,drop = F], rts, ori, alpha = "free")


# DNDEA
res.100  <- dm.dynamic.network(df.3d[,id.x.s1,], df.3d[,id.y.s1,], df.3d[,id.z,, drop = F], 
                                     df.3d[,id.x.s2,], df.3d[,id.y.s2,,drop = F], rts, ori, alpha = c(1, 0, 0))

res.532  <- dm.dynamic.network(df.3d[,id.x.s1,], df.3d[,id.y.s1,], df.3d[,id.z,, drop = F], 
                                     df.3d[,id.x.s2,], df.3d[,id.y.s2,,drop = F], rts, ori, alpha = c(0.5, 0.3, 0.2))

res.free <- dm.dynamic.network(df.3d[,id.x.s1,], df.3d[,id.y.s1,], df.3d[,id.z,, drop = F], 
                                     df.3d[,id.x.s2,], df.3d[,id.y.s2,,drop = F], rts, ori, alpha = "free")

data.frame(No.Split = res.100$eff.sys, Static.Split = res.532$eff.sys, Free.Split = res.free$eff.sys)

