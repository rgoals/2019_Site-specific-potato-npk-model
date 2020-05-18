CoDaDendrogram2 = function(comp, V, equal.height=FALSE, range=c(-4,4),
         show.range=TRUE, n_digits=1,
         group=NULL,
         type="none", conf.level=0.95, conf.method="t", 
         leaf_las = 2,
         pch.col = TRUE, ...) 
######
# A fork of the CoDaDendrogram function in the compositions R package (van den Boogart et al., 2013)
# CoDaDendrogram2 allows to plot statistics on the vertical branches.
# comp: compositional matrix 
# V: orthonormal basis
# equal.heigh: boolean, if TRUE, vertival branches will all have the same height. 
#              if FALSE, the height of vertical bars is the variance of the balance divided by total variance
# range: a vector of 2 numbers, designing the range of all balances. 
#  For automatic range scaled on all data, specify "auto"
#  For automatic range scaled on confidence intervals, specify "autoconf"
# show.range: boolean, show range labels on nodes
# n_digits: the number of digits to show in range
# group: grouping factor for type != "none"
# type: plots over balances, choices "boxplot", "conf" and "density"
# conf.level: confidence level if type is "conf"
# conf.method: if type is "conf", "t" for student (qt) or "normal" for normal (qnorm)
#####

{
  # funct
  ## initialisation
  V <- V[ , order(colSums(V == 0))]
  sbp = t(sign(V))
  signary = data.frame(t(sbp))
  Vzero = abs(V) > 1e-15 # non-zero values in the balance
  nb_parts_in_bal = colSums(abs(signary)) # number of parts implied in each balance
  ord_o = order(nb_parts_in_bal, decreasing = TRUE)
  V_o = V[, ord_o] # reordering the orthog. mat. per nb of parts. eq. to Vo$ilrBase
  comp_id = do.call("order", as.list(signary)) # comp_id: eq. to Vo$order # order(signary[,1], signary[,2], ...)
  bal_o = ilr(comp, V=V_o) # eq. to idtx
  var_bal_o = apply(bal_o,2,var) # variance per ordered bal
  mean_bal_o = apply(bal_o,2,mean) # mean per ordered bal
  signary_o = sign(V_o) # ordered signary (signary is the transposed sbp)
  binary_o = abs(signary_o) # ordered binary (binary is the absolute of the signary)
  
  ## y position of each fulcrum
  if(equal.height) {
    heights = rep(1, ncol(binary_o))
    maxheight = 1
    for (i in 1:ncol(binary_o)) {
      heights[i:ncol(binary_o)] = heights[i:ncol(binary_o)] - sign(binary_o[, i] %*% binary_o[, i:ncol(binary_o)])
    }
    inv_heights = max(abs(heights))+heights + 1 # +1 to elevate lowest from one step, in order to see all leave branches
    scaled_heights = inv_heights/max(inv_heights)
    heights = scaled_heights
    dh = rep(heights[1] - heights[2], times=length(heights))
  } else {
    heights = rep(0, ncol(binary_o))
    heights[1:length(heights)] = max(binary_o[, 1:ncol(binary_o)] %*% var_bal_o[1:ncol(binary_o)])
    maxheight = max(heights)
    for (i in 1:ncol(binary_o)) {
      heights[i:ncol(binary_o)] = heights[i:ncol(binary_o)] - var_bal_o[i] * sign(binary_o[, i] %*% binary_o[, i:ncol(binary_o)])
    }
    dh = var_bal_o
  }  
  
  # horizontal bars' ranges
  if (is.matrix(range)) {
    range=range # obviously
  } else if (is.numeric(range)) { # not neat. change condition here. if any character of whatever not numeric, it supposes auto
    range = as.matrix(range)
    range = range[,rep(1,ncol(bal_o))]
    auto.range = FALSE
  } else if (is.character(range)) {
    if (range=='auto') {
      range = matrix()
      range = apply(bal_o, 2, range)
      auto.range = TRUE
    } else if (range=='autoconf') {
      range_min = rep(NA, ncol(bal_o))
      range_max = range_min
      for (i in 1:nlevels(group)){
        confrange = apply(bal_o[group==levels(group)[i], ], 2, ci, level=conf.level)
        for (j in 1:ncol(bal_o)) {
          range_min[j] = min(range_min[j], confrange[1, j], na.rm=TRUE)
          range_max[j] = max(range_max[j], confrange[3, j], na.rm=TRUE)
        }
      }
      range = rbind(range_min, range_max)
      auto.range = TRUE
    } else {
      range = matrix()
      range = apply(bal_o, 2, range)
      auto.range = TRUE
      print('unknown range type, using auto')
    }
  } else {
    range = matrix()
    range = apply(bal_o, 2, range)
    auto.range = TRUE
    print('unknown range type, using auto')
  }
  
  
  ## x-position of the leaves (parts)
  nodes = matrix(0, ncol = 2, nrow = ncol(binary_o))
  is.leaf = nodes
  for (i in ncol(binary_o):1) {
    nparts = c(sum(signary_o[, i] < 0), sum(signary_o[, i] > 0)) # nb of parts with -1 sign, nb of parts with 1 sign in the ith balance
    for (j in 1:2) {
      if (nparts[j] == 1) { # if there is only one part on the jth side of the ith balance
        take = signary_o[, i] == (-1)^j # which one is it
        take = c(1:length(comp_id))[take] == comp_id # to which comp_id does it correspond?
        nodes[i, j] = c(1:length(comp_id))[take]
        is.leaf[i, j] = 1 # indicate with 1 that the node is a leave (no more sub-balance)
      }
      else { # if there is more than one part on this side of the ith balance
        take = signary_o[, i] == (-1)^j # which one are they?
        take = as.integer(take) * binary_o[, i] # ??? seems to repeat first take...
        take = as.logical(take)
        aux = rep(1, sum(take)) %*% binary_o[take, ]
        aux[1:i] = 0
        i1 = c(1:length(aux))[aux == max(aux)]
        if(!auto.range) nodes[i, j] = approx(x = range[,i], y = nodes[i1,], xout = mean_bal_o[i1], rule = 2)$y
        if(auto.range) nodes[i, j] = approx(x = range[,i1], y = nodes[i1,], xout = mean_bal_o[i1], rule = 2)$y
      }
    }
  }
  
  ## plot
  ### window
  if (equal.height) {
    xlim = c(1, nrow(V))
    ylim = c(0, maxheight+dh[1])
  } else {
    xlim = c(0, nrow(V))
    ylim = c(0, maxheight)
  }
  
  plot.window(xlim = xlim, ylim = ylim)
  if(equal.height) {
    plot(x = xlim, y = ylim, ann = FALSE, bty = "n",
         xaxt = "n", yaxt = "n", type="n")
  } else {
    plot(x = xlim, y = ylim, ann = FALSE, bty = "n", xaxt = "n", type = "n")
  }
  
  ### leaf labels
  axis(side = 1, at = 1:ncol(comp), labels = colnames(comp)[comp_id], las = leaf_las, lty = 0) 
  
  ### horizontal bars
  for (i in 1:(ncol(V))) {
    segments(x0 = nodes[i, 1], y0 = heights[i], x1 = nodes[i, 2], y1 = heights[i])
  }
  
  ### liens vers les feuilles
  for (i in 1:nrow(is.leaf)) {
    for (j in 1:ncol(is.leaf)) {
      if (is.leaf[i, j]) {
        segments(x0 = nodes[i, j], y0 = heights[i], x1 = nodes[i, j], y1 = 0)
      }
    }
  }
  
  ### range of each ilr
  if (show.range) {
    for (i in 1:(ncol(V))) {
      text(x = c(nodes[i, 1], nodes[i, 2]), y = heights[i],
           labels = formatC(c(round(range[1,i], n_digits), round(range[2,i], n_digits)), format='f', digits=n_digits),
           pos=3, cex=0.6)
    }
  }
  
  ### vertical lines
  ### they are placed on mean, but there should be an option to place them on target of median of whatever...
  for (i in 1:length(var_bal_o)) {
    aux = approx(x = range[,i], y = nodes[i, ], xout = mean_bal_o[i], rule = 2)$y
    segments(x0 = aux, y0 = heights[i], x1 = aux, y1 = heights[i] + dh[i])
  }
  
  ### boxplots
  if (type == "boxplot") {
    if (is.null(group)) {
      group = factor(rep(1, times=nrow(comp)))
      ngroup = 1
    } else {
      ngroup = length(levels(group))
    }
    
    dh_gr = matrix(ncol=ngroup, nrow=length(heights))
    for(i in 1:length(heights)) {
      dh_gr[i,] = (ngroup:1)*dh[i]/(ngroup+1)
    }
    
    for (i in 1:length(heights)) {
      approxFunc <- approxfun(x = range[,i], y = nodes[i, ], rule = 2)
      for (j in 1:ngroup) {
        bal_type = bal_o[group == levels(group)[j],i]
        bp = boxplot.stats(bal_type)
        quantile_mapped = approxFunc(bp$stats)
        out_mapped = approxFunc(bp$out)
        segments(x0=quantile_mapped[1], y0=heights[i]+dh_gr[i,j], 
                 x1=quantile_mapped[5], y1=heights[i]+dh_gr[i,j])
        polygon(x=c(quantile_mapped[2], quantile_mapped[4], quantile_mapped[4], quantile_mapped[2]),
                y=c(heights[i]+dh_gr[i,j]-dh[i]/(ngroup+1)/2.5, heights[i]+dh_gr[i,j]-dh[i]/(ngroup+1)/2.5,
                    heights[i]+dh_gr[i,j]+dh[i]/(ngroup+1)/2.5, heights[i]+dh_gr[i,j]+dh[i]/(ngroup+1)/2.5),
                col=j+1)
        segments(x0=quantile_mapped[3], y0=heights[i]+dh_gr[i,j]-dh[i]/(ngroup+1)/2.5, 
                 x1=quantile_mapped[3], y1=heights[i]+dh_gr[i,j]+dh[i]/(ngroup+1)/2.5,lwd=2)
        points(out_mapped, y=rep(heights[i]+dh_gr[i,j], times=length(out_mapped)), cex=1.5/ngroup)
      }
    }
  }
  
  ### confidence intervals
  if (type == "conf") {
    if (is.null(group)) {
      group = factor(rep(1, times=nrow(comp)))
      ngroup = 1
    } else {
      ngroup = length(levels(group))
    }
    
    dh_gr = matrix(ncol=ngroup, nrow=length(heights))
    for(i in 1:length(heights)) {
      dh_gr[i,] = (ngroup:1)*dh[i]/(ngroup+1)
    }
    
    for (i in 1:length(heights)) {
      approxFunc <- approxfun(x = range[,i], y = nodes[i, ], rule = 2)
      for (j in 1:ngroup) {
        bal_type = bal_o[group == levels(group)[j],i]
        conf = ci(bal_type, level=conf.level, method=conf.method)
        conf_mapped = approxFunc(conf)
        segments(x0=conf_mapped[1], y0=heights[i]+dh_gr[i,j], 
                 x1=conf_mapped[3], y1=heights[i]+dh_gr[i,j])
        segments(x0=conf_mapped[1], y0=heights[i]+dh_gr[i,j]+dh[i]/(ngroup+1)/3, 
                 x1=conf_mapped[1], y1=heights[i]+dh_gr[i,j]-dh[i]/(ngroup+1)/3)
        segments(x0=conf_mapped[3], y0=heights[i]+dh_gr[i,j]+dh[i]/(ngroup+1)/3, 
                 x1=conf_mapped[3], y1=heights[i]+dh_gr[i,j]-dh[i]/(ngroup+1)/3)
        if(pch.col) {
          points(x=conf_mapped[2], y=heights[i]+dh_gr[i,j], pch=21, bg=j+1, cex=0.7)
        } else {
          points(x=conf_mapped[2], y=heights[i]+dh_gr[i,j], pch=21, bg='grey35', cex=0.7)
        }
      }
    }
  }
  
  ### density
  if (type == "density") {
    
    if (is.null(group)) {
      group = factor(rep(1, times=nrow(comp)))
      ngroup = 1
    } else {
      ngroup = length(levels(group))
    }
    
    for (i in 1:length(heights)) {
      
      # compute densities
      n=100 # 100 is the number of points in the density
      dens = list()
      ymax=0
      for (j in 1:ngroup) {
        bal_type = bal_o[group == levels(group)[j],i]
        tmp = density(bal_type, n=n)
        dens[[j]] = data.frame(x=tmp$x, y=tmp$y)
        if(ymax < max(dens[[j]]$y)) ymax=max(dens[[j]]$y)
      }
      
      approxFunc_x <- approxfun(x = range[,i], y = nodes[i, ], rule = 2)
      approxFunc_y <- approxfun(x = c(0, ymax*1.1), y = c(heights[i], heights[i]+dh[i]), rule = 2) # *1.1 to avoid that the density top coincide to the upper horizontal bar
      
      for (j in 1:ngroup) {
        dens_mapped = dens[[j]]
        dens_mapped$x = approxFunc_x(dens[[j]]$x)
        dens_mapped$y = approxFunc_y(dens[[j]]$y)
        lines(dens_mapped, col=j+1)
      }
    }
  }
  
}

# home made conf. int. function
ci <- function (x, level=conf.level, method="t") {
  a <- mean(x)
  s <- sd(x)
  n <- length(x)
  
  if (method == "t") {
    error <- qt(1-(1-level)/2, df=n-1)*s/sqrt(n)
  } else {
    error <- qnorm(1-(1-level)/2)*s/sqrt(n)
  }
  
  left <- a-error
  right <- a+error
  out <- c(left, a, right)
  return(out)
}
