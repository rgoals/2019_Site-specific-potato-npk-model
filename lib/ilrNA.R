ilrNA = function(comp, sbp, bal) {
  comp = unclass(comp)
  bal = unclass(bal)
  for (n in 1:nrow(comp)){
    for (p in 1:ncol(comp)) {
      for (q in 1:ncol(bal)) {
        if (sbp[q,p] != 0 & is.na(comp[n,p])) bal[n,q] <- NA
      }
    }
  }
  return(bal)
}