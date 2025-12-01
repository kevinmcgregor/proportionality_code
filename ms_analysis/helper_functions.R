# Helper functions for MS data analysis

getTopVals <- function(v1, v2, cut.val=100) {
  lt <- lower.tri(v1)
  abs.v1 <- abs(v1[lt])
  abs.v2 <- abs(v2[lt])
  sort.v1 <- sort(abs.v1, decreasing = TRUE)
  sort.v2 <- sort(abs.v2, decreasing = TRUE)
  # Cutoff value for each
  cutoff.v1 <- sort.v1[cut.val]
  cutoff.v2 <- sort.v2[cut.val]
  wh.top.v1 <- which(abs.v1>=cutoff.v1)
  wh.top.v2 <- which(abs.v2>=cutoff.v2)
  
  l <- list(top=list(MS=wh.top.v1, Control=wh.top.v2),
            cutoff1=cutoff.v1, cutoff2=cutoff.v2)
  return(l)
}
