# Helper functions for permutation testing

# Permuting individual columns of a count matrix on the log-ratio scale and 
# transforming back to original count scale.
# createPerm <- function(x) {
#   rs <- rowSums(x)
#   lr <- compositions::alr(x)
#   
#   # Permutating on log-ratio scale
#   nr <- NROW(lr)
#   nc <- NCOL(lr)
#   for (i in 1:nc) {
#     lr[,i] <- lr[sample(nr),i]
#   }
#   
#   # Rescaling proportions back to original count scale by row
#   lr.inv <- unclass(compositions::alrInv(lr))
#   x.new <- floor(lr.inv*rs)
#   
#   return(x.new)
# }
createPerm <- function(x) {
  rs <- rowSums(x)
  lr <- compositions::alr(x)
  
  # Permutating on log-ratio scale
  nr <- NROW(lr)
  nc <- NCOL(lr)
  for (i in 1:nr) {
    lr[i,] <- lr[i,sample(nc)]
  }
  
  # Rescaling proportions back to original count scale by row
  lr.inv <- unclass(compositions::alrInv(lr))
  x.new <- floor(lr.inv*rs)
  
  return(x.new)
}

# Run naive variation/rho/phi on zero-imputed data
# Requires "variation_functions.R"
# If pc=NULL then zero-imputation happens
# IF pc is a non-negative value, then that value is added as a pseudo-count
runNaive <- function(x, pc=NULL, lr=c("alr", "clr")) {
  if (any(x==0)) {
    if (is.null(pc)){
      x <- as.matrix(zCompositions::cmultRepl(x, output = "p-counts", z.warning = 0.9999,
                                              suppress.print=TRUE))
    } else {
      x <- x+pc
    }
  }
  v <- naiveVariation(x, lr=lr)
  phi <- naiveVariation(x, type="phi", lr=lr)
  phis <- naiveVariation(x, type="phis", lr=lr)
  rho <- naiveVariation(x, type="rho", lr=lr)
  return(list(v=v,phi=phi,phis=phis,rho=rho,imp.data=x))
}
# Run model estimated variation/rho/phi
# Requires "variation_functions.R" and "mle_multinom_logit_normal.R"
runModelBased <- function(x, impute=FALSE, n.lam=8,
                          lambda.min.ratio=0.1, 
                          lambda.gl=NULL,
                          lr=c("alr","clr")) {
  if (any(x==0) & impute) {
    x <- as.matrix(zCompositions::cmultRepl(x, output = "p-counts",
                                            z.warning = 0.9999,
                                            suppress.print=TRUE))
  }
  
  if (is.null(lambda.gl)) {
    mle <- mlePath(x, tol = 1e-4, tol.nr = 1e-4,
                   n.lambda = n.lam, 
                   lambda.min.ratio = lambda.min.ratio, 
                   n.cores = 8,
                   lr.penalty=lr)
  } else {
    mle <- mlePath(x, tol = 1e-4, tol.nr = 1e-4,
                   lambda.gl = lambda.gl, 
                   n.cores = 8,
                   lr.penalty=lr)
  }
  mu.hat <- mle$est.min$mu
  Sigma.hat <- mle$est.min$Sigma
  
  v <- logitNormalVariation(mu.hat, Sigma.hat, lr=lr)
  phi <- logitNormalVariation(mu.hat, Sigma.hat, type = "phi", lr=lr)
  phis <- logitNormalVariation(mu.hat, Sigma.hat, type = "phis", lr=lr)
  rho <- logitNormalVariation(mu.hat, Sigma.hat, type = "rho", lr=lr)
  
  return(list(v=v,phi=phi, phis=phis, rho=rho, 
              mu=mu.hat, Sigma=Sigma.hat, mle=mle))
}




