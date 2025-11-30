# MLE for multinomial logistic normal dist using EM algorithm from Hoff 2003

# Function to get MLE
mleLR <- function(y, max.iter=10000, max.iter.nr=100, tol=1e-6, tol.nr=1e-6, lambda.gl=0, gamma=0.1) {
  n <- NROW(y)
  k <- NCOL(y)
  ni <- rowSums(y)
  # observed proportions
  #p.hat <- y/rowSums(y)
  # initialize 
  pseudo.count <- 0.1
  v <- unclass(compositions::alr(y+pseudo.count))
  attr(v, "orig") <- NULL
  attr(v, "V") <- NULL
  colnames(v) <- NULL
  mu <- rep(0, k-1)
  mu.old <- mu+1
  
  # Sigma <- cov(v)
  # Sigma.old <- Sigma+0.1
  # i.inv <- array(0, dim=c(k-1,k-1,n))
  
  Sigma.inv <- diag(k-1) #MASS::ginv(cov(v))
  Sigma.inv.old <- Sigma.inv+0.1
  i.inv <- array(0, dim=c(k-1,k-1,n))
  
  count <- 1
  #while(sum(abs((mu-mu.old)/mu.old))>tol & sum(abs((diag(Sigma.inv)-diag(Sigma.inv.old))/diag(Sigma.inv.old)))>tol) {
  while(max(abs((mu-mu.old)/mu.old))>tol & max(abs((diag(Sigma.inv)-diag(Sigma.inv.old))/diag(Sigma.inv.old)))>tol) { 
    if (count%%10==0) {
      cat("Iter=", count, "\n")
      cat("Error=", max(abs((mu-mu.old)/mu.old)), "\n")
      #print(Sigma[1,1])
      print(mu)
      #print(eigen(h)$values)
    }
    if (count>max.iter) stop("Maximum number of iterations reached")
    
    #Sigma.inv <- qr.solve(Sigma)
    
    for (i in 1:n) {
      # cat("i=", i, "\n")
      g <- rep(1, k-1)
      count.nr <- 1
      while(abs(sum(g))>tol.nr) {
        if (count.nr>max.iter.nr) {
          print(g)
          stop("Maximum number of Newton-Raphson iterations reached")
        }
        g <- grad(v[i,], y[i,-k], ni[i], mu, Sigma.inv)
        h <- hess(v[i,], ni[i], Sigma.inv)
        #if (any(is.complex(eigen(h)$values))) browser()
        i.inv[,,i] <- qr.solve(-h)
        # Newton-Raphson step
        v[i,] <- v[i,] + i.inv[,,i]%*%g
        count.nr <- count.nr + 1
      }
    }
    
    mu.old <- mu
    mu <- colMeans(v)
    #Sigma.old <- Sigma
    Sigma.inv.old <- Sigma.inv
    #Sigma.inv <- huge::huge(v, lambda.gl, method="glasso", verbose=FALSE)$icov[[1]]
    gl <- suppressWarnings(glasso::glasso(cov(v), lambda.gl))
    Sigma.inv <- gl$wi
    Sigma <- gl$w
    # Sigma <- matrix(0, k-1, k-1)
    # for (i in 1:n){
    #   # Need to update Hessian to final value of v[i,]
    #   h <- hess(v[i,], ni[i], Sigma.inv)
    #   neg.hess.inv <- qr.solve(-h)
    #   Sigma <- Sigma + (tcrossprod(v[i,]-mu) + neg.hess.inv)/n
    # }
    count <- count+1
    
  }
  
  S <- cov(v)
  
  log.lik <- logLik(v, y, ni, S, Sigma.inv)
  # Currently only using Gaussian part for loglik
  #log.lik <- logLikG(v, S, Sigma.inv)
  
  df <- sum(Sigma.inv[lower.tri(Sigma.inv)]!=0)
  eb <- ebic(log.lik, n, k-1, df, gamma)
  
  # return(list(v=v, mu=mu, Sigma=Sigma, var.v=neg.hess.inv))
  return(list(v=v, mu=mu, Sigma.inv=Sigma.inv, Sigma=Sigma, log.lik=log.lik, 
              ebic=eb, df=df))
}

# Wrapper for mleLR() to help with mclapply
wrapMLE <- function(x) {
  do.call(mleLR, x)
}

# Computes the path for various values of GLasso in mleLR()
mlePath <- function(y, max.iter=10000, max.iter.nr=100, tol=1e-6, tol.nr=1e-6, lambda.gl=NULL,
                    lambda.min.ratio=0.1, n.lambda=1,
                    n.cores=NULL, gamma=0.1) {

k <- NCOL(y)

# Set lambda.gl if it's not specified by the user
  if (is.null(lambda.gl)) {
    pc <- 0.1
    alr.y <- log(y[,-k]+pc) - log(y[,k]+pc)
    d <- NCOL(alr.y)
    S <- cov(alr.y)
    lmax <- max(max(S - diag(d)), -min(S - diag(d)))
    lmin <- lambda.min.ratio*lmax
    lambda.gl <- exp(seq(log(lmin), log(lmax), length.out=n.lambda))
    #lambda.gl <- c(0.0968419,0.1869724,0.3609873,0.6969573,1.3456139,2.5979737,5.0159019,9.6841903)
  }
  
  n.lam <- length(lambda.gl)
  
  m.pars <- vector("list", length=n.lam)
  for (i in 1:n.lam) {
    m.pars[[i]] <- list(y, max.iter, max.iter.nr, tol, tol.nr, lambda.gl[i], gamma)
  }
  
  if (is.null(n.cores)) parallel::detectCores()

  est <- parallel::mclapply(m.pars, wrapMLE, mc.cores = n.cores)

  print("inside_here2")

  ebic.vec <- unlist(lapply(est, function(x){x$ebic}))
  wm <- which.min(ebic.vec)
  
  est.min <- est[[wm]]
  
  return(list(est=est, est.min=est.min,
              lambda.gl=lambda.gl, min.idx=wm, ebic=ebic.vec))
}

# Function to calculate gradient with respect to normal RVs on logit scale
grad <- function(v, y, ni, mu, Sigma.inv) {
  ev <- exp(v)
  c(y - ni*ev/(1+sum(ev)) - Sigma.inv%*%(v-mu))
}

# Function to calculate hessian
hess <- function(v, ni, Sigma.inv) {
  ev <- exp(v)
  p.hat <- ev/(1+sum(ev))
  ni*(tcrossprod(p.hat)-diag(p.hat)) - Sigma.inv
}

# Log-likelihood just for the Gaussian part
logLikG <- function(v, S, invSigma) {
  n <- NROW(v)
  ldet <- determinant(invSigma, logarithm = TRUE)$modulus
  n*0.5*(ldet - sum(diag(S%*%invSigma)))
}

logLik <- function(v, y, ni, S, invSigma) {
  n <- NROW(y)
  n.sp <- NCOL(y)
  rs <- log(rowSums(exp(v)+1))
  ldet <- determinant(invSigma, logarithm = TRUE)$modulus
  sum(y[,-n.sp]*v) - sum(ni*rs) + n*0.5*(ldet - sum(diag(S%*%invSigma)))
}

ebic <- function(l, n, d, df, gamma) {
  -2*l + log(n)*df + 4*gamma*log(d)*df
}

# EBIC plot for 
ebicPlot <- function(fit, xlog=FALSE) {
  if (xlog) {
    x.ax <- log(fit$lambda.gl)
  } else {
    x.ax <- fit$lambda.gl
  }
  plot(x.ax, fit$ebic, type="o", pch=19, cex=0.5,
       col="darkred", xlab="log(lambda)", ylab="EBIC")
}

rmse <- function(x,y) {
  sqrt(mean((x-y)^2))
}

rmse_by_row <- function(x,y) {
  n <- NROW(x)
  r <- rep(0, n)
  for (i in 1:n) {
    r[i] <- rmse(x[i,], y[i,])
  }
  return(r)
}

