#' Function to calculate variation matrix under various models for lattice compositional data.
#'
#' @param counts Write what to put in counts here
#'
#' @return
#' @export
#'
#' @examples
varEst <- function(counts, p.model=c("logitNormal", "dirichlet", "plugin"), type=c("standard","phi","rho"), refColumn=NULL) {
  if (any(counts<0) | any(counts!=floor(counts)) | !is.matrix(counts)) stop("'counts' must be a matrix containing non-negative integers")
  p.model <- match.arg(p.model)
  type <- match.arg(type)
  
  if (is.null(refColumn)&p.model=="logitNormal") { 
    refColumn <- J
  }
  
  if (p.model=="logitNormal") {
    result <- logitNormalVariation(mu, Sigma, lmu, lsigma, type=type)
  } else if (p.model=="dirichlet") {
    result <- dirichletVariation(counts)
  } else {
    result <- pluginVariation(counts)
  }
  
  return(result)
}


logitNormalVariation <- function(mu, Sigma, type=c("standard","phi", "phis","rho"), 
                                 lr=c("alr", "clr")) {
  type <- match.arg(type)
  lr <- match.arg(lr)

  #if (length(lr)>1) stop("lr must be of length 1")
  #if (!lr%in%c("alr","clr")) stop("lr must be either 'alr' or 'clr'")
  
  J <- length(mu)
  ones <- rep(1, J)
  d.S <- diag(Sigma)
  V <- tcrossprod(d.S, ones) + tcrossprod(ones, d.S) - 2*Sigma
  V <- cbind(rbind(V, d.S), c(d.S, 0))
  dimnames(V) <- NULL
  
  if (lr=="alr") {
    if (type=="phi") {
      den <- tcrossprod(d.S, ones)
      V[1:J,1:J] <- V[1:J,1:J]/den
      V[J+1,1:J] <- rep(Inf, J)
      V[,J+1] <- c(rep(1, J), NaN)
    } else if (type=="phis" | type=="rho") {
      # Calculate rho.  If type is "phis" then rho is still calculated
      # but gets transformed to phis later.
      V <- cbind(rbind(2*Sigma/outer(d.S, d.S, "+"), 0), 0)
      V[J+1,J+1] <- 1
    }
  } else {
    H.inv <- qr.solve(diag(J) + matrix(1, J, J))
    Fm <- cbind(diag(J), -1)
    HiF <- H.inv%*%Fm
    Sigma.clr <- crossprod(HiF, Sigma)%*%HiF
    d.Sc <- diag(Sigma.clr)
    if (type=="phi") {
      V <- V/tcrossprod(d.Sc, c(ones, 1))
    } else if (type=="phis" | type=="rho") {
      # Calculate rho.  If type is "phis" then rho is still calculated
      # but gets transformed to phis later.
      V <- 2*Sigma.clr/outer(d.Sc, d.Sc, "+")
    }
  }
  
  if (type=="phis") {
    V <- (1-V)/(1+V)
    if (lr=="alr") {
      V[J+1,1:J] <- rep(1, J)
      V[,J+1] <- c(rep(1, J), 0)
    }
  }
  
  return(V)
}


logitNormalVariationOld <- function(mu, Sigma, type=c("standard","phi", "phis","rho"),
                                 order=c("first", "second")) {
  type <- match.arg(type)
  J <- length(mu)
  
  print(order)

  ones <- rep(1, J)
  d.S <- diag(Sigma)
  V <- tcrossprod(d.S, ones) + tcrossprod(ones, d.S) - 2*Sigma

  if (type=="phi") {
    lv <- logVarTaylor(mu, Sigma, order=order)
    lv.row <- tcrossprod(diag(lv)[-(J+1)], ones)
    V <- V/lv.row
  } else if (type=="phis") {
    # lv <- logVarTaylor(mu, Sigma, order=order)
    # lv.d <- diag(lv)
    # den <- outer(lv.d, lv.d, "+") + l3v
    # V <- V/den
  } else if (type=="rho") {
    lv <- logVarTaylor(mu, Sigma, order=order)
    lv.d <- diag(lv)
    den <- outer(lv.d, lv.d, "+")
    V <- 2*lv/den
  }
   
  return(V)
}

# logitNormalVariationOLD_OLD <- function(mu, Sigma, lmu, lsigma, type=c("standard","phi","rho")) {
#   type <- match.arg(type)
#   J <- length(mu)
# 
#   d.S <- diag(Sigma)
#   ones <- rep(1, J)
#   
#   tcp.dS <- tcrossprod(d.S, ones)
#   Elr <- tcp.dS + t(tcp.dS) - 2*Sigma
#   En <- exp(-lmu+lsigma^2/2)
#   
#   tcp.e.mu.dS <- tcrossprod(exp(-mu+d.S/2), ones)
#   sum.rows <- rowSums(tcp.e.mu.dS+tcp.dS/2-Sigma)
#   
#   Ep <- tcp.e.mu.dS*(1+tcrossprod(sum.rows, ones))
# 
#   V <- Elr + En*(2+Ep+t(Ep))
#   diag(V) <- 0
#   return(V)
# }

dirichletVariation <- function(counts) {
  return(NULL)
}

pluginVariation <- function(counts) {
  return(NULL)
}

#' Title
#'
#' @param counts 
#' @param pseudo.count 
#' @param type 
#' @param use 
#' @param set.inf.na 
#' @param already.log 
#'
#' @return
#' @export
#'
#' @examples
naiveVariation <- function(counts, pseudo.count=0, type=c("standard","phi", "phis", "rho"),
                           lr=c("alr", "clr"), impute.zeros=TRUE, ...) {
  type <- match.arg(type)
  lr <- match.arg(lr)

  if (!is.matrix(counts) | !is.numeric(counts)) stop("counts must be a numeric matrix")
  if (!is.logical(impute.zeros)) stop("impute.zeros must be TRUE or FALSE")
  if (!is.numeric(pseudo.count)) stop("pseudo.count must be numeric")
  if (pseudo.count<0) stop("pseudo.count must be non-negative")
  
  l <- counts
  
  l <- l + pseudo.count
  if (lr=="alr") {
    l <- compositions::alr(l)
  } else {
    l <- compositions::clr(l)
  }

  J <- NCOL(l)
    
  v <- matrix(0,J,J)
  for (i in 1:J) {
    for (j in 1:J){
      if (type=="standard") {
        v[i,j] <- var(l[,i]-l[,j])
      } else if (type=="phi") {
        v[i,j] <- var(l[,i]-l[,j])/var(l[,i])
      } else if (type=="phis") {
        v[i,j] <- var(l[,i]-l[,j])/(var(l[,i]+l[,j]))
      } else if (type=="rho") {
        v[i,j] <- 2*cov(l[,i],l[,j])/(var(l[,i])+var(l[,j]))
      }
    }
  }
  
  return(v)
}

# Function to get MCMC samples from logit-normal model
# K is the number of samples to return. Number of features returned
# is length(mu)+1
# Can also sample from logit-t dist if lr.dist=t.  df is the degrees of freedom
MCSample <- function(mu, Sigma, K=1, df=NULL, lr.dist=c("normal", "t", "tnc")) {
  lr.dist <- match.arg(lr.dist)
  
  if (lr.dist=="normal") {
    x.norm <- mvtnorm::rmvnorm(K, mu, Sigma)
  } else if (lr.dist=="t") {
    x.norm <- tcrossprod(rep(1, K), mu) + 
                mvtnorm::rmvt(K, (df-2)/df*Sigma, df)
  } else {
    nc.adj <- sqrt(df/2)*gamma((df-1)/2)/gamma(df/2)
    x.norm <- mvtnorm::rmvt(K, (df-2)/df*Sigma, df,
                            mu/nc.adj, type="Kshirsagar")
  }
  
  ex <- cbind(exp(x.norm), 1)
  return(ex/rowSums(ex))
}

# Calculating variation using Monte Carlo integration
MCVariation <- function(mu=NULL, Sigma=NULL, x=NULL, K=1e6, 
                        type=c("standard","phi", "phis","rho","logx")) {
  
  if (is.null(mu) & is.null(Sigma) & is.null(x)) {
    stop("If x is missing then mu and Sigma must both not be missing.")
  }
  
  type <- match.arg(type)
  
  if (is.null(x)) {
    x <- MCSample(mu, Sigma, K=K)
  }
  
  l <- log(x)
  J <- NCOL(x)
  
  v <- matrix(0,J,J)
  for (i in 1:J) {
    for (j in 1:J){
      if (type=="standard") {
        v[i,j] <- var(l[,i]-l[,j])
      } else if (type=="phi") {
        v[i,j] <- var(l[,i]-l[,j])/var(l[,i])
      } else if (type=="phis") {
        v[i,j] <- var(l[,i]-l[,j])/(var(l[,i]+l[,j]))
      } else if (type=="rho") {
        v[i,j] <- 2*cov(l[,i],l[,j])/(var(l[,i])+var(l[,j]))
      }
    }
  }
  
  if (type=="logx") v <- cov(l)
  
  return(v)
}

meanLN <- function(mu, Sigma, ind) {
  sd <- diag(Sigma)[-ind]
  sr <- Sigma[ind,-ind]
  si <- Sigma[ind, ind]
  mu.sub <- mu[-ind]
  e <- exp(mu.sub - mu[ind] + 0.5*(sd+si-2*sr))
  1 + exp(-mu[ind]+si/2) + sum(e)
}

meanLOGN <- function(lmu, lsigma) {
  exp(-lmu+lsigma^2/2)
}

logVarMC <- function(mu, Sigma, K=100000) {
  x.norm <- mvtnorm::rmvnorm(K, mu, Sigma)
  ex <- cbind(exp(x.norm), 1)
  x <- ex/rowSums(ex)
  lx <- log(x)
  cov(lx)
}

logVarTaylor <- function(mu, Sigma, transf=c("alr", "clr"), order=c("first","second")) {
  transf <- match.arg(transf)
  order <- match.arg(order)
  
  D <- length(mu)
  ones <- rep(1, D)
  emu <- exp(mu)
  if (transf=="alr") {
    ainv <- emu/(1+sum(emu))
  } else {
    ainv <- emu/sum(emu)
  }
  M <- diag(D)-tcrossprod(ones, ainv)
  t2 <- 0
  if (order=="second") {
    mat <- Sigma%*%(tcrossprod(ainv)-diag(ainv))
    t2 <- sum(diag(mat%*%mat))
    #print(t2)
  }
  M%*%tcrossprod(Sigma, M) + 0.5*t2
}

# Includes reference category
logVarTaylorFull <- function(mu, Sigma, transf=c("alr", "clr"), order=c("first", "second")) {
  transf <- match.arg(transf)
  order <- match.arg(order)
  
  D <- length(mu)
  ones <- rep(1, D+1)
  emu <- exp(mu)
  if (transf=="alr") {
    ainv <- emu/(1+sum(emu))
  } else {
    ainv <- emu/sum(emu)
  }
  #cat("note: function has been changed here!\n")
  M <- rbind(diag(D),0)-tcrossprod(ones, ainv)
  #M <- matrix(1, J+1, J) - tcrossprod(ones, ainv)
  t2 <- 0
  if (order=="second") {
    mat <- Sigma%*%(tcrossprod(ainv)-diag(ainv))
    t2 <- sum(diag(mat%*%mat))
    #print(t2)
  }
  M%*%tcrossprod(Sigma, M) + 0.5*t2
}

# Estimating variance using unscented transformation
# https://users.isy.liu.se/en/rt/fredrik/reports/07SSPut.pdf
logVarUnscented <- function(mu, Sigma, transf=c("alr", "clr"), alpha=1e-3, beta=2, kappa=0) {
  transf <- match.arg(transf)
  
  D <- length(mu)
  lambda <- alpha^2*(D+kappa)-D
  
  ones <- rep(1, 2*D+1)
  ones.short <- rep(1, D)
  
  # svd.Sigma <- svd(Sigma)
  # v <- svd.Sigma$v # This is the transpose of u from Hendeby/Gustafsson paper
  # s <- svd.Sigma$d
  #eig <- eigen(Sigma+D+lambda)
  #eig <- eigen(Sigma+D)
  #sqrt.Sigma <- eig$vectors%*%tcrossprod(diag(sqrt(eig$values)), eig$vectors)
  L <- t(chol(D*Sigma))
  
  x.pm <- tcrossprod(ones, mu)
  # pm.term <- sqrt(D+lambda)*tcrossprod(s, ones.short)*v
  # pm.term <- tcrossprod(s, ones.short)*v
  # x.pm[1:D,] <- x.pm[1:D,] + pm.term
  # x.pm[(D+2):(2*D+1),] <- x.pm[(D+2):(2*D+1),] - pm.term
  x.pm[1:D,] <- x.pm[1:D,] + t(L) #sqrt.Sigma
  x.pm[(D+2):(2*D+1),] <- x.pm[(D+2):(2*D+1),] - L #sqrt.Sigma
  
  #omega.pm <- rep(1/(2*(D+lambda)), 2*D+1)
  #omega.pm[D+1] <- lambda/(D+lambda)
  omega.pm <- rep(1, D+1)
  
  z <- t(apply(x.pm, 1, g))
  mu.z <- colSums(z*omega.pm)
  zmm <- z - mu.z
  
  P <- matrix(0, D+1, D+1)
  for (i in 1:(2*D+1)) {
    #a <- ifelse(i==D+1, (1-alpha^2+beta), 0)
    #P <- P + (omega.pm[i]+a)*tcrossprod(zmm[i,])
  }
  
  P <- cov(z)

  return(P)
}

logVarUnscentedW <- function(mu, Sigma, transf=c("alr", "clr")) {
  transf <- match.arg(transf)
  
  D <- length(mu)
    
  ones <- rep(1, 2*D+1)
  ones.short <- rep(1, D)
  
  # svd.Sigma <- svd(Sigma)
  # v <- svd.Sigma$v # This is the transpose of u from Hendeby/Gustafsson paper
  # s <- svd.Sigma$d
  #eig <- eigen(Sigma+D+lambda)
  #eig <- eigen(Sigma+D)
  #sqrt.Sigma <- eig$vectors%*%tcrossprod(diag(sqrt(eig$values)), eig$vectors)
  w0 <- 0.05
  L <- t(chol(D/(1-w0)*Sigma))
  
  x.pm <- tcrossprod(ones, mu)
  # pm.term <- sqrt(D+lambda)*tcrossprod(s, ones.short)*v
  # pm.term <- tcrossprod(s, ones.short)*v
  # x.pm[1:D,] <- x.pm[1:D,] + pm.term
  # x.pm[(D+2):(2*D+1),] <- x.pm[(D+2):(2*D+1),] - pm.term
  x.pm[1:D,] <- x.pm[1:D,] + t(L) #sqrt.Sigma
  x.pm[(D+2):(2*D+1),] <- x.pm[(D+2):(2*D+1),] - L #sqrt.Sigma
  
  omega.pm <- rep((1-w0)/(2*(D)), 2*D+1)
  omega.pm[D+1] <- w0
    
  z <- t(apply(x.pm, 1, g))
  mu.z <- colSums(z*omega.pm)
  zmm <- z - mu.z
  
  P <- matrix(0, D+1, D+1)
  for (i in 1:(2*D+1)) {
    P <- P + omega.pm[i]*tcrossprod(zmm[i,])
  }
  
  return(P)
}

# Transformation of interest: g(x)=log(alr.inv(x))
g <- function(x) {
  ls <- log(1+sum(exp(x)))
  p1 <- x - ls
  c(p1, -ls)
}
