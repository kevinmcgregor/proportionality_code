# New version of sim_marginal.R where I use proportionality metrics on LR-transformed
# data as suggested by reviewer.
# Does the marginal variation matrix element estimate match up with approximation?

library(ggplot2)
library(dplyr)
source("proportionality_code/mle_multinom_logit/mle_multinom_logit_normal.R")
source("proportionality_code/mle_multinom_logit/variation_functions.R")

# Read in existing simulated data file?
read.existing <- TRUE

if (!read.existing) {
  #set.seed(8963140)
  set.seed(89631) #set.seed(896314)
  # Testing
  library(mvtnorm)
  n <- 1000
  lmu <- c(4,3,3,2,2,1) # c(5,4,3,2,1,0.5)
  lsigma <- diag(c(4,3,3,2,2,1)) # diag(c(4,3,3,2,1,1))
  lsigma[lower.tri(lsigma)] <- lsigma[upper.tri(lsigma)] <- 0.8
  #lsigma <- rWishart(1, 7, diag(lmu)/2)[,,1]
  k <- length(lmu)
  
  # Read depths simulated from log normal
  # Mean/var of read depths on log-scale
  ln.mean <- 9
  ln.sd <- 0.6
  
  # True variation matrix
  v.true <- logitNormalVariation(lmu, lsigma)
  phi.true <-  logitNormalVariation(lmu, lsigma, "phi")
  rho.true <- logitNormalVariation(lmu, lsigma, type = "rho")
  phi.clr.true <-  logitNormalVariation(lmu, lsigma, "phi", lr="clr")
  rho.clr.true <- logitNormalVariation(lmu, lsigma, type = "rho", lr="clr")
  
  # True CLR variance matrix
  H.inv <- qr.solve(diag(k) + matrix(1, k, k))
  Fm <- cbind(diag(k), -1)
  HiF <- H.inv%*%Fm
  lsigma.clr <- crossprod(HiF, lsigma)%*%HiF
  
  n2 <- 2000
  n.rep <- 1000
  yi.arr <- array(0, dim=c(n2, length(lmu)+1, n.rep))
  p.arr <- array(0, dim=c(n2, length(lmu)+1, n.rep))
  ni.mat <- matrix(0, n.rep, n2)
  for (i in 1:n.rep) {
    if (i%%10==0) cat("i =", i, "\n")
    ni.mat[i,] <- rlnorm(n2, ln.mean, ln.sd) #rlnorm(n2, 9, 0.6)
    p.l <- rmvnorm(n2, lmu, lsigma)
    pl.exp <- cbind(exp(p.l), 1)
    p.arr[,,i] <- pl.exp/rowSums(pl.exp)
    yi.arr[,,i] <- mc2d::rmultinomial(n2, ni.mat[i,], p.arr[,,i])
    if (any(yi.arr[,,i]==0)) {
      yi.arr[,,i] <- as.matrix(zCompositions::cmultRepl(yi.arr[,,i], output = "p-counts", z.warning = 0.9999,
                                                        suppress.print=TRUE))
    }
  }
  save.image("data/bias_simulation/new_sim_marginal_data.RData")
} else {
  load("data/bias_simulation/new_sim_marginal_data.RData")
}


bijk <- function(ni.mat, p.arr, l, m, i) {
  e.pt <- mean(1/ni.mat[,i]*(1/p.arr[i,l,]+1/p.arr[i,m,]))
  v1.pt <- 0.25*var(1/ni.mat[,i]*(1/p.arr[i,l,]-p.arr[i,m,]))
  cv.pt <- cov(log(p.arr[i,l,])-log(p.arr[i,m,]), 1/ni.mat[,i]*(1/p.arr[i,l,]-p.arr[i,l,]))
  return(e.pt+v1.pt-cv.pt)
}

bij <- function(ni.mat, p.arr, l, i) {
  D <- dim(p.arr)[2]
  p <- p.arr[i,l,]
  lg.p <- colMeans(log(p.arr[i,,]))
  p.inv <- 1/p.arr[i,,]
  rm.p.inv <- colMeans(p.inv)
  ni <- ni.mat[,i]
  v <- var(0.5/ni*(p.inv[l,]-rm.p.inv))
  cv <- cov(log(p)-lg.p, 0.5/ni*(p.inv[l,]-rm.p.inv))
  ev <- (((D-1)/D)^2*mean(p.inv[l,])+1/(D^2)*sum(rowMeans(p.inv[-l,])))*mean(1/ni)
  return(v-2*cv+ev)
}

b <- var.lr.dist <- phi.dist <- rho.dist <-
  phi.clr.dist <- rho.clr.dist <- array(0, dim=c(k+1, k+1, n.rep))
for (l in 2:(k+1)) {
  cat("l =", l, "\n")
  for (m in 1:(l-1)) {
    #cat("m =", m, "\n")
    for (i in 1:n.rep) {
      b[l,m,i] <- bijk(ni.mat, p.arr, l, m, i)
      
      # Empirical ALR
      var.lr.dist[l,m,i] <- var(log(yi.arr[,l,i]) - log(yi.arr[,m,i]))
      phi.dist[l,m,i] <- var.lr.dist[l,m,i]/var(log(yi.arr[,l,i]) - log(yi.arr[,k+1,i]))
      nm <- 2*cov(log(yi.arr[,l,i])-log(yi.arr[,k+1,i]),log(yi.arr[,m,i])-log(yi.arr[,k+1,i]))
      dn <- (var(log(yi.arr[,l,i])-log(yi.arr[,k+1,i]))+var(log(yi.arr[,m,i])-log(yi.arr[,k+1,i])))
      rho.dist[l,m,i] <- nm/dn
      
      # Empirical CLR
      log.gmeans <- rowMeans(log(yi.arr[,,i]))
      phi.clr.dist[l,m,i] <- var.lr.dist[l,m,i]/var(log(yi.arr[,l,i]) - log.gmeans)
      nm.clr <- 2*cov(log(yi.arr[,l,i])-log.gmeans,log(yi.arr[,m,i])-log.gmeans)
      dn.clr <- (var(log(yi.arr[,l,i])-log.gmeans)+var(log(yi.arr[,m,i])-log.gmeans))
      rho.clr.dist[l,m,i] <- nm.clr/dn.clr
    }
  }
}

var.clr.dist <- matrix(0, n.rep, k+1)
for (l in 1:(k+1)) {
  cat("l =", l, "\n")
  for (i in 1:n.rep) {
    log.gmeans <- rowMeans(log(yi.arr[,,i]))
    var.clr.dist[i,l] <- var(log(yi.arr[,l,i]) - log.gmeans)
  }
}

# Creating a few matrices to help calculate biases
ones <- rep(1, k)
ones.long <- rep(1, k+1)
d.S <- diag(lsigma)
oS <- rbind(cbind(tcrossprod(d.S, ones), 0), 0)
# Same for CLR
d.Sc <- diag(lsigma.clr)

# Bias corrected values
v.corrected <- var.lr.dist-b
phi.corrected <- rho.corrected <-
  phi.clr.corrected <- rho.clr.corrected <- array(0, dim=c(k+1,k+1,n.rep))
for (i in 1:n.rep) {
  if (i%%10==0) cat("i =", i, "\n")
  br.cur <- b[k+1,1:k,i]
  b.ref <- tcrossprod(b[k+1,,i], ones.long)
  phi.corrected[,,i] <- phi.dist[,,i] - (oS*b[,,i] - v.true*b.ref)/(oS*(oS+b.ref))
  for (l in 2:(k+1)) {
    bij.c <- bij(ni.mat, p.arr, l, i)
    for (m in 1:(l-1)) {
      bik.c <- bij(ni.mat, p.arr, m, i)
      num <- v.true[l,m]*(br.cur[l]+br.cur[m]) - (d.S[l]+d.S[m])*b[l,m,i]
      den <- d.S[l]^2 + d.S[m]^2 + (d.S[l]+d.S[m])*(br.cur[l]+br.cur[m]) +
              2*d.S[l]*d.S[m]
      rho.corrected[l,m,i] <- rho.dist[l,m,i] - num/den
      phi.clr.corrected[l,m,i] <- phi.clr.dist[l,m,i] - 
              (d.Sc[l]*b[l,m,i]-v.true[l,m]*bij.c)/(d.Sc[l]*(d.Sc[l]+bij.c))
      
      num.clr <- v.true[l,m]*(bij.c+bik.c) - b[l,m,i]*(d.Sc[l]+d.Sc[m])
      den.clr <- (d.Sc[l]+d.Sc[m])*(d.Sc[l]+bij.c+d.Sc[m]+bik.c)
      rho.clr.corrected[l,m,i] <- rho.clr.dist[l,m,i] - num.clr/den.clr
    }
  }
}

var.clr.corrected <- matrix(0, n.rep, k+1)
for (i in 1:n.rep) {
  for (l in 1:(k+1)) {
    bij.c <- bij(ni.mat, p.arr, l, i)
    var.clr.corrected[i,l] <- var.clr.dist[i,l] - bij.c
  }
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Bias plot for main manuscript file ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
el <- c(3, 1)

# ALR
vb <- var.lr.dist[el[1],el[2],]-v.true[el[1],el[2]]
vc <- v.corrected[el[1],el[2],]-v.true[el[1],el[2]]
phb <- phi.dist[el[1],el[2],]-phi.true[el[1],el[2]]
phc <- phi.corrected[el[1],el[2],]-phi.true[el[1],el[2]]
rhb <- rho.dist[el[1],el[2],]-rho.true[el[1],el[2]]
rhc <- rho.corrected[el[1],el[2],]-rho.true[el[1],el[2]]

# CLR
# Don't need v for CLR since it's the same as ALR
phb.clr <- phi.dist[el[1],el[2],]-phi.true[el[1],el[2]]
phc.clr <- phi.clr.corrected[el[1],el[2],]-phi.clr.true[el[1],el[2]]
rhb.clr <- rho.clr.dist[el[1],el[2],]-rho.clr.true[el[1],el[2]]
rhc.clr <- rho.clr.corrected[el[1],el[2],]-rho.clr.true[el[1],el[2]]

df.alr <- data.frame(Error=c(vb,phb,rhb,vc,phc,rhc),
                 type=rep(c("Not corrected", "Corrected"), each=n.rep*3),
                 metric=rep(rep(c("v", "phi", "rho"), each=n.rep), 2),
                 LR=rep("ALR", n.rep*6))
df.alr$type <- factor(df.alr$type, levels=c("Not corrected", "Corrected"))
df.alr$metric <- factor(df.alr$metric, levels=c("v", "phi", "rho"),
                    labels=c("v", expression(phi), expression(rho)))

df.clr <- data.frame(Error=c(vb,phb.clr,rhb.clr,vc,phc.clr,rhc.clr),
          type=rep(c("Not corrected", "Corrected"), each=n.rep*3),
          metric=rep(rep(c("v", "phi", "rho"), each=n.rep), 2),
          LR=rep("CLR", n.rep*6))

df <- rbind(df.alr, df.clr)
df$type <- factor(df$type, levels=c("Not corrected", "Corrected"))
df$metric <- factor(df$metric, levels=c("v", "phi", "rho"),
                    labels=c("v", expression(phi), expression(rho)))
df$LR <- factor(df$LR, levels=c("ALR", "CLR"))

p <- ggplot(df, aes(x=type, y=Error, fill=type)) +
  facet_grid(rows=vars(LR), cols=vars(metric), labeller=label_parsed,
             scales="free_y") +
  scale_fill_manual(values=c("#61D04F", "#2297E6")) +
  geom_boxplot() +
  theme_bw() +
  theme(strip.background=element_rect(fill="white"),
        strip.text.x = element_text(angle = 0, size = 16),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title=element_text(size=16,face="bold"),
        legend.key.size = unit(1.5, 'cm'),
        legend.text = element_text(size=12)) +
  geom_hline(yintercept = 0, color="red", linetype="dashed")
ggsave("proportionality_code/bias_sim/plots/main_bias_plot_new.pdf",
       p, height=6, width=10)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Bias plots all elements - supplementary ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n.plt <- 7

pdf("proportionality_code/bias_sim/plots/v_bias_plot_new.pdf",
    height=10, width=12)
par(mar=c(0.5,3,0.5,1))
layout(matrix(c(rep(1,n.plt-1), 2:((n.plt-1)^2+1)), n.plt, n.plt-1, byrow = TRUE),
       heights=c(0.5, rep(1, n.plt-1)))
plot.new()
text(0.5,0.5,"Bias in v",cex=4,font=2)
for (l in 1:n.plt) {
  for (m in 1:n.plt) {
    if (l>m) {
      bp <- cbind(var.lr.dist[l,m,]-v.true[l,m],
                  v.corrected[l,m,]-v.true[l,m])
      boxplot(bp, col=3:4, ylim=c(-0.6,0.6), xaxt="n")
      abline(h=0, lwd=2, lty=2, col=2)
    } else if (l==2 & m==n.plt-1) {
      plot(-100, -100, xaxt="n", yaxt="n", col="white", bty="n")
      legend("top", fill=c(3,4), legend=c("Empirical", "Bias corrected"), 
             cex=2, bty="n")
    } else if (l==1 | m==n.plt) {
      # Do nothing
    } else {
      plot.new()
    }
  }
}
dev.off()

# Phi
# Since ref category for Phi after ALR is Inf, removing
# that row/column
pdf("proportionality_code/bias_sim/plots/phi_bias_plot_new.pdf",
    height=10, width=12)
par(mar=c(0.5,3,0.5,1))
n.plt.tmp <- n.plt-1
layout(matrix(c(rep(1,n.plt.tmp-1), 2:((n.plt.tmp-1)^2+1)), n.plt.tmp, n.plt.tmp-1, byrow = TRUE),
       heights=c(0.5, rep(1, n.plt.tmp-1)))
plot.new()
text(0.5,0.5,expression(paste("Bias in ", phi)),cex=4,font=2)
for (l in 1:n.plt.tmp) {
  for (m in 1:n.plt.tmp) {
    if (l>m) {
      bp <- cbind(phi.dist[l,m,]-phi.true[l,m],
                  phi.corrected[l,m,]-phi.true[l,m])
      boxplot(bp, col=3:4, ylim=c(-1,1), xaxt="n")
      abline(h=0, lwd=2, lty=2, col=2)
    } else if (l==2 & m==n.plt.tmp-1) {
      plot(-100, -100, xaxt="n", yaxt="n", col="white", bty="n")
      legend("top", fill=c(3,4), legend=c("Empirical", "Bias corrected"), 
             cex=2, bty="n")
    } else if (l==1 | m==n.plt.tmp) {
      # Do nothing
    } else {
      plot.new()
    }
  }
}
dev.off()

# Rho
# Since ref category for Rho after ALR is Inf, removing
# that row/column
pdf("proportionality_code/bias_sim/plots/rho_bias_plot_new.pdf",
    height=10, width=12)
par(mar=c(0.5,3,0.5,1))
layout(matrix(c(rep(1,n.plt.tmp-1), 2:((n.plt.tmp-1)^2+1)), n.plt.tmp, n.plt.tmp-1, byrow = TRUE),
       heights=c(0.5, rep(1, n.plt.tmp-1)))
plot.new()
text(0.5,0.5,expression(paste("Bias in ", rho)),cex=4,font=2)
for (l in 1:n.plt.tmp) {
  for (m in 1:n.plt.tmp) {
    if (l>m) {
      bp <- cbind(rho.dist[l,m,]-rho.true[l,m],
                  rho.corrected[l,m,]-rho.true[l,m])
      boxplot(bp, col=3:4, ylim=c(-0.3,0.3), xaxt="n")
      abline(h=0, lwd=2, lty=2, col=2)
    } else if (l==2 & m==n.plt.tmp-1) {
      plot(-100, -100, xaxt="n", yaxt="n", col="white", bty="n")
      legend("top", fill=c(3,4), legend=c("Empirical", "Bias corrected"), 
             cex=2, bty="n")
    } else if (l==1 | m==n.plt.tmp) {
      # Do nothing
    } else {
      plot.new()
    }
  }
}
dev.off()

# Phi - CLR
pdf("proportionality_code/bias_sim/plots/phi_clr_bias_plot_new.pdf",
    height=10, width=12)
par(mar=c(0.5,3,0.5,1))
layout(matrix(c(rep(1,n.plt-1), 2:((n.plt-1)^2+1)), n.plt, n.plt-1, byrow = TRUE),
       heights=c(0.5, rep(1, n.plt-1)))
plot.new()
text(0.5,0.5,expression(paste("Bias in ", phi)),cex=4,font=2)
for (l in 1:n.plt) {
  for (m in 1:n.plt) {
    if (l>m) {
      bp <- cbind(phi.clr.dist[l,m,]-phi.clr.true[l,m],
                  phi.clr.corrected[l,m,]-phi.clr.true[l,m])
      boxplot(bp, col=3:4, ylim=range(bp), xaxt="n")
      abline(h=0, lwd=2, lty=2, col=2)
    } else if (l==2 & m==n.plt-1) {
      plot(-100, -100, xaxt="n", yaxt="n", col="white", bty="n")
      legend("top", fill=c(3,4), legend=c("Empirical", "Bias corrected"), 
             cex=2, bty="n")
    } else if (l==1 | m==n.plt) {
      # Do nothing
    } else {
      plot.new()
    }
  }
}
dev.off()

# Rho - CLR
pdf("proportionality_code/bias_sim/plots/rho_clr_bias_plot_new.pdf",
    height=10, width=12)
par(mar=c(0.5,3,0.5,1))
layout(matrix(c(rep(1,n.plt-1), 2:((n.plt-1)^2+1)), n.plt, n.plt-1, byrow = TRUE),
       heights=c(0.5, rep(1, n.plt-1)))
plot.new()
text(0.5,0.5,expression(paste("Bias in ", rho)),cex=4,font=2)
for (l in 1:n.plt) {
  for (m in 1:n.plt) {
    if (l>m) {
      bp <- cbind(rho.clr.dist[l,m,]-rho.clr.true[l,m],
                  rho.clr.corrected[l,m,]-rho.clr.true[l,m])
      boxplot(bp, col=3:4, ylim=c(-0.1,0.1), xaxt="n")
      abline(h=0, lwd=2, lty=2, col=2)
    } else if (l==2 & m==n.plt-1) {
      plot(-100, -100, xaxt="n", yaxt="n", col="white", bty="n")
      legend("top", fill=c(3,4), legend=c("Empirical", "Bias corrected"), 
             cex=2, bty="n")
    } else if (l==1 | m==n.plt) {
      # Do nothing
    } else {
      plot.new()
    }
  }
}
dev.off()


# Bias plot for CLR variance
pdf("proportionality_code/bias_sim/plots/variance_clr_bias_plot_new.pdf",
height=10, width=12)
par(mar=c(0.5,3,0.5,1))
layout(matrix(c(rep(1,n.plt-1), 2, 3:(n.plt+3-1)), 2, n.plt, byrow = TRUE),
       heights=c(0.5, rep(1, n.plt-1)))
plot.new()
text(0.5,0.5,expression(paste("Bias in CLR variance")),cex=4,font=2)
plot(-100, -100, xaxt="n", yaxt="n", col="white", bty="n")
legend("top", fill=c(3,4), legend=c("Empirical", "Bias corrected"), 
       cex=2, bty="n")
for (i in 1:n.plt) {
  boxplot(var.clr.dist[,i]-d.Sc[i], var.clr.corrected[,i]-d.Sc[i],
          col=3:4, xaxt="n", ylim=c(-1,1))
  abline(h=0, lwd=2, lty=2, col=2)
}
dev.off()

