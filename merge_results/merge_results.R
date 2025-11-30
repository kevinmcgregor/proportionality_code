# Calculating performance and merging simulation results
library(Bolstad2)
library(PRROC)

source("/home/kdm147/research/variation/R/prec_recall_functions.R")

# Base input directory
bdir <- "/scratch/kdm147/sim_sparse_rho/"

# n/p subdirectories
np.dir <- c("n100_p50", "n200_p50", "n500_p50", "n100_p100", "n200_p100",
       	  "n500_p100", "n100_p200", "n200_p200", "n500_p200")

#np.dir <- c("n500_p50", "n500_p100", "n500_p200")

#np.dir <- "n100_p50"

# Sparsity probability subdirectories
#sp.dir <- c("p0_0.5")
sp.dir <- c("p0_0.15", "p0_0.5", "p0_0.85")
n.spar <- length(sp.dir)

# gamma (from EBIC) subdirectories
gamma.dir <- "gamma_0.1"

# Subdirectories to loop over
subdir <- expand.grid(np.dir, sp.dir, gamma.dir, stringsAsFactors = FALSE)
subdir <- apply(subdir, 1, function(x){paste(x, collapse="/")})

n.run <- length(subdir)

# Base of input filename
bname <- "var"

# Number of simulation reps
n.rep <- 50
# Number of thresholds for prec/recall
n.thresh <- 500
# p-value cutoff for Pearson correlation test
p.cut <- 0.05

aupr.naive <- aupr.est <- matrix(NA, n.rep, n.run)
maxf.naive <- maxf.est <- matrix(NA, n.rep, n.run)
prc.naive <- prc.est <- array(NA, dim=c(n.rep, n.run, n.thresh, 2))
mse.naive <- mse.est <- matrix(0, n.rep, n.run)

pr.naive <- pr.est <- array(NA, dim=c(n.rep, n.run, 2))
f.naive <- f.est <- matrix(NA, n.rep, n.run)

for (r in 1:n.run) {
    cat("r =", r, "/", n.run, "\n")
    for (s in 1:n.rep) {
    if (s%%10==0) cat("\t s =", s, "/", n.rep, "\n")
	       
    fn <- paste0(bdir, "/", subdir[r], "/", bname, s, ".RData")

	  tryCatch({
      fn <- paste0(bdir, "/", subdir[r], "/", bname, s, ".RData")
      
      load(fn)
    
      prc.naive[s,r,,] <- getPRC(rho.naive, rho.true, n.thresh)
      prc.est[s,r,,] <- getPRC(rho.est[-n.g, -n.g], rho.true, n.thresh)
    
    	maxf.naive[s,r] <- max_f_score(prc.naive[s,r,,])
    	maxf.est[s,r]	<- max_f_score(prc.est[s,r,,])
    
    	# Prec/recall
      lt <- lower.tri(rho.true)
      true.class <- ifelse(rho.true[lt]==0, 0, 1)
    
    	# TMP - don't need the -n.g if calculating rho.est using Sigma.hat directly 
    	pred.est <- ifelse(rho.est[-n.g, -n.g][lt]==0, 0, 1)
    	#pred.est <- ifelse(rho.est[lt]==0, 0, 1)
    
    	ct.naive <- cor_test_matrix(y.alr)
    	pred.naive <- ifelse(ct.naive[lt]>=p.cut, 0, 1)
    
    	pr.naive[s,r,] <- getPR(pred.naive, true.class)
    	pr.est[s,r,] <- getPR(pred.est, true.class)
    	f.naive[s,r] <- f_score(pr.naive[s,r,])
    	f.est[s,r] <- f_score(pr.est[s,r,])
    
      aupr.naive[s,r] <- suppressWarnings(sintegral(prc.naive[s,r,,2], prc.naive[s,r,,1])$int)
      aupr.est[s,r] <- suppressWarnings(sintegral(prc.est[s,r,,2], prc.est[s,r,,1])$int)
            
      rt <- rho.true
      re <- rho.est[-n.g,-n.g]
      mse.naive[s,r] <- mean((rho.naive[lt]-rt[lt])^2)
      mse.est[s,r] <- mean((re[lt]-rt[lt])^2)}, error=function(x){}, warning=function(x){})
	
     }
}

save.dir <- "/home/kdm147/research/variation/sim_sparse_rho/results"
ofile <- paste0(save.dir, "/merged_sim_sparse_rho.RData")
print(ofile)
save.image(ofile)

rnames <- rep(np.dir, length(sp.dir))

f.out <- cbind(colMeans(f.est, na.rm=TRUE), colMeans(f.naive, na.rm=TRUE))
colnames(f.out) <- c("F-score Est", "F-score Naive")
rownames(f.out) <- rnames
print(f.out)
cat("\n\n")

e.out <- cbind(colMeans(mse.est, na.rm=TRUE), colMeans(mse.naive, na.rm=TRUE))
colnames(e.out) <- c("MSE Est", "MSE Naive")
rownames(e.out) <- rnames
print(e.out)
