# Calculating performance and merging simulation results
library(Bolstad2)
library(PRROC)

source("/home/kdm147/research/variation/R/prec_recall_functions.R")

# Base input directory
bdir <- "/scratch/kdm147/sim_sparse_rho/"

# n/p subdirectories
np.dir <- c("n100_p50", "n200_p50", "n500_p50", "n100_p100", "n200_p100",
       	  "n500_p100", "n100_p200", "n200_p200", "n500_p200")

#np.dir <- c("n500_p200")


# Sparsity probability subdirectories
sp.dir <- c("p0_0.15", "p0_0.5", "p0_0.85")
#sp.dir <- c("p0_0.85")
n.spar <- length(sp.dir)

# gamma (from EBIC) subdirectories
gamma.dir <- "gamma_0.1"

# Subdirectories to loop over
subdir <- expand.grid(np.dir, sp.dir, gamma.dir, stringsAsFactors = FALSE)
subdir <- apply(subdir, 1, function(x){paste(x, collapse="/")})

n.run <- length(subdir)

# Base of input filename
bname <- "var_clr"

# Number of simulation reps
n.rep <- 50
# Number of thresholds for prec/recall
n.thresh <- 100
# p-value cutoff for Pearson correlation test
p.cut <- 0.05

aupr.naive <- aupr.est <- aupr.plugin <- matrix(NA, n.rep, n.run)
prc.naive <- prc.est <- prc.plugin <- array(NA, dim=c(n.rep, n.run, n.thresh, 2))

mse.naive <- mse.est <- matrix(NA, n.rep, n.run)

pr.naive <- pr.est <- array(NA, dim=c(n.rep, n.run, 2))
f.naive <- f.est <- matrix(NA, n.rep, n.run)
f2.naive <- f2.est <- matrix(NA, n.rep, n.run)
fh.naive <- fh.est <- matrix(NA, n.rep, n.run)

for (r in 1:n.run) {
    cat("r =", r, "/", n.run, "\n")
    for (s in 1:n.rep) {
    	if (s%%10==0) cat("\t s =", s, "/", n.rep, "\n")
	       
    	fn <- paste0(bdir, "/", subdir[r], "/", bname, s, ".RData")

	tryCatch({
		load(fn)

		lt <- lower.tri(rho.true)
                b.true.lt <- ifelse(rho.true[lt]!=0, 1, 0)

                #aupr.naive[s,r] <- pr.curve(abs(rho.naive[lt]), weights.class0=b.true.lt)$auc.integral
                #aupr.est[s,r] <- pr.curve(abs(rho.est[-n.g, -n.g][lt]), weights.class0=b.true.lt)$auc.integral
                #aupr.plugin[s,r] <- pr.curve(abs(rho.plugin[-n.g, -n.g][lt]), weights.class0=b.true.lt)$auc.integral


		# TMP - trying another lambda value for rho est
		#idx <- ifelse(mle.sim$min.idx>=3, mle.sim$min.idx, mle.sim$min.idx+6)
		#Sigma.clr.hat <- mle.sim$est[[idx]]$Sigma.clr
		#d.Sc.hat <- diag(Sigma.clr.hat)
		#rho.est <- 2*Sigma.clr.hat/outer(d.Sc.hat, d.Sc.hat, "+")

		true.class <- ifelse(rho.true[lt]==0, 0, 1)
		pred.est <- ifelse(rho.est[lt]==0, 0, 1)

		y.clr <- compositions::clr(y.no0)
		ct.naive <- cor_test_matrix(y.clr)
		pred.naive <- ifelse(ct.naive[lt]>=p.cut, 0, 1)

        	pr.naive[s,r,] <- getPR(pred.naive, true.class)
        	pr.est[s,r,] <- getPR(pred.est, true.class)
		f.naive[s,r] <- f_score(pr.naive[s,r,])
                f.est[s,r] <- f_score(pr.est[s,r,])
        	f2.naive[s,r] <- f_score(pr.naive[s,r,], beta=2)
        	f2.est[s,r] <- f_score(pr.est[s,r,], beta=2)
		fh.naive[s,r] <- f_score(pr.naive[s,r,], beta=0.5)
                fh.est[s,r] <- f_score(pr.est[s,r,], beta=0.5)

		lt2 <- lower.tri(rho.naive)
		rt <- rho.true[-n.g,-n.g]
		re <- rho.est[-n.g,-n.g]
		mse.naive[s,r] <- mean((rho.naive[lt2]-rt[lt2])^2)
		mse.est[s,r] <- mean((re[lt2]-rt[lt2])^2)

		#prc.naive[s,r,,] <- getPRC(rho.naive, rho.true[-n.g,-n.g], n.thresh)
        	#prc.est[s,r,,] <- getPRC(rho.est[-n.g,-n.g], rho.true[-n.g,-n.g], n.thresh)

		
		#aupr.naive[s,r] <- suppressWarnings(sintegral(prc.naive[s,r,,2], prc.naive[s,r,,1])$int)
		#aupr.est[s,r] <- suppressWarnings(sintegral(prc.est[s,r,,2], prc.est[s,r,,1])$int)
	}, error=function(x){}, warning=function(x){})

   }
}

save.dir <- "/home/kdm147/research/variation/sim_sparse_rho/results"
ofile <- paste0(save.dir, "/merged_sim_sparse_rho_clr.RData")
print(ofile)
save.image(ofile)

rnames <- rep(np.dir, length(sp.dir))

#m.out <- cbind(colMeans(aupr.est, na.rm=TRUE), colMeans(aupr.naive, na.rm=TRUE))
#colnames(m.out) <- c("AUPR Est", "AUPR Naive")
#rownames(m.out) <- rnames
#print(m.out)
#cat("\n\n")

f.out <- cbind(colMeans(f.est, na.rm=TRUE), colMeans(f.naive, na.rm=TRUE))
colnames(f.out) <- c("F-score Est", "F-score Naive")
rownames(f.out) <- rnames
print(f.out)
cat("\n\n")

f2.out <- cbind(colMeans(f2.est, na.rm=TRUE), colMeans(f2.naive, na.rm=TRUE))
colnames(f2.out) <- c("F2-score Est", "F2-score Naive")
rownames(f2.out) <- rnames
print(f2.out)
cat("\n\n")

fh.out <- cbind(colMeans(fh.est, na.rm=TRUE), colMeans(fh.naive, na.rm=TRUE))
colnames(fh.out) <- c("F0.5-score Est", "F0.5-score Naive")
rownames(fh.out) <- rnames
print(fh.out)
cat("\n\n")

e.out <- cbind(colMeans(mse.est, na.rm=TRUE), colMeans(mse.naive, na.rm=TRUE))
colnames(e.out) <- c("MSE Est", "MSE Naive")
rownames(e.out)	<- rnames
print(e.out)