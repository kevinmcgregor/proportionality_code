# Run the various variation estimation methods

library(compositions)

# Read in command line arguments
args <- commandArgs(trailingOnly = TRUE)
sim_id <- as.numeric(args[1])
lmr <- as.numeric(args[2])
gamma <- as.numeric(args[3])
pr.0 <- as.numeric(args[4])
indir <- args[5]

# Load simulated data on which to run the methods
idr <- paste0(indir, "/", "p0_", pr.0)
input_file <- paste0(idr, "/sim_data", sim_id, ".RData")
load(input_file)

# Scripts
source("/home/kdm147/research/variation/R/variation_functions.R")
source("/home/kdm147/research/variation/R/mle_multinom_logit/mle_multinom_logit_normal.R")

# Did we simulate from t-dist?
tsim <- exists("nc.adj")

# Get sample size
n.samp <- NROW(y)

n.lam <- 8
n.cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
if (tsim & n.g>50 & n.g <= 101) {
   # This is the p100, tdist case
   	lam <- seq(0.7, 1.6, length.out=n.lam)
	mle.sim <- mlePath(y, tol=1e-4, tol.nr=1e-4, lambda.gl=lam, gamma = gamma, n.cores=n.cores)
} else if (n.samp==500 & pr.0==0.85 & n.g<=50) {
  #This is the p50 case
  lam <- seq(0.1, 0.28, length.out=n.lam)
  mle.sim <- mlePath(y, tol=1e-4, tol.nr=1e-4, lambda.gl=lam, gamma = gamma, n.cores=n.cores)
} else if (n.samp==500 & pr.0==0.85 & n.g<=100) {
  # p100 Gaussian case
  lam <- seq(0.8, 1.5, length.out=n.lam)
  mle.sim <- mlePath(y, tol=1e-4, tol.nr=1e-4, lambda.gl=lam, gamma = gamma, n.cores=n.cores)
} else if (n.samp==500 & pr.0==0.85) {
  lam <- seq(0.28, 1, length.out=n.lam)
  mle.sim <- mlePath(y, tol=1e-4, tol.nr=1e-4, lambda.gl=lam, gamma = gamma, n.cores=n.cores)
} else if (n.samp==100 & n.g<=50 & pr.0==0.5) {
  lam <- seq(0.06, 0.4, length.out=n.lam)
  mle.sim <- mlePath(y, tol=1e-4, tol.nr=1e-4, lambda.gl=lam, gamma = gamma, n.cores=n.cores)
} else {
  mle.sim <- mlePath(y, tol=1e-4, tol.nr=1e-4, n.lambda = n.lam, lambda.min.ratio = lmr, gamma = gamma,
 n.cores=n.cores)
}

if (!tsim) {
   if (pr.0==0.85) {
      mu.hat <- mle.sim$est.min$mu
      Sigma.hat <- mle.sim$est.min$Sigma
   } else {
     mu.hat <- mle.sim$est[[2]]$mu
      Sigma.hat <- mle.sim$est[[2]]$Sigma
   }
} else {
  if (pr.0==0.85) {
      mu.hat <- mle.sim$est[[4]]$mu
      Sigma.hat <- mle.sim$est[[4]]$Sigma
   } else if (n.g>50 & n.g <= 101) {
     # This is the p100 case
      mu.hat <- mle.sim$est.min$mu
      Sigma.hat <- mle.sim$est.min$Sigma
   } else {
      mu.hat <- mle.sim$est.min$mu
      Sigma.hat <- mle.sim$est.min$Sigma
   }
}


# Naive versions
v.naive <- naiveVariation(y.no0)
phi.naive <- naiveVariation(y.no0, type="phi")
phis.naive <- naiveVariation(y.no0, type="phis")
rho.naive <- naiveVariation(y.no0, type="rho")

# Estimated versions
v.est <- logitNormalVariation(mu.hat, Sigma.hat)
phi.est <- logitNormalVariation(mu.hat, Sigma.hat, type="phi")
phis.est <- logitNormalVariation(mu.hat, Sigma.hat, type="phis")
rho.est <- logitNormalVariation(mu.hat, Sigma.hat, type="rho")

# Estimating mu/Sigma of ALR-transformed data, then plugging into
# log.p approximation
y.alr <- log(y.no0[,-n.g]) - log(y.no0[,n.g])
mu.emp <- colMeans(y.alr)
Sigma.emp <- cov(y.alr)
rho.plugin <- logitNormalVariation(mu.emp, Sigma.emp, type="rho")

# Create output folder if it doesn't already exist
od <- paste0(idr, "/gamma_", gamma)
cmd <- paste0("mkdir -p ", od)
system(cmd)

savefile <- paste0(od, "/var", sim_id, ".RData")
save.image(savefile)