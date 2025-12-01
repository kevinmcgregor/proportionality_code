# Run the various variation estimation methods
# CLR version

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
input_file <- paste0(idr, "/sim_data_clr", sim_id, ".RData")
load(input_file)

# Scripts
source("/home/kdm147/research/variation/R/variation_functions.R")
source("/home/kdm147/research/variation/R/mle_multinom_logit/mle_multinom_logit_normal.R")

# Run MLE method
n.lam <- 8
n.cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
#mle.sim <- mlePath(y, tol=1e-4, tol.nr=1e-4, n.lambda = n.lam, lambda.min.ratio = lmr, gamma = gamma, n.cores=8, lr="clr")

if (pr.0==0.85 & n==500 & n.g<=50) {
   mle.sim <- mlePath(y, tol=1e-5, tol.nr=1e-5, lambda.gl=seq(0.0005, 0.0025, length.out=8), gamma = gamma, n.cores=8, lr="clr")
} else if (pr.0==0.85 & n==500 & n.g>100) {
  mle.sim <- mlePath(y, tol=1e-6, tol.nr=1e-6, lambda.gl=seq(0.01, 0.03, length.out=8), gamma = gamma, n.cores=8, lr="clr")
} else if (pr.0==0.85 & (n!=500 | n>50)) {
  mle.sim <- mlePath(y, tol=1e-5, tol.nr=1e-5, lambda.gl=seq(0.01, 0.04, length.out=8), gamma = gamma, n.cores=8, lr="clr")
} else {
  mle.sim <- mlePath(y, tol=1e-4, tol.nr=1e-4, lambda.gl=seq(0.0001, 0.2, length.out=8), gamma = gamma, n.cores=8, lr="clr"
)
}

# DON'T RUN
if (FALSE) {
if (pr.0==0.85) {
  mu.hat <- mle.sim$est[[5]]$mu
  Sigma.hat <- mle.sim$est[[5]]$Sigma
  Sigma.clr.hat <- mle.sim$est[[5]]$Sigma.clr
} else {
  mu.hat <- mle.sim$est.min$mu
  Sigma.hat <- mle.sim$est.min$Sigma
  Sigma.clr.hat <- mle.sim$est.min$Sigma.clr
}
}

mu.hat <- mle.sim$est.min$mu
Sigma.hat <- mle.sim$est.min$Sigma
Sigma.clr.hat <- mle.sim$est.min$Sigma.clr

# Naive versions
v.naive <- naiveVariation(y.no0, lr="clr")[-n.g,-n.g]
phi.naive <- naiveVariation(y.no0, type="phi", lr="clr")[-n.g,-n.g]
phis.naive <- naiveVariation(y.no0, type="phis", lr="clr")[-n.g,-n.g]
rho.naive <- naiveVariation(y.no0, type="rho", lr="clr")[-n.g,-n.g]

# Estimated versions
v.est <- logitNormalVariation(mu.hat, Sigma.hat, lr="clr")
phi.est <- logitNormalVariation(mu.hat, Sigma.hat, type="phi", lr="clr")
phis.est <- logitNormalVariation(mu.hat, Sigma.hat, type="phis", lr="clr")
# rho.est <- logitNormalVariation(mu.hat, Sigma.hat, type="rho", lr="clr")
d.Sc.hat <- diag(Sigma.clr.hat)
rho.est <- 2*Sigma.clr.hat/outer(d.Sc.hat, d.Sc.hat, "+")

# Create output folder if it doesn't already exist
od <- paste0(idr, "/gamma_", gamma)
cmd <- paste0("mkdir -p ", od)
system(cmd)

savefile <- paste0(od, "/var_clr", sim_id, ".RData")
save.image(savefile)