# Simulate from Multinom LN model using parameters from single-cell dataset

source("/home/kdm147/research/variation/R/variation_functions.R")
source("/home/kdm147/research/variation/R/mle_multinom_logit/mle_multinom_logit_normal.R")

# Read in command line arguments
args <- commandArgs(trailingOnly = TRUE)
sim_id <- as.numeric(args[1])
ni_ratio <- as.numeric(args[2])
n <- as.numeric(args[3])
pr.0 <- as.numeric(args[4])
infile <- args[5]
outdir <- args[6]

cat("sim_id =", sim_id, "\n")
cat("ni_ratio =", ni_ratio, "\n")
cat("n =", n, "\n")
cat("pr.0 =", pr.0, "\n")
cat("infile =", infile, "\n")
cat("outdir =", outdir, "\n")

# Load in MLE data
load(infile)

# Number of genes in original dataset
n.g <- length(mle$est.min$mu)+1

mu <- mle$est.min$mu
Sigma <- mle$est.min$Sigma

# Mapping between zero probs in Cholesky decomp and
# zero probs in Sigma
if (n.g==50) {
   if (pr.0==0.15) pr.0.chol <- 0.61
   if (pr.0==0.5) pr.0.chol <- 0.8
   if (pr.0==0.85) pr.0.chol <-	0.93
} else if (n.g==100) {
   if (pr.0==0.15) pr.0.chol <- 0.7
   if (pr.0==0.5) pr.0.chol <- 0.86
   if (pr.0==0.85) pr.0.chol <- 0.945
} else if (n.g==200) {
   if (pr.0==0.15) pr.0.chol <- 0.79
   if (pr.0==0.5) pr.0.chol <- 0.9
   if (pr.0==0.85) pr.0.chol <- 0.96
} else {
   stop("Invalid entry for pr.0")
}


# Getting Cholesky decomp and setting values to zero
L0 = chol(Sigma)
sp0 = rbinom(choose(n.g-1, 2), 1, prob=pr.0.chol)
L0[upper.tri(L0)] = (1-sp0)*L0[upper.tri(L0)]

#Getting final values of sparse Sigma (meaning rho will be sparse as well)
Sigma = crossprod(L0)
d.S <- diag(Sigma)

# "True" values of phi, phis, and rho.
ones <- rep(1, n.g-1)
# V
v.true <- tcrossprod(d.S, ones) + tcrossprod(ones, d.S) - 2*Sigma
# Phi
lv.row <- tcrossprod(d.S, ones)
phi.true <- v.true/lv.row
# Rho
den.rho <- outer(d.S, d.S, "+")
rho.true <- 2*Sigma/den.rho
# Phis
phis.true <- (1-rho.true)/(1+rho.true)

# Simulating data
x.logit <- mvtnorm::rmvnorm(n, mu, Sigma)
xl.exp <- cbind(exp(x.logit), 1)
x <- xl.exp/rowSums(xl.exp)

rd <- rowSums(dat.ss)/ni_ratio
ln.mean <- mean(log(rd))
ln.var <- var(log(rd))
ni <- rlnorm(n, ln.mean, ln.var)

y <- mc2d::rmultinomial(n, ni, x)

# Filtering out genes with too many zeros
thresh <- 0.2
pzero <- apply(y, 2, function(x){mean(x==0)})
incl <- pzero<thresh
y <- y[,incl]

cat("dim(phi.true)", dim(phi.true), "\n")
cat("sum(incl) =", sum(incl), "\n")
cat("length(incl) =", length(incl), "\n")

# Subsetting true values
v.true <- v.true[incl[-n.g], incl[-n.g]]
phi.true <- phi.true[incl[-n.g], incl[-n.g]]
phis.true <- phis.true[incl[-n.g], incl[-n.g]]
rho.true <- rho.true[incl[-n.g], incl[-n.g]]

# Number of genes kept after filtering
n.g <- NCOL(y)

# Bayesian imputation of zeros
y.no0 <- as.matrix(zCompositions::cmultRepl(y, output = "p-count"))

# Create output folder if it doesn't already exist
odir <- paste0(outdir, "/p0_", pr.0)
cmd <- paste0("mkdir -p ", odir)
system(cmd)

fout <- paste0(odir, "/", "sim_data", sim_id, ".RData")
save.image(fout)