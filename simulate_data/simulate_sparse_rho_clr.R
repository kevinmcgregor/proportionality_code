# Simulate from Multinom LN model using parameters from single-cell dataset
# Making zeros on CLR scale

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
Sigma.orig <- mle$est.min$Sigma

# Calculating CLR covariance
Fm <- cbind(diag(n.g-1), -1)
H <- diag(n.g-1) + 1
Hinv <- qr.solve(H)
clr.Sigma <- crossprod(Hinv%*%Fm, Sigma.orig)%*%Hinv%*%Fm

# Mapping between zero probs in Cholesky decomp and
# zero probs in Sigma
if (n.g==50) {
   if (pr.0==0.15) pr.0.chol <- 0.61
   if (pr.0==0.5) pr.0.chol <- 0.8
   if (pr.0==0.85) pr.0.chol <- 0.93
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
L0 = chol(clr.Sigma+0.1*diag(n.g))
sp0 = rbinom(choose(n.g, 2), 1, prob=pr.0.chol)
L0[upper.tri(L0)] = (1-sp0)*L0[upper.tri(L0)]

# Making diagonals smaller in CLR so resulting rho values will be
# further from zero
diag(L0) <- diag(L0)/5

#Getting final values of sparse CLR Sigma (meaning rho will be sparse as well)
clr.Sigma = crossprod(L0)
d.clr.S <- diag(clr.Sigma)

# Getting back to original Sigma scale
Sigma <- Fm%*%tcrossprod(clr.Sigma, Fm)

# "True" values of phi, phis, and rho.
v.true <- logitNormalVariation(mu, Sigma, type="standard", lr="clr")
phi.true <- v.true/tcrossprod(d.clr.S, rep(1, n.g))
rho.true <- 2*clr.Sigma/outer(d.clr.S, d.clr.S, "+")
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
if (any(y==0)) {
   y.no0 <- as.matrix(zCompositions::cmultRepl(y, output = "p-count"))
} else {
  y.no0 <- y
}

# Create output folder if it doesn't already exist
odir <- paste0(outdir, "/p0_", pr.0)
cmd <- paste0("mkdir -p ", odir)
system(cmd)

fout <- paste0(odir, "/", "sim_data_clr", sim_id, ".RData")
save.image(fout)