# Saving MLE fit from single-cell dataset to use in simulations

source("~/research/variation/R/mle_multinom_logit/mle_multinom_logit_normal.R")

load("/project/def-kdm147/kdm147/single_cell/large_cleaned_data_G1.RData")

# First argument for this script is the number of features to subset
# Second argument contains output directory
args <- commandArgs(trailingOnly=TRUE)
n.ss <- as.numeric(args[1])
lmr <- as.numeric(args[2])
outdir <- args[3]

dat.ss <- dat[,1:n.ss]

n.cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
mle <- mlePath(dat.ss, tol=1e-4, tol.nr=1e-4, n.lambda = 8, lambda.min.ratio = lmr, n.cores=n.cores)

outfile <- paste0(outdir, "/", "mle_results_p", n.ss, ".RData")
save(mle, dat.ss, file=outfile)
