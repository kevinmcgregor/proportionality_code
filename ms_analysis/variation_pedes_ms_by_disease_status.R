# Variation estimation in MS pedes data, split up by disease status
source("proportionality_code/mle_multinom_logit/variation_functions.R")
source("proportionality_code/mle_multinom_logit/mle_multinom_logit_normal.R")

load("data/ms/ms_pedes_cleaned.RData")

n.g <- NCOL(ms.genus)

lambda.vec <- exp(seq(-2, 1, length.out=8))

# Running model-based variation on MS samples
ms.genus.ms <- ms.genus[covars$Diagnosis=="MS",]
mle.ms <- mlePath(ms.genus.ms, tol = 1e-4, tol.nr = 1e-4, lambda.gl=lambda.vec, n.cores=8)
save(mle.ms, file="data/ms/ms_pedes_mle_ms_ebic.RData")

# Running model-based variation on control samples
ms.genus.control <- ms.genus[covars$Diagnosis=="Controls",]
mle.control <- mlePath(ms.genus.control, tol = 1e-4, tol.nr = 1e-4, lambda.gl=lambda.vec, n.cores=8)
save(mle.control, file="data/ms/ms_pedes_mle_control_ebic.RData")