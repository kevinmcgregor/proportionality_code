# Variation estimation in MS pedes data, split up by disease status
# Zero imputing first to check change from original results
source("~/research/variation/R/variation_functions.R")
source("~/research/variation/R/mle_multinom_logit/mle_multinom_logit_normal.R")

load("/project/def-kdm147/kdm147/heather_armstrong/data_2024/ms_pedes_cleaned.RData")

ms.no0 <- as.matrix(zCompositions::cmultRepl(ms.genus, output = "p-counts", z.warning = 0.9999,
                                             suppress.print=TRUE))

#n.g <- NCOL(ms.genus.top)
n.g <- NCOL(ms.no0)

lambda.vec <- exp(seq(-2, 1, length.out=8))

# Running model-based variation on MS samples
ms.genus.ms <- ms.no0[covars$Diagnosis=="MS",]
mle.ms <- mlePath(ms.genus.ms, tol = 1e-4, tol.nr = 1e-4, lambda.gl=lambda.vec, n.cores=8)
save(mle.ms, file="/home/kdm147/research/heather_armstrong/output_2024/ms_pedes_mle_zi_ms_ebic.RData")

# Running model-based variation on control samples
ms.genus.control <- ms.no0[covars$Diagnosis=="Controls",]
mle.control <- mlePath(ms.genus.control, tol = 1e-4, tol.nr = 1e-4, lambda.gl=lambda.vec, n.cores=8)
save(mle.control, file="/home/kdm147/research/heather_armstrong/output_2024/ms_pedes_mle_zi_control_ebic.RData")