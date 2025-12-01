# Running MLE on cell cycle associated genes
source("~/research/variation/R/variation_functions.R")
source("~/research/variation/R/helper_functions.R")
source("~/research/variation/R/mle_multinom_logit/mle_multinom_logit_normal.R")

# Load in cleaned data
load("/project/def-kdm147/kdm147/all_cycles_single_cell/cleaned_data_all_stages.RData")

#model.est.g1 <- runModelBased(sc.g1, impute = FALSE)
#save(model.est.g1, file="data/mle_g1_GO_filter0.2.RData")

#model.est.g2m <- runModelBased(sc.g2m, impute = FALSE)
#save(model.est.g2m, file="data/mle_g2m_GO_filter0.2.RData")

model.est.s <- runModelBased(sc.s, impute = FALSE)
save(model.est.s, file="/project/def-kdm147/kdm147/all_cycles_single_cell/mle_s_GO_filter0.2.RData")

