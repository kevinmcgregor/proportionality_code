# Analyzing results of variation estimation in 
# MS pedes data (2024 version)
# VERSION WHERE ZERO INFLATION WAS DONE
# split up by disease status
# ANAYLSIS DONE ON GRAHAM
# Best lambda selected using EBIC

library(ggplot2)
library(ggVennDiagram)
library(compositions)
library(kableExtra)
source("~/Dropbox/Documents/research/code/variation_estimation/variation_functions.R")
source("~/Dropbox/Documents/research/code/variation_estimation/single_cell/R/helper_functions.R")
source("~/Dropbox/Documents/research/code/variation_estimation/mle_multinom_logit/mle_multinom_logit_normal.R")
source("/Users/kevin/Dropbox/Documents/research/collaborators/heather_armstrong/data/pedes_ms/R/helper_functions.R")

load("/Users/kevin/Dropbox/Documents/research/collaborators/heather_armstrong/data/pedes_ms_2024/cleaned_data/ms_pedes_cleaned.RData")
n.g <- NCOL(ms.genus)

# Load model-based variation on MS samples
ms.genus.ms <- ms.genus[covars$Diagnosis=="MS",]
load("/Users/kevin/Library/CloudStorage/Dropbox/Documents/research/collaborators/heather_armstrong/data/pedes_ms_2024/from_graham/ms_pedes_mle_zi_control_ebic.RData")

# Load model-based variation on control samples
ms.genus.control <- ms.genus[covars$Diagnosis=="Controls",]
load("/Users/kevin/Library/CloudStorage/Dropbox/Documents/research/collaborators/heather_armstrong/data/pedes_ms_2024/from_graham/ms_pedes_mle_zi_ms_ebic.RData")

# Parameter estimates
mu.ms <- mle.ms$est.min$mu
Sigma.ms <- mle.ms$est.min$Sigma
mu.control <- mle.control$est.min$mu
Sigma.control <- mle.control$est.min$Sigma

g.names <- colnames(ms.genus)

# MLE-based variation (ALR)
v.ms <- logitNormalVariation(mu.ms, Sigma.ms, lr = "alr")
phi.ms <- logitNormalVariation(mu.ms, Sigma.ms, type = "phi", lr = "alr")
phis.ms <- logitNormalVariation(mu.ms, Sigma.ms, type = "phis", lr = "alr")
rho.ms <- logitNormalVariation(mu.ms, Sigma.ms, type = "rho", lr = "alr")
v.c <- logitNormalVariation(mu.control, Sigma.control, lr = "alr")
phi.c <- logitNormalVariation(mu.control, Sigma.control, type = "phi", lr = "alr")
phis.c <- logitNormalVariation(mu.control, Sigma.control, type = "phis", lr = "alr")
rho.c <- logitNormalVariation(mu.control, Sigma.control, type = "rho", lr = "alr")

rownames(rho.ms) <- colnames(rho.ms) <- colnames(ms.genus.ms)
rownames(rho.c) <- colnames(rho.c) <- colnames(ms.genus.control)

# Imputing zeros for naive variation
ms.all.no.zero <- as.matrix(zCompositions::cmultRepl(ms.genus, output = "p-counts", z.warning = 0.9999,
                                                     suppress.print=TRUE))
ms.no0 <- ms.all.no.zero[covars$Diagnosis=="MS",]
control.no0 <- ms.all.no.zero[covars$Diagnosis=="Controls",]

# Naive variation
library(propr)
v.naive.ms <- naiveVariation(ms.no0, lr = "alr")
phi.naive.ms <- naiveVariation(ms.no0, type="phi", lr = "alr")
phis.naive.ms <- naiveVariation(ms.no0, type="phis", lr = "alr")
rho.naive.ms <- naiveVariation(ms.no0, type="rho", lr = "alr")
v.naive.c <- naiveVariation(control.no0, lr = "alr")
phi.naive.c <- naiveVariation(control.no0, type="phi", lr = "alr")
phis.naive.c <- naiveVariation(control.no0, type="phis", lr = "alr")
rho.naive.c <- naiveVariation(control.no0, type="rho", lr = "alr")

# Top differences MS vs control
cv <- 100
rho.diff <- rho.ms - rho.c
rho.diff.tvals <- getTopVals(rho.diff, 1)

rho.tvals <- getTopVals(rho.ms, rho.c, cut.val = 100)

# Getting list of top genus pairs
is.top.ms <- abs(rho.ms)>=rho.tvals$cutoff1
is.top.c <- abs(rho.c)>=rho.tvals$cutoff2

is.top.diff <- abs(rho.diff)>=rho.diff.tvals$cutoff1

top.ms <- top.c <- matrix("", cv, 2)
top.diff <- matrix("", cv, 2)
top.ms.idx <- top.c.idx <- matrix(0, cv, 2)
top.rho.ms <- top.rho.c <- top.rho.diff <- rep(0, cv)
top.rho.diff.both <- matrix(0, cv, 2)
count1 <- count2 <- count1.5 <- count2.5 <- count.diff <- 1
for (i in 1:(n.g-1)) {
  for (j in 1:i) {
    if (is.top.ms[i,j] & i!=j) {
      top.ms[count1,] <- c(g.names[i], g.names[j])
      top.rho.ms[count1] <- rho.ms[i,j]
      count1 <- count1 + 1
    }
    if (is.top.c[i,j] & i!=j) {
      top.c[count2,] <- c(g.names[i], g.names[j])
      top.rho.c[count2] <- rho.c[i,j]
      count2 <- count2 + 1
    }
    if (is.top.diff[i,j] & i!=j) {
      top.diff[count.diff,] <- c(g.names[i], g.names[j])
      top.rho.diff[count.diff] <- rho.diff[i,j]
      top.rho.diff.both[count.diff,] <- c(rho.ms[i,j], rho.c[i,j])
      count.diff <- count.diff + 1
    }
  }
}

table(c(top.ms))
table(c(top.c))

intersect(top.ms, top.c)

# Finding names of top pairs in common
ntopms <- paste0(top.ms[,1], top.ms[,2])
ntopc <- paste0(top.c[,1], top.c[,2])
top.common <- intersect(ntopms, ntopc)


# Making table of top findings in ms and controls
top.ms.vals <- top.c.vals <- top.diff.vals <- rep(0, cv)
for (i in 1:cv) {
  top.ms.vals[i] <- rho.ms[top.ms[i,1], top.ms[i,2]]
  top.c.vals[i] <- rho.c[top.c[i,1], top.c[i,2]]
  top.diff.vals[i] <- rho.diff[top.diff[i,1], top.diff[i,2]]
}


# Table for top differences MS - Control
diff.out.tab <- data.frame(gsub("_", "", top.diff),
                           round(top.rho.diff.both, 4),
                           round(top.diff.vals,4))
col.range <- range(diff.out.tab[,5])
diff.out.tab <- diff.out.tab[order(diff.out.tab[,5],decreasing = TRUE),]
names(diff.out.tab) <- c("Genus 1", "Genus 2", "$\\rho_{MS}$", "$\\rho_{C}$", "$\\rho_{MS}-\\rho_{C}$")
rownames(diff.out.tab) <- NULL
diff.spec.tab <- diff.out.tab
diff.spec.tab[3:5] <- lapply(diff.spec.tab[3:5], function(x){
  cell_spec(x, bold = T, 
            color = spec_color(x, end = 0.9,
                               scale_from=col.range),
            format="latex")
})
diff.t <- kbl(diff.spec.tab, escape = F, align = c("|c","c","c","c","c|"), format="latex") %>%
  kable_classic_2("striped", full_width = F)
write(diff.t, "/Users/kevin/Dropbox/Documents/research/collaborators/heather_armstrong/data/pedes_ms_2024/tables/diff_table_zero_imputed.txt")



# Table for top MS pairs
ms.out.tab <- data.frame(gsub("_", "", top.ms), round(top.ms.vals,4))
ms.out.tab <- ms.out.tab[order(ms.out.tab[,3],decreasing = TRUE),]
names(ms.out.tab) <- c("Genus 1", "Genus 2", "$\\rho_{MS}$")
rownames(ms.out.tab) <- NULL
ms.spec.tab <- ms.out.tab
ms.spec.tab[3] <- cell_spec(ms.spec.tab[[3]], bold = T, 
                             color = spec_color(ms.spec.tab[[3]], end = 0.9,
                                                scale_from=col.range),
                             format="latex")
ms.t <- kbl(ms.spec.tab, escape = F, align = c("|c","c","c|"), format="latex") %>%
  kable_classic("striped", full_width = F)
write(ms.t, "/Users/kevin/Dropbox/Documents/research/collaborators/heather_armstrong/data/pedes_ms_2024/tables/ms_table_zero_imputed.txt")


# Table for top control pairs
c.out.tab <- data.frame(gsub("_", "", top.c), round(top.c.vals,4))
c.out.tab <- c.out.tab[order(c.out.tab[,3],decreasing = TRUE),]
names(c.out.tab) <- c("Genus 1", "Genus 2", "$\\rho_{C}$")
rownames(c.out.tab) <- NULL
c.spec.tab <- c.out.tab
c.spec.tab[3] <- cell_spec(c.spec.tab[[3]], bold = T, 
                          color = spec_color(c.spec.tab[[3]], end = 0.9,
                                             scale_from=col.range),
                           format="latex")
c.t <- kbl(c.spec.tab, escape = F, align = c("|c","c","c|"), format="latex") %>%
  kable_classic("striped", full_width = F)
write(c.t, "/Users/kevin/Dropbox/Documents/research/collaborators/heather_armstrong/data/pedes_ms_2024/tables/control_table_zero_imputed.txt")



# Reduced table for top differences MS - Control (Reduced for main text)
n.incl <- 25
cut.diff <- sort(abs(diff.out.tab[,5]), decreasing = TRUE)[n.incl]
diff.out.reduced <- diff.out.tab[abs(diff.out.tab[,5])>=cut.diff,]
diff.out.reduced <- diff.out.reduced[order(diff.out.reduced[,5],decreasing = TRUE),]
rownames(diff.out.reduced) <- NULL
names(diff.out.reduced) <- c("Genus 1", "Genus 2", "$\\rho_{MS}$", "$\\rho_{C}$", "$\\rho_{MS}-\\rho_{C}$")
diff.spec.reduced <- diff.out.reduced
diff.spec.reduced[3:5] <- lapply(diff.spec.reduced[3:5], function(x){
  cell_spec(x, bold = T, 
            color = spec_color(x, end = 0.9, scale_from=col.range),
            format="latex")
})
diff.t.reduced <- kbl(diff.spec.reduced, escape = F, align = c("|c","c","c","c","c|"), format="latex") %>%
  kable_classic_2("striped", full_width = F)
write(diff.t.reduced, "/Users/kevin/Dropbox/Documents/research/collaborators/heather_armstrong/data/pedes_ms_2024/tables/diff_table_reduced_zero_imputed.txt")

save(diff.out.reduced, diff.out.tab, file="/Users/kevin/Dropbox/Documents/research/collaborators/heather_armstrong/data/pedes_ms_2024/tables/diff_zero_imputed_obj.RData")

