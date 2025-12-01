# Running MLE on cell cycle associated genes
library(ggVennDiagram)
library(ggplot2)
library(kableExtra)
library(ggpubr)


source("proportionality_code/mle_multinom_logit/variation_functions.R")
source("proportionality_code/mle_multinom_logit/helper_functions.R")
source("proportionality_code/mle_multinom_logit/mle_multinom_logit_normal.R")

#setwd("/Users/kevin/Dropbox/Documents/research/code/variation_estimation/single_cell")

# Load in data
load("data/single_cell/cleaned_data_all_stages.RData")
load("data/single_cell/mle_g1_GO_filter0.05.RData")
load("data/single_cell/mle_g2m_GO_filter0.05.RData")
load("data/single_cell/mle_s_GO_filter_0.8to1.6.RData")

g.names <- colnames(sc.g1)

# Consolidating rho
rho.all.cycles <- list(g1=model.est.g1$rho[-n.gene,-n.gene],
                       g2m=model.est.g2m$rho[-n.gene,-n.gene],
                       s=model.est.s$rho[-n.gene,-n.gene])
lt <- lower.tri(rho.all.cycles$g1)

# Getting naive metrics
rho.all.cycles$g1.naive <- runNaive(sc.g1)$rho
rho.all.cycles$g2m.naive <- runNaive(sc.g2m)$rho
rho.all.cycles$s.naive <- runNaive(sc.s)$rho

# Rho differences between cell cycles
rho.diff <- list()
rho.diff$g1.g2m <- rho.all.cycles$g1 - rho.all.cycles$g2m
rho.diff$g1.s <- rho.all.cycles$g1 - rho.all.cycles$s
rho.diff$g2m.s <- rho.all.cycles$g2m - rho.all.cycles$s
rho.diff$g1.g2m.naive <- rho.all.cycles$g1.naive - rho.all.cycles$g2m.naive
rho.diff$g1.s.naive <- rho.all.cycles$g1.naive - rho.all.cycles$s.naive
rho.diff$g2m.s.naive <- rho.all.cycles$g2m.naive - rho.all.cycles$s.naive

# Creating list of difference pairs to feed into mapply
r.pairs <- list(list(rho.all.cycles$g1, rho.all.cycles$g2m),
                list(rho.all.cycles$g1, rho.all.cycles$s),
                list(rho.all.cycles$g2m, rho.all.cycles$s),
                list(rho.all.cycles$g1.naive, rho.all.cycles$g2m.naive),
                list(rho.all.cycles$g1.naive, rho.all.cycles$s.naive),
                list(rho.all.cycles$g2m.naive, rho.all.cycles$s.naive))

n.top <- 100
n.cycle <- length(rho.all.cycles)
is.top <- is.top.diff <- vector("list", length=n.cycle)
names(is.top) <- names(rho.all.cycles)
names(is.top.diff) <- names(rho.diff)
for (i in 1:n.cycle) {
  cat("i =", i, "/", n.cycle, "\n")
  
  # Top for rho
  abs.rho <- abs(rho.all.cycles[[i]][lt])
  sort.vals <- sort(abs.rho, decreasing = TRUE)
  ct <- sort.vals[n.top]
  is.top[[i]] <- (abs(rho.all.cycles[[i]])>=ct)
  
  # Top for rho difference
  abs.rho.diff <- abs(rho.diff[[i]][lt])
  sort.vals.diff <- sort(abs.rho.diff, decreasing = TRUE)
  ct.diff <- sort.vals.diff[n.top]
  is.top.diff[[i]] <- (abs(rho.diff[[i]])>=ct.diff)
}

# Function to get names of top rho values
getNamesTop <- function(it, r, r.pair=NULL) {
  diff = !is.null(r.pair)
  if (diff) {
    top <- data.frame(gene1=rep("", n.top),
                      gene2=rep("", n.top),
                      rho1=rep(0, n.top),
                      rho2=rep(0, n.top),
                      rho=rep(0, n.top))
  } else {
    top <- data.frame(gene1=rep("", n.top),
                      gene2=rep("", n.top),
                      rho=rep(0, n.top))
  }
  
  count <- 1
  for (i in 1:(n.gene-1)) {
    for (j in 1:i) {
      if (it[i,j] & i!=j) {
        top$gene1[count] <- g.names[i]
        top$gene2[count] <- g.names[j]
        top$rho[count] <- r[i,j]
        
        if (diff) {
          top$rho1[count] <- r.pair[[1]][i,j]
          top$rho2[count] <- r.pair[[2]][i,j]
        }
        
        count <- count + 1
      }
    }
  }
  return(top)
}

top.rho <- mapply(getNamesTop, is.top, rho.all.cycles, SIMPLIFY = FALSE)
top.rho.diff <- mapply(getNamesTop, is.top.diff, rho.diff, r.pair=r.pairs,
                       SIMPLIFY = FALSE)

# Intersection of model based and naive
findIntersect <- function(s1, s2) {
  g.pairs1 <- paste(s1$gene1, s1$gene2)
  g.pairs2 <- paste(s2$gene1, s2$gene2)
  int <- intersect(g.pairs1, g.pairs2)
  return(list(gp1=g.pairs1, gp2=g.pairs2, int=int))
}

# Intersection for top rho values (model vs naive)
int <- list()
int$g1 <- findIntersect(top.rho$g1, top.rho$g1.naive)
int$g2m <- findIntersect(top.rho$g2m, top.rho$g2m.naive)
int$s <- findIntersect(top.rho$s, top.rho$s.naive)

# Intersection for top differences (model vs naive)
g1.g2m.int <- findIntersect(top.rho.diff$g1.g2m, top.rho$g1.g2m.naive)
g1.s.int <- findIntersect(top.rho.diff$g1.s, top.rho$g1.s.naive)
g2m.s.int <- findIntersect(top.rho.diff$g2m.s, top.rho.diff$g2m.s.naive)

# Venn diagrams within each cell cycle
plot.dir <- "proportionality_code/single_cell/plots/"
for (i in 1:(n.cycle/2)) {
  cyc <- names(int)[i]
  ti <- paste0("Top 100 gene pairs: ", toupper(cyc))
  p <- ggVennDiagram(int[[i]][-3], 
      category.names = c("Model based", "Empirical"), set_size=8,
      label_color = "lightblue") + 
    theme(legend.position = "none") + labs(title=ti) +
    theme(plot.title = element_text(size=22, hjust=0.5)) +
    scale_fill_distiller(palette = "Reds", direction = 1) +
    coord_flip() # Horizontal Venn diagram
  p.name <- paste0(plot.dir, "/venn_", cyc, "_model_vs_naive.pdf")
  ggsave(p.name, p, device="pdf", height=8, width=12)
}

# Adding associated gene names
addAssoc <- function(df, diff=FALSE) {
  an1 <- annotations[df$gene1,]$AssociatedGeneName
  an2 <- annotations[df$gene2,]$AssociatedGeneName
  if (diff) {
    ret <- data.frame(gene1=df$gene1, an1=an1, gene2=df$gene2, an2=an2,
               rho1=round(df$rho1, 4),
               rho2=round(df$rho2, 4),
               rho=round(df$rho, 4)) 
  } else {
    ret <-data.frame(gene1=df$gene1, an1=an1, gene2=df$gene2, an2=an2,
               rho=round(df$rho, 4)) 
  }
  return(ret)
}

top.rho.out <- lapply(top.rho, addAssoc)
top.rho.diff.out <- lapply(top.rho.diff, addAssoc, diff=TRUE)

# Top rho in each cell cycle
tab.dir <- "proportionality_code/single_cell/tables/"
for (i in 1:n.cycle) {
  out.tab <- top.rho.out[[i]]
  out.tab <- out.tab[order(out.tab[,5],decreasing = TRUE),]
  cyc <- toupper(names(top.rho.out)[i])
  names(out.tab) <- c("ID 1", "Assoc. gene 1", "ID 2", "Assoc. gene 2", 
                      paste0("$\\rho_{", cyc, "}$"))
  rownames(out.tab) <- NULL
  out.tab[,2] <- toupper(out.tab[,2])
  out.tab[,4] <- toupper(out.tab[,4])
  spec.tab <- out.tab
  spec.tab[5] <- cell_spec(spec.tab[[5]], bold = T, 
                           color = spec_color(spec.tab[[5]], end = 0.9,
                                              scale_from=c(-0.7,0.7)),
                           format="latex")
  t <- kbl(spec.tab, escape = F, align = c("|c","c","c","c","c|"), 
           format = "latex") %>%
    kable_classic("striped", full_width = F)
  f.name <- paste0(tab.dir, "/top_rho_", cyc, ".txt")
  write(t, f.name)
}


# Top differences between cell cycle pairs
tab.dir <- "proportionality_code/single_cell/tables/"
comp.names <- list(c("G1", "G2M"),
                   c("G1", "S"),
                   c("G2M", "S"))
for (i in 1:length(comp.names)) {
  out.tab <- top.rho.diff.out[[i]]
  out.tab <- out.tab[order(out.tab[,7],decreasing = TRUE),]
  names(out.tab) <- c("ID 1", "Assoc. gene 1", "ID 2", "Assoc. gene 2",
      paste0("$\\rho_{", comp.names[[i]][1], "}$"),
      paste0("$\\rho_{", comp.names[[i]][2], "}$"),
      paste0("$\\rho_{", comp.names[[i]][1], "}-\\rho_{", comp.names[[i]][2], "}$"))
  # "$\\rho_{1}-\\rho_{2}$"
  rownames(out.tab) <- NULL
  out.tab[,2] <- toupper(out.tab[,2])
  out.tab[,4] <- toupper(out.tab[,4])
  spec.tab <- out.tab
  spec.tab[5:7] <- lapply(spec.tab[5:7], function(x){
                      cell_spec(x, bold = T, 
                      color = spec_color(x, end = 0.9, scale_from=c(-0.7,0.7)),
                      format="latex")
                   })
  t <- kbl(spec.tab, escape = F, align = c("|c","c","c","c","c","c","c|"), format = "latex") %>%
    kable_classic("striped", full_width = F)
  f.name <- paste0(tab.dir, "/diff_", paste0(comp.names[[i]][1],comp.names[[i]][2]), 
                   ".txt")
  write(t, f.name)
}

# Subset table of top differences in cell cycle pairs
sub.dir <- "proportionality_code/single_cell/tables/subset/"
n.ss <- 25
for (i in 1:length(comp.names)) {
  out.tab <- top.rho.diff.out[[i]][,-c(1,3)]
  ct <- sort(abs(out.tab[,5]), decreasing = TRUE)[n.ss]
  out.tab <- out.tab[abs(out.tab[,5])>=ct,]
  out.tab <- out.tab[order(out.tab[,5],decreasing = TRUE),]
  names(out.tab) <- c("Gene 1", "Gene 2",
                      paste0("$\\rho_{", comp.names[[i]][1], "}$"),
                      paste0("$\\rho_{", comp.names[[i]][2], "}$"),
                      paste0("$\\rho_{", comp.names[[i]][1], "}-\\rho_{", comp.names[[i]][2], "}$"))
  rownames(out.tab) <- NULL
  out.tab[,1] <- toupper(out.tab[,1])
  out.tab[,2] <- toupper(out.tab[,2])
  spec.tab <- out.tab
  spec.tab[3:5] <- lapply(spec.tab[3:5], function(x){
                    cell_spec(x, bold = T, 
                    color = spec_color(x, end = 0.9, scale_from=c(-0.7,0.7)),
                    format="latex")
                 })
  t <- kbl(spec.tab, escape = F, align = c("|c","c","c","c","c|"), format = "latex") %>%
    kable_classic("striped", full_width = F)
  f.name <- paste0(sub.dir, "/diff_subset_", paste0(comp.names[[i]][1],comp.names[[i]][2]), 
                   ".txt")
  write(t, f.name)
}


# Comparing all values of model vs empirical estimates
# of rho in cell stages.

stg <- c("G1", "G2M", "S")
rho.plot <- vector("list", length(stg))
names(rho.plot) <- stg
count <- 1
for (s in stg) {
  df.rho <- data.frame(naive=rho.all.cycles[[count+3]][lt], 
                     model=rho.all.cycles[[count]][lt])
  rho.plot[[s]] <- ggplot(df.rho, aes(x=naive, y=model, color=col.scale)) +
    geom_point(colour="#009933", show.legend=FALSE, size=0.25) +
    ggtitle(s) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5)) +
    #scale_colour_gradient(low="green", high="red", breaks=brks,
    #                      labels=round(brk.vals, 6)) + 
    geom_abline(intercept=0, slope=1) +
    geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.1) +
    geom_vline(xintercept=0, linetype="dashed", color = "black", size=0.1) +
    xlab(expression(rho~(Empirical))) +
    ylab(expression(rho~(Model-based)))
  count <- count+1
}

combined.rho.plots <- ggarrange(rho.plot[[1]], rho.plot[[2]], rho.plot[[3]],
          ncol=3, nrow=1)

ggsave("proportionality_code/single_cell/plots/rho_plot_single_cell.pdf", 
       combined.rho.plots, height=6, width=12)



