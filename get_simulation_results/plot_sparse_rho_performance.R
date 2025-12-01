# Make plots for variation matrix estimation accuracy
library(ggplot2)
library(gridExtra)
library(ggpubr)

plot.dir <- "/Users/kevin/Dropbox/Documents/research/code/variation_estimation/single_cell/plot_sparse_rho_performance/plots/"

# Load in results from Graham/Cedar clusters
load("/Users/kevin/Dropbox/Documents/research/code/variation_estimation/single_cell/from_graham/sim_results/merged_sim_sparse_rho_clr.RData")
f.est.clr <- f.est
f.naive.clr <- f.naive
mse.est.clr <- mse.est
mse.naive.clr <- mse.naive
load("/Users/kevin/Dropbox/Documents/research/code/variation_estimation/single_cell/from_graham/sim_results/merged_sim_sparse_rho.RData")


# Specifying sparsity probs, sample sizes, and numbers of features
nsp <- factor(rep(rep(paste0("Spar=",c(0.15,0.5,0.85)), each=9),2))
ns <- factor(rep(rep(paste0("n=", c(100,200,500)), 9), 2))
nf <- factor(rep(rep(rep(paste0("Features=", c(50,100,200)), each=3), 3),
            levels=paste0("Features=", c(50,100,200)),2))
lr <- factor(rep(paste0("LR=",c("ALR", "CLR")), each=9*3*2))

# List containing metrics to plot for
metric <- "f.est"
# n.metric <- length(metrics)
# r.list.naive <- list(aupr.naive)
# r.list <- list(aupr.est)

# Binding results to data frame for ggplot
# df <- data.frame(val=c(aupr.naive, aupr.est),
#                  method=rep(c("Empirical", "Model-based"), each=n.rep*n.run),
#                  spar=rep(nsp, each=n.rep),
#                  n=rep(ns, each=n.rep),
#                  J=rep(nf, each=n.rep))

# # Plotting results
# x.ax <- ""
# y.ax <- "AUPR"
# title <- ""
# aupr_plot <- ggplot(df, aes(n, val, fill=method)) + geom_boxplot() +
#   facet_grid(J~spar, scales = "free_y") +
#   theme_bw() +
#   #theme(legend.position = "bottom",text = element_text(size=20), axis.text.x = element_text(size=16,angle=0,vjust=1,hjust=0.5)) +
#   theme(legend.position = "none", plot.title = element_text(size=16, hjust = 0.5), strip.background = element_rect(
#     color="black", fill="skyblue2", size=1, linetype="solid"),
#     strip.text = element_text(size = 12),
#     axis.text=element_text(size=14),
#     axis.title = element_text(size = 14)) +
#   labs(title=title) +
#   #ylim(0,5) +
#   #geom_hline(yintercept = 0, color="darkblue", linetype=2) +
#   xlab(x.ax) + ylab(y.ax) + cowplot::panel_border() + scale_fill_brewer(palette="Set1")
# 
# name <- paste0(plot.dir, "/perf_", metric, ".pdf")
# ggsave(name, aupr_plot, device = "pdf", height = 6, width=10)
# 


# Barplots... but they don't seem good for comparison in this data

# Calculating summary stats of aupr for each set of simulation scenarios
# med <- c(apply(aupr.naive, 2, median),
#          apply(aupr.plugin, 2, median),
#          apply(aupr.est, 2, median))
# low <- c(apply(aupr.naive, 2, function(x){quantile(x,probs=0.25,na.rm=TRUE)}),
#          apply(aupr.plugin, 2, function(x){quantile(x,probs=0.25,na.rm=TRUE)}),
#          apply(aupr.est, 2, function(x){quantile(x,probs=0.25,na.rm=TRUE)}))
# high <- c(apply(aupr.naive, 2, function(x){quantile(x,probs=0.75,na.rm=TRUE)}),
#           apply(aupr.plugin, 2, function(x){quantile(x,probs=0.75,na.rm=TRUE)}),
#           apply(aupr.est, 2, function(x){quantile(x,probs=0.75,na.rm=TRUE)}))

med <- c(apply(f.naive, 2, median, na.rm=TRUE),
         apply(f.est, 2, median, na.rm=TRUE),
         apply(f.naive.clr, 2, median, na.rm=TRUE),
         apply(f.est.clr, 2, median, na.rm=TRUE))
low <- c(apply(f.naive, 2, function(x){quantile(x,probs=0.25,na.rm=TRUE)}),
         apply(f.est, 2, function(x){quantile(x,probs=0.25,na.rm=TRUE)}),
         apply(f.naive.clr, 2, function(x){quantile(x,probs=0.25,na.rm=TRUE)}),
         apply(f.est.clr, 2, function(x){quantile(x,probs=0.25,na.rm=TRUE)}))
high <- c(apply(f.naive, 2, function(x){quantile(x,probs=0.75,na.rm=TRUE)}),
          apply(f.est, 2, function(x){quantile(x,probs=0.75,na.rm=TRUE)}),
          apply(f.naive.clr, 2, function(x){quantile(x,probs=0.75,na.rm=TRUE)}),
          apply(f.est.clr, 2, function(x){quantile(x,probs=0.75,na.rm=TRUE)}))

med.mse <- c(apply(mse.naive, 2, median, na.rm=TRUE),
         apply(mse.est, 2, median, na.rm=TRUE),
         apply(mse.naive.clr, 2, median, na.rm=TRUE),
         apply(mse.est.clr, 2, median, na.rm=TRUE))
low.mse <- c(apply(mse.naive, 2, function(x){quantile(x,probs=0.25,na.rm=TRUE)}),
         apply(mse.est, 2, function(x){quantile(x,probs=0.25,na.rm=TRUE)}),
         apply(mse.naive.clr, 2, function(x){quantile(x,probs=0.25,na.rm=TRUE)}),
         apply(mse.est.clr, 2, function(x){quantile(x,probs=0.25,na.rm=TRUE)}))
high.mse <- c(apply(mse.naive, 2, function(x){quantile(x,probs=0.75,na.rm=TRUE)}),
          apply(mse.est, 2, function(x){quantile(x,probs=0.75,na.rm=TRUE)}),
          apply(mse.naive.clr, 2, function(x){quantile(x,probs=0.75,na.rm=TRUE)}),
          apply(mse.est.clr, 2, function(x){quantile(x,probs=0.75,na.rm=TRUE)}))


# df <- data.frame(med=med, low=low, high=high,
#                  spar=rep(nsp,3), n=rep(ns,3), J=rep(nf,3),
#                  method=rep(c("Empirical", "Plugin", "Model-based"), each=n.run))

df <- data.frame(med=med, low=low, high=high, 
                spar=rep(nsp,2), n=rep(ns,2), J=rep(nf,2),
                Method=rep(rep(c("Empirical", "Model-based"), each=n.run), 2),
                lr=rep(lr, 2))
df.mse <- data.frame(med=med.mse, low=low.mse, high=high.mse, 
                 spar=rep(nsp,2), n=rep(ns,2), J=rep(nf,2),
                 Method=rep(rep(c("Empirical", "Model-based"), each=n.run), 2),
                 lr=rep(lr, 2))


#df$method <- factor(df$method, levels=c("Empirical", "Plugin", "Model-based"))
df$Method <- factor(df$Method, levels=c("Empirical", "Model-based"))
df$J <- factor(df$J, levels=paste0("Features=", c(50,100,200)))

df.mse$Method <- factor(df.mse$Method, levels=c("Empirical", "Model-based"))
df.mse$J <- factor(df.mse$J, levels=paste0("Features=", c(50,100,200)))


# Plotting results for ALR
x.ax <- ""
y.ax <- "F-Score"
title <- "ALR"
alr_plot <- ggplot(df[df$lr=="LR=ALR",], aes(n, med, fill=Method)) +
  geom_bar(stat="identity", position=position_dodge(0.4), width = 0.4) +
  geom_errorbar(aes(ymin=low, ymax=high), width=.2,
                position=position_dodge(0.4)) +
  facet_grid(J~spar, scales = "free_y") +
  theme_bw() +
  #theme(legend.position = "bottom",text = element_text(size=20), axis.text.x = element_text(size=16,angle=0,vjust=1,hjust=0.5)) +
  theme(legend.position = "bottom", plot.title = element_text(size=16, hjust = 0.5), strip.background = element_rect(
    color="black", fill="skyblue2", size=1, linetype="solid"),
    strip.text = element_text(size = 12),
    axis.text=element_text(size=14),
    axis.title = element_text(size = 18)) +
  #labs(title=title) +
  #ylim(0.5,1) +
  #geom_hline(yintercept = 0, color="darkblue", linetype=2) +
  xlab(x.ax) + ylab(y.ax) + cowplot::panel_border() + scale_fill_brewer(palette="Set2")

y.ax.mse <- "MSE"
alr_mse_plot <- ggplot(df.mse[df.mse$lr=="LR=ALR",], aes(n, med, fill=Method)) +
  geom_bar(stat="identity", position=position_dodge(0.4), width = 0.4) +
  geom_errorbar(aes(ymin=low, ymax=high), width=.2,
                position=position_dodge(0.4)) +
  facet_grid(J~spar, scales = "free_y") +
  theme_bw() +
  #theme(legend.position = "bottom",text = element_text(size=20), axis.text.x = element_text(size=16,angle=0,vjust=1,hjust=0.5)) +
  theme(legend.position = "bottom", plot.title = element_text(size=16, hjust = 0.5), strip.background = element_rect(
    color="black", fill="skyblue2", size=1, linetype="solid"),
    strip.text = element_text(size = 12),
    axis.text=element_text(size=14),
    axis.title = element_text(size = 18)) +
  #labs(title=title) +
  #ylim(0.5,1) +
  #geom_hline(yintercept = 0, color="darkblue", linetype=2) +
  xlab(x.ax) + ylab(y.ax.mse) + cowplot::panel_border() + scale_fill_brewer(palette="Set2")

alr_combined <- ggarrange(alr_plot, alr_mse_plot, ncol=2, nrow=1, 
          common.legend = TRUE, legend="bottom")
alr_combined_n <- annotate_figure(alr_combined, 
          top = text_grob("ALR Performance Metrics", face = "bold", size = 14))

name <- paste0(plot.dir, "/perf_alr.pdf")
ggsave(name, alr_combined_n, device = "pdf", height = 8, width=16)



# Plotting results for CLR
x.ax <- ""
y.ax <- "F-Score"
title <- "CLR"
clr_plot <- ggplot(df[df$lr=="LR=CLR",], aes(n, med, fill=Method)) +
  geom_bar(stat="identity", position=position_dodge(0.4), width = 0.4) +
  geom_errorbar(aes(ymin=low, ymax=high), width=.2,
                position=position_dodge(0.4)) +
  facet_grid(J~spar, scales = "free_y") +
  theme_bw() +
  #theme(legend.position = "bottom",text = element_text(size=20), axis.text.x = element_text(size=16,angle=0,vjust=1,hjust=0.5)) +
  theme(legend.position = "bottom", plot.title = element_text(size=16, hjust = 0.5), strip.background = element_rect(
    color="black", fill="skyblue2", size=1, linetype="solid"),
    strip.text = element_text(size = 12),
    axis.text=element_text(size=14),
    axis.title = element_text(size = 18)) +
  #labs(title=title) +
  #ylim(0.5,1) +
  #geom_hline(yintercept = 0, color="darkblue", linetype=2) +
  xlab(x.ax) + ylab(y.ax) + cowplot::panel_border() + scale_fill_brewer(palette="Set2")

y.ax.mse <- "MSE"
clr_mse_plot <- ggplot(df.mse[df.mse$lr=="LR=CLR",], aes(n, med, fill=Method)) +
  geom_bar(stat="identity", position=position_dodge(0.4), width = 0.4) +
  geom_errorbar(aes(ymin=low, ymax=high), width=.2,
                position=position_dodge(0.4)) +
  facet_grid(J~spar, scales = "free_y") +
  theme_bw() +
  #theme(legend.position = "bottom",text = element_text(size=20), axis.text.x = element_text(size=16,angle=0,vjust=1,hjust=0.5)) +
  theme(legend.position = "bottom", plot.title = element_text(size=16, hjust = 0.5), strip.background = element_rect(
    color="black", fill="skyblue2", size=1, linetype="solid"),
    strip.text = element_text(size = 12),
    axis.text=element_text(size=14),
    axis.title = element_text(size = 18)) +
  #labs(title=title) +
  #ylim(0.5,1) +
  #geom_hline(yintercept = 0, color="darkblue", linetype=2) +
  xlab(x.ax) + ylab(y.ax.mse) + cowplot::panel_border() + scale_fill_brewer(palette="Set2")

clr_combined <- ggarrange(clr_plot, clr_mse_plot, ncol=2, nrow=1, 
                          common.legend = TRUE, legend="bottom")
clr_combined_n <- annotate_figure(clr_combined, 
          top = text_grob("CLR Performance Metrics", face = "bold", size = 14))

name <- paste0(plot.dir, "/perf_clr.pdf")
ggsave(name, clr_combined_n, device = "pdf", height = 8, width=16)

