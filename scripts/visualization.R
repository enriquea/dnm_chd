library(ggplot2)
library(ggrepel)
library(ggpubr)
library(dplyr)

# read results from IHW
ihw_results <- read.delim("data/ihw_results_meta_analysis_dnm_cnv_15122019.tsv", 
                          stringsAsFactors=FALSE)

#-------------------------------------------------------------------------------
# 1. manhattan plot with ggplot2 and ggrepel for annotation
#-------------------------------------------------------------------------------

dat <- ihw_results

# prepare data for plotting
dat %>%
  filter(chromosome %in% 1:22,
         p_comb < 1) -> dat

dat$chromosome <- factor(dat$chromosome, levels = 1:22)
dat$logP <- -log10(dat$p_comb)

pos <- position_jitter(w = 0.4, h = 0.4, seed = 123)

y_threshold <- max(subset(dat, p_comb_adj < 0.05)$p_comb)

# pos <- position_jitter(width = 0.5, seed = 1)
p <- ggplot(dat, aes(chromosome, logP, colour = as.factor(chromosome)), label = gene) +
  geom_point(pos = pos, size=2.0, alpha=.8) +
  geom_hline(yintercept = -log10(y_threshold) - 0.3, 
             colour='red', alpha=.5, linetype='dashed') +
  scale_x_discrete(limits=c(1:22)) +
  scale_color_manual(breaks = seq(1,22,1), 
                     values= rep(c('deepskyblue2', 'darkblue'), 11)) +
  ylim(c(0.1, 30)) +
  theme_bw() +
  xlab('Chromosome') +
  ylab(expression(~-log[10](italic(p)))) +
  theme(panel.border = element_rect(fill = NA, colour = "black", size = 0.2), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position='none',
        axis.text=element_text(size=12),
        axis.title=element_text(size=14))

p <- p + ggrepel::geom_text_repel(data = subset(dat, p_comb_adj < 0.05), 
                                  aes(label = gene),
                                  size = 2.5,
                                  direction = 'both',
                                  segment.size = 0.2,
                                  ylim = c(-log10(y_threshold) + 0.2, NA))

# export plot
ggsave(plot = p,
       filename = 'manhattan_plot_Pcomb_06012020.pdf',
       device = 'pdf',
       width = 10,
       height = 6)

#-------------------------------------------------------------------------------
# 2. quantile expression vs. enrichment (P)
#-------------------------------------------------------------------------------

heg <- read.delim("data/heart_gene_mean_expression_percentile_13112019.tsv")


# select fields
heg %>% select(gene, rank_pct_expression, dev_stage) -> heg

# merge with gene p-values
df <- merge(ihw_results, heg, by = 'gene')

# transform/sort dev_stage levels
df %>%
  mutate(dev_stage=as.character(dev_stage)) %>%
  mutate(dev_stage=ifelse(dev_stage=='infant','infant/adult', dev_stage)) %>%
  mutate(dev_stage=factor(dev_stage, levels = c('development','maturation','infant/adult'))) -> df

p <- ggplot(subset(df, p_comb_adj < 0.05), 
            aes(x=rank_pct_expression, 
                y=oe_lof_upper, 
                color = dev_stage, 
                label=gene)) + 
       geom_point(alpha=0.8) +
       geom_vline(xintercept = 0.75, linetype=2, size=0.3) +
       geom_hline(yintercept = 0.30, linetype=2, size=0.3) +
       xlab('Percentile expression') +
       ylab('LOEUF') +
       scale_y_reverse() +
       theme_bw() +
       theme(panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             legend.position = "none") +
       ggrepel::geom_text_repel(size=2, segment.size = 0.2) + 
       facet_grid(dev_stage ~ .)


# export plot
ggsave(plot = p,
       filename = 'sig_genes_vs_heart_expression_06012020.pdf',
       device = 'pdf',
       width = 6,
       height = 10)


#-------------------------------------------------------------------------------
# 3. Boxplot p-value cutoff (Pcomb) vs. expression on the heart
#-------------------------------------------------------------------------------

heg <- read.delim("data/heart_gene_mean_expression_percentile_13112019.tsv")


# select fields
heg %>% select(gene, mean_expression, dev_stage) -> heg

# merge with gene p-values
heg <- merge(heg, ihw_results, by = 'gene')

# filter out genes with rpkm < 1 and convert to log
heg$logME <- log2(heg$mean_expression)

# annotate p-value bins
heg %>% 
  dplyr::mutate(pval_bins=case_when(p_comb < 1e-5 ~ "<1e-5",
                                    p_comb < 1e-3 ~ "<1e-3",
                                    p_comb < 0.05 ~ "<0.05",
                                    p_comb >= 0.05 ~ ">=0.05")) -> df

# reshape dataframe to include "all genes" as a independent group 
df <- subset(df, select = c("pval_bins", "logME", "dev_stage"))
df_copy <- df
df_copy$pval_bins <- 'all genes'
df <- rbind(df, df_copy)

# sort p-value bins levels
df$pval_bins <- factor(df$pval_bins, levels = c("<1e-5",
                                                "<1e-3",
                                                "<0.05",
                                                ">=0.05",
                                                "all genes"))

# transform/sort dev_stage levels
df %>%
  mutate(dev_stage=as.character(dev_stage)) %>%
  mutate(dev_stage=ifelse(dev_stage=='infant','infant/adult', dev_stage)) %>%
  mutate(dev_stage=factor(dev_stage, levels = c('development','maturation','infant/adult'))) -> df

# boxplot with pairwise comparision against 'all genes' (Wilcox test)
p1 <- ggboxplot(data = df,
               x = 'pval_bins',
               y = 'logME',
               fill = 'pval_bins',
               legend = 'none',
               alpha=0.7) +
  rotate_x_text(angle = 45) +
  xlab('p-value cutoff') +
  ylab('log2(mean expression)') +
  stat_compare_means(label = "p.signif",
                     label.y = 15.5,
                     method = "wilcox.test",
                     ref.group = "all genes") # pairwise comparison against 'all genes'

# split plot by development stage
p1 <- p1 + ggplot2::facet_wrap(. ~ dev_stage)

# export plot
ggsave(plot = p1,
       filename = 'genes_pcomb_stratified_heart_expression_06012020.pdf',
       device = 'pdf',
       width = 10,
       height = 5)



#-------------------------------------------------------------------------------
# 4. Boxplot p-value cutoff (Pdnm) vs. LOEUF
#-------------------------------------------------------------------------------

# plot DNM enrichment

# annotate p-value bins
ihw_results %>% 
  dplyr::mutate(pval_bins=case_when(min_p < 1e-6 ~ "<1e-6",
                                    min_p >= 1e-6 & min_p < 1e-3 ~ "[1e-6,1e-3)",
                                    min_p >= 0.001 & min_p < 0.05 ~ "[0.001,0.05)",
                                    min_p >= 0.05 & min_p <= 1 ~ "[0.05,1]")) -> df

# reshape dataframe to include all genes as group 
df <- subset(df, select = c("pval_bins", "oe_lof_upper"))
df_copy <- df
df_copy$pval_bins <- 'all genes'
df <- rbind(df, df_copy)

# sort p-value bins levels
df$pval_bins <- factor(df$pval_bins, levels = c("<1e-6",
                                                "[1e-6,1e-3)",
                                                "[0.001,0.05)",
                                                "[0.05,1]",
                                                "all genes"))

# boxplot with pairwise comparision against 'all genes' (Wilcox test)
p <- ggboxplot(data = df,
               x = 'pval_bins',
               y = 'oe_lof_upper',
               fill = 'pval_bins',
               legend = 'none') +
  rotate_x_text(angle = 45) +
  xlab('p-value bins') +
  ylab('LOEUF') +
  stat_compare_means(label = "p.signif",
                     label.y = 2.1,
                     method = "wilcox.test",
                     ref.group = "all genes") # pairwise comparison against 'all genes'

# export plot
ggsave(plot = p,
       filename = 'DNM_geneset_enrichment_LOEUF_05012020.pdf',
       device = 'pdf',
       width = 8,
       height = 6)



#-------------------------------------------------------------------------------
# 5. Sccatter plot p-value (Pcomb) vs. Heart-DEG (r2)
#-------------------------------------------------------------------------------


ddg <- read.delim("data/DEG_heart_vs_kidney_bonferroni_0.01_03012020.tsv")

# merge with gene p-values
ddg <- merge(ddg, ihw_results, by = 'gene')

# reshape dataframe to include "all genes" as a independent group 
ddg %>% 
  select("gene", "R2", "p_value_adj", "p_comb", "p_comb_adj") %>%
  mutate(logP=-log10(p_comb),
         significant_genes=ifelse(p_comb_adj < 0.05, T, F),
         Heart_DEG=ifelse(R2 >= 0.50, T, F)) %>%
  filter(p_value_adj < 0.001) ->df
  


# scatter plot
p <- ggscatter(data = df,
               x = 'R2',
               y = 'logP',
               color = 'Heart_DEG',
               label = "gene",
               label.select = subset(df, significant_genes==T)$gene,
               repel = TRUE,
               font.label = list(color="Heart_DEG",size=8),
               legend="none",
               alpha=0.7) +
     geom_vline(xintercept = 0.50, linetype=2, size=0.3) +
     xlab('R2') +
     ylab('-log10(P)')

p <- ggpar(p, palette = "grey")

# export plot
ggsave(plot = p,
       filename = 'DEG_HeartKidney_vs_Pcomb_07012020.pdf',
       device = 'pdf',
       width = 8,
       height = 6)



#-------------------------------------------------------------------------------
# 6. Sccatter corrplot p-value (Pcomb) vs. LOEUF (r2)
#-------------------------------------------------------------------------------

ihw_results %>%
  filter(!is.na(p_comb) | !is.na(oe_lof_upper)) %>%
  mutate(logP = -log10(p_comb)) -> df

p <- ggscatter(df,
               x = "oe_lof_upper",
               y = "logP",
               xlab = "LOEUF",
               color = 'grey50',
               ylab = "-log10(P)",
               add = 'reg.line') + 
      stat_cor(label.x = 1)

# export plot
ggsave(plot = p,
       filename = 'corrplot_LOEUF_vs_Pcomb_07012020.pdf',
       device = 'pdf',
       width = 5,
       height = 5)

