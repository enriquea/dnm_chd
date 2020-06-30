library(dplyr)
library(readr)
library(tidyverse)
library(IHW)

## Combine p-values and apply IHW with LOEUF contraint gene metrics.

## @author: Enrique Audain
## @date: 022020

### read/parse the gene lof metrics file (gnomad)
gene_lof <- read.table(file = "resources/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz",
                       header = T,
                       sep = '\t',
                       stringsAsFactors = F)

gene_lof %>% 
  select(gene,
         gene_id,
         chromosome,
         start_position,
         end_position,
         oe_lof,
         oe_lof_lower,
         oe_lof_upper,
         oe_lof_upper_bin) -> gene_lof

### read/parse the results from denovoWEST
dnw_files <- list.files("analysis/denovowest/run_10122019", 
                        full.names = T, 
                        pattern = "*.txt")

for(file in dnw_files){
  dnw <- read_delim(file, delim = '\t')
  dnw <- rbind(dnw)
}

dnw %>% 
    rename(gene = symbol, p_dnw = `p-value`) %>%
    select(gene, p_dnw) -> dnw


### read/parse the results from mupit
dne_mupit <- read.table("analysis/sifrim_jin_mupit_results_10122019.tsv",
                        sep = '\t',
                        header = T,
                        stringsAsFactors = F)

dne_mupit %>% 
  select(hgnc, p_lof) %>%
  rename(gene = hgnc, p_mupit = p_lof) %>%
  filter(!p_mupit == 0) -> dne_mupit

### read/parse the results from CNV-deletion permutation test
cnv_file_path <- "data/gene-wise_CNV_del_permutation_maf_0.01_10122019.tsv"
cnv <- read.delim(cnv_file_path, stringsAsFactors = F)

cnv %>% 
  rename(p_cnv = p_value, gene = symbol) %>%
  select(gene, p_cnv) -> cnv

### merge set
list(gene_lof, cnv, dnw, dne_mupit) %>% 
    reduce(left_join, by = "gene") -> genes

### annotate missing p-value as 1
genes %>%
  mutate(p_cnv = if_else(is.na(p_cnv), 1, p_cnv),
         p_dnw = if_else(is.na(p_dnw), 1, p_dnw),
         p_mupit = if_else(is.na(p_mupit), 1, p_mupit)) -> genes

### compute minimal p-value for DNM enrichment test
genes %>% mutate(min_p = pmin(p_dnw, p_mupit)) -> genes

### compute Fisher combine p-value for DNM/CNV-DEL
genes['p_comb'] <- apply(genes[, c('p_cnv','min_p')], 
                         MARGIN = 1, 
                         FUN = function(x) metap::sumlog(c(x[1], x[2]))$p)


### Bonferroni correction with IHW

# remove entries where covariates is NA
genes %>% filter(!is.na(oe_lof_upper_bin)) -> genes

# apply IHW
ihw_res <- IHW::ihw(p_comb ~ as.factor(oe_lof_upper_bin),  
               data = genes, 
               alpha = 0.05, 
               adjustment_type = 'bonferroni', 
               m_groups=table(genes$oe_lof_upper_bin))

# save p-values
genes['p_comb_ihw_adj_bonferroni'] <- IHW::adj_pvalues(ihw_res)


### FDR correction with IHW

# apply IHW
ihw_res <- IHW::ihw(p_comb ~ as.factor(oe_lof_upper_bin),  
               data = genes, 
               alpha = 0.05, 
               adjustment_type = 'BH', 
               m_groups=table(genes$oe_lof_upper_bin))

# save p-values
genes['p_comb_ihw_adj_fdr'] <- IHW::adj_pvalues(ihw_res)

# export results
write_delim(genes,
            'analysis/ihw_results_meta_analysis_dnm_cnv_15122019.tsv',
            delim = '\t')

