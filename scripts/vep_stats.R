### Compute stats from VEP annotations

library(readr)
library(dplyr)

## Read VEP'ed dataset
vep_df <- read.table(file = 'data/dnm_jin_sifrim_veped_02122019.txt',
                     sep = '\t',
                     header = T,
                     stringsAsFactors = F)

# getting the first consequence hit
vep_df$Consequence <- unlist(lapply(vep_df$Consequence, 
                             function(x) as.character(unlist(strsplit(x, split = ',', fixed = T)))[1]))

## MUPIT VEP consequences classification

# loss-of-function 
LOF_CQ = c("stop_gained", 
           "splice_acceptor_variant", 
           "splice_donor_variant", 
           "frameshift_variant", 
           "initiator_codon_variant", 
           "start_lost", 
           "conserved_exon_terminus_variant")

# protein altering variant
MISSENSE_CQ = c("missense_variant", 
                "stop_lost", 
                "inframe_deletion", 
                "inframe_insertion", 
                "coding_sequence_variant",
                "protein_altering_variant")

# silent
SYNONYMOUS_CQ = c("synonymous_variant")

# number of variants per consequence group
l <- lapply(list(LOF_CQ, MISSENSE_CQ, SYNONYMOUS_CQ), 
            function(x, y=vep_df$Consequence)  length(y[y %in% x]))

# generate gene table with consequence stats
vep_df %>%
  group_by(SYMBOL) %>% 
  summarise(n_ptv = sum(Consequence %in% LOF_CQ),
            n_missense = sum(Consequence %in% MISSENSE_CQ),
            n_syn = sum(Consequence %in% SYNONYMOUS_CQ)) %>%
  mutate(n_nonsyn = n_ptv + n_missense) -> df

write.table(df,
            file='data/gene_variant_table_counts_02122019.tsv',
            row.names = F,
            quote = F,
            sep = '\t')



