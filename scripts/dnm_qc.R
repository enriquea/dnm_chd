library(dplyr)
library(reshape2)

# Evalute the DNM rates between cohort

# read sample info (Study and variants)
sample_info <- read.delim('raw_data/DNM_Sifrim_Jin_2017.txt',
                          stringsAsFactors = F)

# parse/format df
sample_info %>% 
  mutate(vep_variant_key=paste(CHR,POS,REF,ALT,sep = ':')) %>%
  select(vep_variant_key, ID, Study) -> sample_info

# read/parse VEPed variant table
vep_annotated <- read.delim("data/dnm_jin_sifrim_formatted_for_denovoWEST_122019.tsv", 
                            stringsAsFactors=FALSE,
                            na.strings = '-')


vep_annotated %>% 
  mutate(vep_variant_key=paste(chrom,pos,ref,alt, sep = ":")) -> vep_annotated

# merge data
df_merged <- merge.data.frame(sample_info, vep_annotated, by = 'vep_variant_key')

# annotate consequence group
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

# annoate csq group and keep only rare variant (MAF 1% gnomad)
df_merged %>%
  mutate(csq_group=case_when(cq %in% LOF_CQ ~ 'ptv',
                             cq %in% MISSENSE_CQ ~ 'miss',
                             cq %in% SYNONYMOUS_CQ ~ 'syn')) %>%
  filter(maf < 0.01) %>%
  select(Study, ID, csq_group) -> df_merged

# generate count table summary per study
df_merged %>%
  group_by(Study) %>% 
  dplyr::summarise(PTV=sum(csq_group=="ptv", na.rm = T), 
                   Missense=sum(csq_group=="miss", na.rm = T), 
                   Synonymous=sum(csq_group=='syn', na.rm = T)) -> df_count

df_count <- dcast(melt(df_count, id.vars = c("Study")), 
                       variable ~ ..., value.var = "value")


# add total number of trios
df_count$n_trios_jin <- length(unique(subset(df_merged, Study=='Jain_2017')$ID))
df_count$n_trios_sifrim <- length(unique(subset(df_merged, Study=='Sifrim_2016')$ID))

# Compare rates (Poisson test) per variant class
p_vector <- list()
rates <- list()
total_trios_jin <- 1784
total_trios_sifrim <- 704

for(i in 1:nrow(df_count)){
   v1 <- c(df_count[i, 'Jain_2017'], df_count[i, 'Sifrim_2016'])
   v2 <- c(total_trios_jin, total_trios_sifrim)
   poi <- poisson.test(v1, v2)
   p_vector[[i]] <- as.numeric(poi$p.value)
   rates[[i]] <- as.numeric(poi$estimate)
}

df_count$rates <- rates
df_count$p_poisson <- p_vector

# export results
openxlsx::write.xlsx(df_count, file = 'analysis/qc_dnm_08012020.xlsx')


  


