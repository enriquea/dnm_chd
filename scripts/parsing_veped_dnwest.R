library(readr)
library(dplyr)

## Simple script for parsing VEPed tsv file (from web ensembl VEP tool)

# @author enriquea
# @date 122019

# read VEP annotated file
vep_annotated <- read.delim("dnm_jin_sifrim_veped_02122019.txt", 
                            stringsAsFactors=FALSE,
                            na.strings = '-')

# getting variant information
variant_key <- lapply(vep_annotated$Uploaded_variation, 
                      function(x) unlist(strsplit(x, split = '_')))

ref_alt <- lapply(variant_key, 
                  function(x) unlist(strsplit(x[3], split = '/')))

refs <- vapply(ref_alt, function(x) x[[1]], character(1))

alts <- vapply(ref_alt, function(x) x[[2]], character(1))

# building input for denovoWEST tool
vep_annotated %>% 
  mutate(alt = Allele,
         altprop_child = NA,
         chrom = unlist(lapply(variant_key, function(x) x[1])),
         consequence = Consequence, 
         cq = Consequence,
         hgnc_id = HGNC_ID,
         id = 'None',
         maf = ifelse(is.na(gnomAD_AF), 0, gnomAD_AF),
         pos = as.numeric(unlist(lapply(variant_key, function(x) x[2]))),
         prob = NA,
         raw = CADD_RAW,
         ref = unlist(lapply(ref_alt, function(x) x[1])),
         score = CADD_PHRED,
         score_mpc = MPC,
         study = 'None',
         symbol = SYMBOL,
         constrained = 'True') %>%
    select(alt, altprop_child, chrom, consequence, cq, hgnc_id, id, maf, pos, 
           prob, raw, ref, score, score_mpc, study, symbol, constrained) -> dnm
  

# keep only the most severe variant consequnece
dnm$consequence <- unlist(lapply(as.character(dnm$consequence), 
                                 function(x) unlist(strsplit(x, split = ',', fixed = T))[1]))
dnm$cq <- dnm$consequence

# remove entries with missing position information
dnm <- dnm[!is.na(dnm$pos),]

## Annotate missense constraint information (needed for denovoWEST)
## Criteria:
## 1. missense and stop gained variants are annotated as constrained if fell
## in a constrained gene region (o/e < 0.6).
## 2. others variant consequences are annotated as contrained if missense contraint
## evidence was found for the gene (missense z-score > 3).

miss_cons_genes <- unique(union(subset(miss_table, mis_z > 3)['gene'], 
                                miss_cons['gene']))

annotate_constraint <- function(gene, pos, consequence=NULL){
         cons <- miss_cons
         if (gene %in% cons$gene & consequence %in% c('stop_gained', 'missense_variant')){
                  cons %>% filter(gene == gene,
                                  pos >= genomic_start,
                                  pos <= genomic_end,
                                  obs_exp <= 0.6) -> cons
            return(ifelse(nrow(cons) > 0, 'True', 'False'))
         } 
         else if(gene %in% miss_cons_genes$gene & !consequence %in% c('stop_gained', 'missense_variant')) 
         {
             return('True')
         } else {
           return('False')
         }
}

dnm$constrained <- apply(dnm, 1, function(x) annotate_constraint(gene = x['symbol'], 
                                                                 pos = x['pos'],
                                                                 consequence = x['cq']))

# export table
write.table(dnm, file = "~/tools/DeNovoWEST/dnm_jin_sifrim_formatted_for_denovoWEST.tsv", 
            row.names = FALSE, 
            sep = '\t', 
            quote = FALSE)


#### export TSV DNM file for denovonear ####

# formatted file to run denovonear tool
dnm %>% mutate(gene_name = as.character(symbol),
               chr = as.character(chrom),
               pos = pos,
               consequence = as.character(consequence),
               snp_or_indel = ifelse(nchar(as.character(ref)) == 1 & 
                                       nchar(as.character(alt)) == 1, 
                                     "DENOVO-SNP", "DENOVO-INDEL")) %>%
        dplyr::select(c(gene_name, chr, pos, consequence, snp_or_indel)) -> dnm_denovonear

write.table(dnm_denovonear, file = "~/tools/denovonear/dnm_jin_sifrim_formatted_for_denovonear.tsv", 
            row.names = FALSE, 
            sep = '\t', 
            quote = FALSE)
