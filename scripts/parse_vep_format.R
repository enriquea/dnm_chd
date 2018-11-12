# Simple script to format SNP data for VEP annotation

library(readr)
library(plyr)


dnm_janin_sifrim <- read_delim("analysis/dnm/dnm_janin_sifrim_formatted.tsv",
                               "\t", 
                               escape_double = FALSE, 
                               trim_ws = TRUE)


df <- as.data.frame(dnm_janin_sifrim)

# The default (VEP) format is a simple whitespace-separated format (columns may be separated by space or tab characters), 
# containing five required columns plus an optional identifier column:
#   
#   1. chromosome - just the name or number, with no 'chr' prefix
#   2. start
#   3. end (start + length(ref) - 1) !
#   4. allele - pair of alleles separated by a '/', with the reference allele first
#   5. strand - defined as + (forward) or - (reverse).
#   6. identifier - this identifier will be used in VEP's output. 
#      If not provided, VEP will construct an identifier from the given coordinates and alleles.

df_vep <- data.frame(chromosome = df$chrom,
                     start = df$start_pos,
                     end = df$start_pos + nchar(df$ref_allele) - 1,  # -1 is important
                     allele = paste(df$ref_allele, df$alt_allele, sep = '/'))

write.table(df_vep, 
            file = 'analysis/dnm/dnm_janin_sifrim_formatted_for_vep.txt',
            sep = ' ',
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)


### Format DNM file to run MUPIT tool

# expected columns 
#
# "person_id": ["temp"],
# "chrom": ["16"],
# "start_pos": [89348744],
# "end_pos": [89348744],
# "ref_allele": ["A"],
# "alt_allele": ["G"],
# "hgnc": ["ANKRD11"],
# "consequence": ["missense_variant"],
# "study_code": [None],
# "publication_doi": [None],
# "study_phenotype": [None],
# "type": ["snv"] 

# parsing output vep
vep_annotated <- read.delim("analysis/dnm/sifrim_janin_vep_reannotated_2018.txt", stringsAsFactors=FALSE)

# retriving variant information
variant_key <- lapply(vep_annotated$Uploaded_variation, function(x) unlist(strsplit(x, split = '_')))

ref_alt <- lapply(variant_key, function(x) unlist(strsplit(x[3], split = '/')))

refs <- vapply(ref_alt, function(x) x[[1]], character(1))

alts <- vapply(ref_alt, function(x) x[[2]], character(1))

# build df

dnm <- data.frame(person_id = 'None',
                  chrom = unlist(lapply(variant_key, function(x) x[1])),
                  start_pos = unlist(lapply(variant_key, function(x) x[2])),
                  end_pos = unlist(lapply(variant_key, function(x) x[2])),
                  ref_allele = unlist(lapply(ref_alt, function(x) x[1])),
                  alt_allele = unlist(lapply(ref_alt, function(x) x[2])),
                  hgnc = vep_annotated$SYMBOL,
                  consequence = vep_annotated$Consequence,
                  study_code = 'None',
                  publication_doi = 'None',
                  study_phenotype = 'None',
                  type = ifelse(nchar(refs) != nchar(alts), "indel", "snv"))

# parse cosequence column
dnm$consequence <- unlist(lapply(as.character(dnm$consequence), function(x) unlist(strsplit(x, split = ',', fixed = T))[1]))


write.table(dnm, file = "analysis/dnm/dnm_janin_sifrim_formatted_30042018.tsv", 
            row.names = FALSE, 
            sep = '\t', 
            quote = FALSE)
