# Simple script to format SNP data for VEP annotation

library(readr)
library(plyr)


dnm_jin_sifrim <- read_delim("analysis/dnm/dnm_jin_sifrim_formatted.tsv",
                               "\t", 
                               escape_double = FALSE, 
                               trim_ws = TRUE)


df <- as.data.frame(dnm_jin_sifrim)

# The default (VEP) format is a simple whitespace-separated format (columns may be separated by space or tab characters), 
# containing five required columns plus an optional identifier column:
#   
#   1. chromosome - just the name or number, without 'chr' prefix
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
            file = 'data/dnm_jin_sifrim_formatted_for_vep.txt',
            sep = ' ',
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)


### DNM input file format to run MUPIT tool

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
vep_annotated <- read.delim("data/dnm_jin_sifrim_veped_02122019.txt", 
                            stringsAsFactors=FALSE,
                            na.strings = '-')

# retriving variant information
variant_key <- lapply(vep_annotated$Uploaded_variation, 
                      function(x) unlist(strsplit(x, split = '_')))

ref_alt <- lapply(variant_key, 
                  function(x) unlist(strsplit(x[3], split = '/')))

refs <- vapply(ref_alt, function(x) x[[1]], character(1))

# alts <- vapply(ref_alt, function(x) x[[2]], character(1))

# build dnm df
dnm <- data.frame(person_id = 'None',
                  chrom = unlist(lapply(variant_key, function(x) x[1])),
                  start_pos = unlist(lapply(variant_key, function(x) x[2])),
                  end_pos = unlist(lapply(variant_key, function(x) x[2])),
                  ref_allele = unlist(lapply(ref_alt, function(x) x[1])),
                  alt_allele = vep_annotated$Allele,
                  hgnc = vep_annotated$SYMBOL,
                  consequence = vep_annotated$Consequence,
                  study_code = 'None',
                  publication_doi = 'None',
                  study_phenotype = 'None')

# annotate type of variant (snv or indel)
dnm$type <- ifelse(nchar(as.character(dnm$alt_allele)) == 1 
                   & nchar(as.character(dnm$ref_allele)) == 1, 
                   "snv", "indel")

# parse consequence column (get most severe consequence)
dnm$consequence <- unlist(lapply(as.character(dnm$consequence), 
                                 function(x) unlist(strsplit(x, split = ',', fixed = T))[1]))

# remove variant with missing symbol/position
dnm %>% filter(!is.na(start_pos), !is.na(hgnc)) -> dnm

# export table
write.table(dnm, file = "analysis/dnm_jin_sifrim_mupit_input_04122019.tsv", 
            row.names = FALSE, 
            sep = '\t', 
            quote = FALSE)
