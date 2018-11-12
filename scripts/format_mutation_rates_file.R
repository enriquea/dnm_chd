#### Formatting mutation rate file ####

library(rtracklayer)
library(GenomicRanges)
library(geneSynonym)
library(dplyr)
library(plyr)

# reading mutation rate file
gene_mutation_rates <- read_excel('resourses/gene_mutation_rates_ng3050s2.xls',
                                  sheet = "mutation_probabilities")

# preparing reference GENCODE v19
genecode19 <- rtracklayer::readGFF(filepath = 'resourses/gencode.v19.annotation.gtf.gz')
genecode19 <- GenomicRanges::as.data.frame(genecode19)
genes_references <- unique(genecode19$gene_name)

# getting missing genes from mutation-rate list
query_genes <- gene_mutation_rates$gene
query_missing <- query_genes[!query_genes %in% genes_references]

# getting all possible aliases for missing genes
query_missing_aliases <- humanSyno(query_missing)
query_missing_aliases = lapply(query_missing_aliases, function(x) unique(unlist(x)))

# mapping aliases to GENCODE v19 (get first hit!!!)
map_to_gencode19 <- lapply(query_missing_aliases, function(x) x[x %in% genes_references][1])

# replace NA (missing aliases) with original names
original_names <- names(map_to_gencode19)
formatted_list <- list()
for(nm in original_names){
    formatted_list[nm] <- ifelse(is.na(map_to_gencode19[nm]), nm, map_to_gencode19[nm])
}

## getting final reference file

# replace symbol from mutation rate with GENCODE v19 aliases
temp <- gene_mutation_rates
temp$gene_alias <- unlist(replace(temp$gene, temp$gene %in% names(formatted_list), formatted_list))

# remove duplicated rows (gene_aliases)
temp <- temp[!duplicated(temp$gene_alias),]
genecode19 <- unique(data.frame(chr=genecode19$seqid, gene=genecode19$gene_name))
merged <- merge.data.frame(temp, genecode19, by.x = 'gene_alias', by.y = 'gene')

rates <- merged[,c('all', 'syn', 'mis', 'non', 'splice_site', "frameshift")]
rates <- as.data.frame(sapply(rates, function(x) 10**as.numeric(x)))

final_rates_ref <- cbind(hgnc=merged$gene_alias, chrom=gsub(merged$chr, replacement = '', pattern = 'chr', fixed = T), rates)

# export table
write.table(final_rates_ref, 'mutation_rates_formatted_31072018.tsv', sep = '\t', quote = F, row.names = F)
