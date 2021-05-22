require(phyloseq)
require(tidyverse)

OTU <- readRDS("dada_table.rds")
TAX <- readRDS("dada_tax.rds")

ps <- phyloseq(otu_table(OTU, taxa_are_rows=FALSE),
               tax_table(TAX))

dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

saveRDS(ps, file="dada_physeq.rds")

ps %>% sample_sums %>% 
