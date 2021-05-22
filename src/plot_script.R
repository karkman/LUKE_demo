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

tibble(sample=sample_names(ps), sum=sample_sums(ps)) %>% \
	ggplot(aes(x=sample, y=sum)) + 
	geom_bar(stat="identity") + theme_classic() + 
	coord_flip() + scale_y_continuous(expand=c(0,0)) + 
	labs(y="Library size", x="")
