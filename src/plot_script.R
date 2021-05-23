require(phyloseq)
require(tidyverse)
require(patchwork)

OTU <- readRDS("dada_table.rds")
TAX <- readRDS("dada_tax.rds")

metadata <- read.table("metadata.txt", sep=",", header=FALSE, row.names=1)
colnames(metadata) <- c("Study", "Plot", "Treatment", "Season")

ps <- phyloseq(otu_table(OTU, taxa_are_rows=FALSE),
               tax_table(TAX), sample_data(metadata))

dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

saveRDS(ps, file="dada_physeq.rds")

ps.s <- subset_samples(ps, Season=="s")
ps.a <- subset_samples(ps, Season=="a")
tmp.s <- metaMDS(otu_table(ps.s))
tmp.a <- metaMDS(otu_table(ps.a))
tmp.s <- data.frame(tmp.s$points, Treatment=sample_data(ps.s)$Treatment)
tmp.a <- data.frame(tmp.a$points, Treatment=sample_data(ps.a)$Treatment)
p1 <- ggplot(tmp.s, aes(MDS1, MDS2, fill=Treatment)) + geom_point(pch=22, size=5) + theme_classic() + theme(panel.border=element_rect(color="black", fill="NA", size=1)) + scale_y_continuous(limits=c(-0.35, 0.35)) + scale_x_continuous(limits=c(-0.35, 0.35))
p2 <- ggplot(tmp.a, aes(MDS1, MDS2, fill=Treatment)) + geom_point(pch=22, size=5) + theme_classic() + theme(panel.border=element_rect(color="black", fill="NA", size=1)) + scale_y_continuous(limits=c(-0.35, 0.35)) + scale_x_continuous(limits=c(-0.35, 0.35))
p1 + p2 + plot_layout(guides = 'collect') & theme(legend.position='bottom')
