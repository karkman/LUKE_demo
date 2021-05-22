require(phyloseq)
require(tidyverse)

OTU <- readRDS("dada_table.rds")
TAX <- readRDS("dada_tax.rds")

metadata <- read.table("metadata.txt", sep="\t", header=TRUE)

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


vegdist(otu_table(ps))
p1 <- plot_ordination(subset_samples(ps, Season=="s"), ordinate(subset_samples(ps, Season=="s"), "NMDS", "bray"), color="Treatment") + geom_point(size=5, pch=15) + theme_classic()
p2 <- plot_ordination(subset_samples(ps, Season=="a"), ordinate(subset_samples(ps, Season=="a"), "NMDS", "bray"), color="Treatment") + geom_point(size=5, pch=15) + theme_classic()
p1+p2 +  plot_layout(guides = 'collect') & theme(legend.position='bottom')

#p2 <- plot_ordination(subset_samples(ps, Season=="a"), ordinate(subset_samples(ps, Season=="a"), "NMDS", "bray"), color="Treatment") + geom_point(size=5, pch=21) + theme_classic() + scale_color_manual(values=c("darkgreen", "yellow", "pink", "grey90"))

ps.s <- subset_samples(ps, Season=="s")
ps.a <- subset_samples(ps, Season=="a")
tmp.s <- metaMDS(otu_table(ps.s))
tmp.a <- metaMDS(otu_table(ps.a))
tmp.s <- data.frame(tmp.s$points, Treatment=sample_data(ps.s)$Treatment)
tmp.a <- data.frame(tmp.a$points, Treatment=sample_data(ps.a)$Treatment)
p1 <- ggplot(tmp.s, aes(MDS1, MDS2, fill=Treatment)) + geom_point(pch=22, size=5) + theme_classic() + theme(panel.border=element_rect(color="black", fill="NA", size=1)) + scale_y_continuous(limits=c(-0.35, 0.35)) + scale_x_continuous(limits=c(-0.35, 0.35))
p2 <- ggplot(tmp.a, aes(MDS1, MDS2, fill=Treatment)) + geom_point(pch=22, size=5) + theme_classic() + theme(panel.border=element_rect(color="black", fill="NA", size=1)) + scale_y_continuous(limits=c(-0.35, 0.35)) + scale_x_continuous(limits=c(-0.35, 0.35))
p1 + p2 + plot_layout(guides = 'collect') & theme(legend.position='bottom')
