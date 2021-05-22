#!/usr/bin/env Rscript

require(dada2)
#require(tidyverse)

path="./trimmed/"
#arg = commandArgs(trailingOnly=TRUE)

fnFs <- sort(list.files(path, pattern="_16S.1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_16S.2.fastq.gz", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

plotQualityProfile(fnFs[1:12])
plotQualityProfile(fnRs[1:12])

filtFs <- file.path(path, "filtered", paste0(sample.names, "_16S.1.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_16S.2.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(280,260),
              maxN=0, maxEE=c(1,1), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE)

errF <- learnErrors(filtFs, multithread=TRUE, nbases=1e+07)
errR <- learnErrors(filtRs, multithread=TRUE, nbases=1e+07)

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

dadaFs <- dada(filtFs, err=errF, multithread=TRUE, pool=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE, pool=TRUE)

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

seqtab <- makeSequenceTable(mergers)
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
taxa <- assignTaxonomy(seqtab.nochim, "/users/antkark/scratch/databases/SILVA138/training_set.138_SSURef_NR99.fa.gz", multithread=TRUE)
taxa <- addSpecies(taxa, "/users/antkark/scratch/databases/SILVA138/species_assignment.138_SSURef_NR99.fa.gz")

saveRDS(seqtab.nochim, file="dada_table.rds")
saveRDS(taxa, file="dada_tax.rds")
