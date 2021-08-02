# Make phyloseq objects that only include fungal ITS2 sequences ####
# ITS2 sequences were extracted using ITSx version 1.1b1
# Sequence alignments were carried out using VSEARCH v2.8.5_linux_x86_64 (https://github.com/torognes/vsearch)

# Load packages
library(ggthemes)
library(PERFect)
library(lulu)
library(phyloseq)
library(Biostrings)
library(dada2)
library(tidyverse)
library(magrittr)
library(here)
# Consult "code/functions.R" for descriptions of function application, inputs, and outputs
source(here('code', 'functions.R'))
 
# Load inputs ####
meta <- readRDS(here('data', 'compile', 'canopy.meta.rds'))
seq.tab <- readRDS(here('output', 'denoise', 'seq.tab.rds'))

# Create output directories ####
out.path <- here('output', 'compile')
unlink(out.path, recursive = T)
dir.create(out.path)
scratch.path <- here('output', 'scratch')
unlink(scratch.path, recursive = T)
dir.create(scratch.path)

# Extract the ITS2 region from all sequences, keeping only full ITS2 sequences ####
otu.names <- paste0('otu.', 1:ncol(seq.tab))
seqs <- getSequences(seq.tab) %>% dada2::rc() %>% DNAStringSet()
names(seqs) <- otu.names
colnames(seq.tab) <- otu.names
writeXStringSet(seqs, file = here(scratch.path, 'tmp.fasta'), width = 600)

itsx.flags <- paste('-i output/scratch/tmp.fasta',
                    '-t "fungi,tracheophyta"',
                    '--preserve T',
                    '--complement F',
                    '--summary T',
                    '--cpu 10',
                    '--multi_thread T',
                    '-o output/scratch/ITSx',
                    '--only_full T',
                    '-E 1e-2')
system2('ITSx', args = itsx.flags)

# Read in ITS2 sequences from ITSx output and filter out very short reads
seqs <- readDNAStringSet(here(scratch.path, 'ITSx.ITS2.fasta')) %>% .[.@ranges@width > 50]

# Remove OTUs removed by ITSx
seq.tab %<>% .[, names(seqs)]

# Collapse identical OTUs
colnames(seq.tab) <- seqs %>% as.character() %>% unname()
seq.tab %<>% collapseNoMismatch()

# Remove chimeras
seq.tab %<>% removeBimeraDenovo(method = 'consensus',
                                multithread = TRUE,
                                verbose = TRUE)

# Reassign OTU names
otu.names <- paste0('otu.', 1:ncol(seq.tab))
seqs <- getSequences(seq.tab) %>% DNAStringSet()
names(seqs) <- otu.names
colnames(seq.tab) <- otu.names

# Remove non-fungal sequences by aligning against the UNITE all-eukaryote release (singletons included) ####
# 5 fungal taxa were added retroactively to the all-euk and fungal releases
writeXStringSet(seqs, file = here(scratch.path, 'tmp.fasta'), width = 600)
vsearch.flags <- paste('--usearch_global output/scratch/tmp.fasta',
                       '--db data/tax/sh_general_release_dynamic_s_all_04.02.2020_et_al.fasta',
                       '--id 0.50',
                       '--userout output/scratch/unite.matches.txt',
                       '--userfields query+target+id',
                       '--notmatched output/scratch/no.match.fasta',
                       '--maxhits 5',
                       '--maxaccepts 500 --maxrejects 0'
)
system2('vsearch', args = vsearch.flags)

# Read in search results and calculate the proportion of top matches that are fungal
unite.matches <- read.table(here(scratch.path, 'unite.matches.txt'), as.is = T) %>%
  group_by(V1) %>%
  summarise(is.fungi = ifelse(sum(!grepl('k__Fungi', V2)) == 0, 1,
                            sum(grepl('k__Fungi', V2))/sum(!grepl('k__Fungi', V2))))

# Keep ASVs that have greater than 0.5 fungal matches to the UNITE euk release
fungal <- unite.matches$V1[unite.matches$is.fungi > 0.5] %>%
  c(readDNAStringSet(here(scratch.path, 'no.match.fasta')) %>% names())

seqs %<>% .[.@ranges@NAMES %in% fungal]
seq.tab %<>% .[, fungal]

# LULU post-clustering ASV curation ####
# Make an OTU table with samples as columns
otu.tab <- seq.tab %>% t() %>% as.data.frame()

# Make a matchlist by aligning all remaining sequences against each other
writeXStringSet(seqs, file = here(scratch.path, 'tmp.fasta'), width = 600)

vsearch.flags <- paste('--usearch_global output/scratch/tmp.fasta',
                       '--db output/scratch/tmp.fasta',
                       '--self',
                       '--id 0.95',
                       '--strand plus',
                       '--iddef 1',
                       '--userout output/scratch/matchlist.txt',
                       '--userfields query+target+id',
                       '--maxaccepts 0',
                       '--query_cov 0.9',
                       '--maxhits 10')
system2('vsearch', args = vsearch.flags)

matchlist <- read.table(here(scratch.path, 'matchlist.txt'),
                        header = FALSE,
                        as.is = TRUE,
                        stringsAsFactors = FALSE)

# Make a curation object and pull out the merged OTU table
# minimum_match = 97, ensuring that we are at least collapsing ASVs that might normally be called as the same
# OTU, minimizing post-clustering false positive (child == parent) errors. However, sequence identities of 
# intragenomic variants of the ITS region are known to vary by as much as 87%. This means that not all true 
# intragenomic ASVs will collapse into a single OTU--i.e. we take on board more false negatives
# (child != parent).
curation <- lulu(otu.tab, matchlist, minimum_match = 97)

# move log, save curation object, reorient output OTU table, and keep parent sequences
system2('mv', args = c('lulu.log*', 'output/compile'))
here(out.path, 'lulu.curation.rds') %>% saveRDS(curation, .)
lulu.otu.tab <- curation$curated_table %>% t() %>% as.data.frame()
seqs %<>% .[.@ranges@NAMES %in% curation$curated_otus]

# Taxonomy prediction ####
tax <- assignTaxonomy(seqs,
                      here('data', 'tax', 'sh_general_release_dynamic_04.02.2020_et_al.fasta'),
                      multithread = T)
rownames(tax) <- names(seqs)

# Make, clean up, and split a phyloseq object ####
phy <- phyloseq(otu_table(lulu.otu.tab, taxa_are_rows = F), 
                sample_data(meta),
                tax_table(tax),
                refseq(seqs))

# Fix taxonomy table text
phy %<>% parse.tax()

# Remove twig extraction test samples and remove empty taxa without reads
phy %<>% subset_samples(control != 'Test')
phy %<>% filter_taxa(function(x) { sum(x) > 0 }, TRUE)

# Remove samples without any reads
keep.reads <- sample_sums(phy) > 0
phy %<>% prune_samples(keep.reads, .)

# Set aside control samples and identify potential contaminants
controls <- subset_samples(phy, control != 'No') %>%
  filter_taxa(function(x) { sum(x) > 0 }, TRUE)
control.sums <- controls@otu_table %>% data.frame() %>% colSums()

# Remove control samples and separate needle and twig samples
phy %<>% subset_samples(control == 'No')
needle <- subset_samples(phy, organ == 'Needle') %>%
  filter_taxa(function(x) { sum(x) > 0 }, TRUE)

# remove OTUs with more than 1% of study reads in controls
needle %<>% remove.contam(control.sums, prop = 0.01)

twig <- subset_samples(phy, organ == 'Twig') %>%
  subset_samples(tree != 'DISCOVERY' & tree != 'DT_NEIGHBOR' & tree != 'PC17') %>%
  filter_taxa(function(x) { sum(x) > 0 }, TRUE)
twig %<>% remove.contam(control.sums, prop = 0.01)

# After separating needle and twig samples, retain the original phyloseq object in order to analyze both sample sets together.
# This requires that we drop needle data for three trees missing twig data.
phy %<>% subset_samples(organ == 'Needle' | organ == 'Twig') %>%
  subset_samples(tree != 'DISCOVERY' & tree != 'DT_NEIGHBOR' & tree != 'PC17') %>%
  filter_taxa(function(x) { sum(x) > 0 }, TRUE)
phy %<>% remove.contam(control.sums, prop = 0.01)

# PERFect OTU filtering ####
perfect.needle <- perfect.filt(needle, 1000, 'needle')
perfect.twig <- perfect.filt(twig, 1000, 'twig')
perfect.both <- perfect.filt(phy, 1000, 'both')

# Assemble output phyloseq list ####
phy.list <- list(lulu.n = needle, 
                 perf.n = perfect.needle,
                 lulu.t = twig,
                 perf.t = perfect.twig,
                 lulu.b = phy,
                 perf.b = perfect.both)

# remove temporary files
unlink(scratch.path, recursive = T)
here(out.path, 'phy.list.rds') %>% saveRDS(phy.list, .)