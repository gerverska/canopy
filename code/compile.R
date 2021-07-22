# Make phyloseq objects that only include fungal ITS2 sequences ####
# ITS2 sequences were extracted using ITSx version 1.1b1
# Sequence alignments were carried out using vsearch v2.8.5_linux_x86_64 (https://github.com/torognes/vsearch)

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

# Functions ####
# Remove otus with >1% of experimental reads found in the pcr and extraction controls
remove.contam <- function(phy, controls, prop) {
  
  # Get the current named list of otu read totals
  otu.tab <- phy@otu_table %>% data.frame()
  otu.names <- colnames(otu.tab)
  
  # Identify otus found in the negative controls
  contam.sums <- otu.tab[, colnames(otu.tab) %in% names(controls)] %>% colSums()
  
  # Calculate a ratio of the number of reads (in negative controls):(in non-control samples) for each otu
  controls <- controls[names(controls) %in% names(contam.sums)]
  
  contam.prop <- controls / contam.sums
  
  # Get the names of the otus with ratio > prop and remove them from the input phyloseq object
  contam.names <- contam.prop[contam.prop > prop] %>% names()
  
  keep <- otu.names[!(otu.names %in% contam.names)]
  
  phy %<>% prune_taxa(keep, .)
  
  keep <- sample_sums(phy) > 0
  prune_samples(keep, phy)
  
}

# Remove noise otus which contribute least to the total covariance of a dataset
perfect.filt <- function(phy, k, file.prefix) {
  
  # For reproducible filtering
  set.seed(666)
  
  otu.tab <- phy@otu_table %>% data.frame()
  
  # Even though we will ultimately use permutation filtering, the authors mention that performance is improved with Order = 'pvals'.
  # This requires that we go ahead and obtain simulataneous PERFect output first. 
  sim <- PERFect_sim(X = otu.tab)
  # Of note here: by default, an alpha of 0.1 is used. Taxa with p-values greater than 0.1 are filtered out.
  # For this reason, filtering is somewhat conservative, which would ideally allow us to retain real and somewhat rare taxa.
  perm <- PERFect_perm(X = otu.tab, Order = 'pvals', pvals_sim = sim, algorithm = 'full', k = k)
  tab.out <- perm$filtX
  
  # Save PERFect output
  here(out.path, paste0(file.prefix, '.', 'perfect.perm.rds')) %>% saveRDS(perm, .)
  
  # Make p-value plots and save to output
  perm %<>% pvals_Plots(otu.tab)
  perm <- perm$plot + ggtitle(paste0('Permutation filtering', '-', file.prefix)) + scale_color_colorblind()
  here(out.path, paste0(file.prefix, '.', 'perfect.perm.pvals.pdf')) %>% ggsave(., perm, dpi = 300)
  
  # Update the phyloseq object
  keep <- tab.out %>% colnames()
  phy %<>% prune_taxa(keep, .)
  otu_table(phy) <- otu_table(tab.out, taxa_are_rows = F)
  keep <- sample_sums(phy) > 0
  prune_samples(keep, phy)
  
}

# Convert ugly UNITE taxonomy to figure-ready taxonomy. Heads up--it's nasty!
parse.tax <- function(phy) {
  
  phy@tax_table@.Data %<>% parse_taxonomy_greengenes()
  
  tax.tab <- phy@tax_table %>% data.frame()
  tax.tab$Genus_species <- ifelse(!is.na(tax.tab$Species), paste0(tax.tab$Genus, ' ', tax.tab$Species),
                                   ifelse(!is.na(tax.tab$Genus), paste0(tax.tab$Genus, ' sp.'),
                                          ifelse(!is.na(tax.tab$Family), paste0(tax.tab$Family, ' sp.'),
                                                 ifelse(!is.na(tax.tab$Order), paste0(tax.tab$Order, ' sp.'),
                                              ifelse(!is.na(tax.tab$Class), paste0(tax.tab$Class, ' sp.'),
ifelse(!is.na(tax.tab$Phylum), paste0(tax.tab$Phylum, ' sp.'),                                                                      ifelse(!is.na(tax.tab$Kingdom), paste0(tax.tab$Kingdom, ' sp.'), NA)))))))
  phy@tax_table <- tax.tab %>% as.matrix() %>% tax_table()
  
  phy
  
}
 
# Load inputs ####
meta <- readRDS(here('data', 'compile', 'full.hja.meta.rds'))
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