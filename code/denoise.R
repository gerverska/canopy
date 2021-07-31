# Denoise sequences and identify amplicon sequence variants (ASVs) ####

# Load packages ####
library(Biostrings)
library(phyloseq)
library(ShortRead)
library(dada2)
library(tidyverse)
library(magrittr)
library(here)

# Make path for output ####
out.path <- here('output', 'denoise')
unlink(out.path, recursive = T)
dir.create(out.path)

# Read in fwd and rev reads ####
in.path <- here('output', 'trim')
path.fwd <- in.path %>% list.files(pattern = '.R1.fq.gz', full.names = T) %>% sort 
path.rev <- in.path %>% list.files(pattern = '.R2.fq.gz', full.names = T) %>% sort 

# Make file names / paths for trimmed and filtered files ####
filt.path <- here(out.path, 'filt')
filt.fwd <- path.fwd %>% gsub(in.path, filt.path, .)
filt.rev <- path.rev %>% gsub(in.path, filt.path, .)

# Preview read quality before trimming ####
qual.samps <- sample(1:length(path.fwd), 12)
qual.fwd <- plotQualityProfile(path.fwd[qual.samps]) + ggtitle('Quality profiles fwd')
here(out.path, 'qual.fwd.pdf') %>% ggsave(width=12, height=9)
qual.rev <- plotQualityProfile(path.rev[qual.samps]) + ggtitle('Quality profiles rev')
here(out.path, 'qual.rev.pdf') %>% ggsave(width=12, height=9)
  
# Trim and quality filter ####
trim <- filterAndTrim(path.fwd, filt.fwd,
                      path.rev, filt.rev,
                      maxEE = c(2,2),
                      multithread = TRUE)
  
# Update list of trimmed file paths to exclude samples with no reads passing filters ####
filt.fwd <- list.files(filt.path, pattern = 'R1.fq.gz', full.names = T)
filt.rev <- list.files(filt.path, pattern = 'R2.fq.gz', full.names = T)

# Check quality of trimmed and filtered reads ####
qual.samps <- sample(1:length(filt.fwd), 12)
qual.filt.fwd <- plotQualityProfile(filt.fwd[qual.samps]) + ggtitle('Qual profiles fwd filtered')
here(out.path, 'qual.fwd.filtered.pdf') %>% ggsave(width = 12, height = 9)
qual.filt.rev <- plotQualityProfile(filt.rev[qual.samps]) + ggtitle('Qual profiles rev filtered')
here(out.path, 'qual.rev.filtered.pdf') %>% ggsave(width = 12, height = 9)
  
# Dereplicate ####
derep.fwd <- derepFastq(filt.fwd, verbose = F)
derep.rev <- derepFastq(filt.rev, verbose = F)

# Trim names of derep objects ####
names(derep.fwd) %<>% gsub('.R1.fq.gz', '', .)
names(derep.rev) %<>% gsub('.R2.fq.gz', '', .)
  
# Learn errors. 1e09 bases will take a long time. ####
err.fwd <- learnErrors(filt.fwd, multithread = TRUE, nbases = 1e09, randomize = T) 
err.plot.fwd <- plotErrors(err.fwd, nominalQ = TRUE) + ggtitle(paste('Forward reads error model'))
here(out.path, 'errMod.fwd.pdf') %>% ggsave(width = 5, height = 5)
err.rev <- learnErrors(filt.rev, multithread = TRUE, nbases = 1e09, randomize = T)
err.plot.rev <- plotErrors(err.rev, nominalQ = TRUE) + ggtitle(paste('Reverse reads error model'))
here(out.path, 'errMod.rev.pdf') %>% ggsave(width = 5, height = 5)
 
# Denoise ####
dada.fwd <- dada(derep.fwd, err = err.fwd, multithread = TRUE, pool = 'pseudo')
dada.rev <- dada(derep.rev, err = err.rev, multithread = TRUE, pool = 'pseudo')
  
# Merge reads and make sequence table ####
merged <- mergePairs(dada.fwd, derep.fwd, dada.rev, derep.rev, trimOverhang = T)
seq.tab <- merged %>% makeSequenceTable()

# Make summary report ####
get.n <- function(x) sum(getUniques(x))
trim.summary <- trim %>% data.frame %>% rownames_to_column('sample') 
trim.summary$sample %<>% strsplit(., '.', fixed = T) %>% sapply(., `[`, 1)

track <- cbind(sapply(dada.fwd, get.n),
               sapply(dada.rev, get.n),
               sapply(merged, get.n)) %>%
  data.frame %>% rownames_to_column('sample')

track$sample %<>% strsplit(., '.', fixed = T) %>% sapply(., `[`, 1)
track %<>% left_join(trim.summary, .)
colnames(track) <- c('sample', 'input', 'filtered', 'denoised.fwd', 'denoised.rev', 'merged')
here(out.path, 'dada.summary.rds') %>% saveRDS(track, .)

# Save output ####
here(out.path, 'seq.tab.rds') %>% saveRDS(seq.tab, .)