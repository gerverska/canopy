# Single out sequences belonging to needles that have undergone LULU post-clustering ####
# and PERFect filtering, and attach point cloud data to each phyloseq object.

# To do:
# > determine whether removal of low-depth samples creates different results
# > determine whether PC derived from PCA of height, depth, and closure works better for interpretation [ANSWER: It doesn't really seem like it...]

# Load packages ####
library(phyloseq)
library(tidyverse)
library(magrittr)
library(here)
# Consult "code/functions.R" for descriptions of function applications, inputs, and outputs
source(here('code', 'functions.R')
       )

# Set output directory ####
out.path <- here('data', 'compile')
dir.create(out.path, recursive = T)

# Read in phyloseq list ####
phy.list <- readRDS(here(out.path, 'phy.list.rds'))

# Load crown data for study trees and randomly sampled trees ####
study.crowns <- readRDS(here(out.path, 'study.crown.rds'))
min.zmed <- study.crowns$zmed %>% min()
rand.crowns <- readRDS(here(out.path, 'trad.random.crown.rds')) %>%
  select(-exp.top, -top.diff) %>% filter(zmed >= min.zmed)
rand.crowns$x <- NA
rand.crowns$y <- NA
rand.crowns$z <- NA
all.crowns <- rbind(study.crowns, rand.crowns)

# Examining correlations between crown variables for sampled trees ####
# Added in response to Reviewer #1 comment
pairs(study.crowns[, c('height', 'depth', 'closure')])
lm(closure ~ height, study.crowns) %>% summary()
# Strong correlation between depth and closure suggests removal of depth
lm(closure ~ depth, study.crowns) %>% summary()
lm(depth ~ height, study.crowns) %>% summary()

# Fit loess curve and assign exposure groups ####
# all.loess <- loess.index(all.crowns, 'closure ~ height')
# all.loess$group <- ifelse(all.loess$residuals > 0, 'Closed', 'Open')
# resid.crowns <- all.loess[1:71, ]
# resid.crowns$tree %<>% str_replace('DT_NEIGHBOR', 'DT NEIGHBOR')

k.crowns <- kmean.it(all.crowns, var = 'closure')
k.crowns$group <- ifelse(k.crowns$group == 1, 'Open', 'Closed')
k.crowns$tree %<>% str_replace('DT_NEIGHBOR', 'DT NEIGHBOR')

perf.n <- joiner(phy.list$perf.n)
lulu.n <- joiner(phy.list$lulu.n)

# Process phyloseq objects ####
a1.tree.ht <- perf.n %>% subset_samples(age == 'A1') %>%
  subset_samples(sampleID != 'P03_05_A' & sampleID != 'P02_01_E' & sampleID != 'P02_04_A' & sampleID != 'P02_05_E') %>%
  filt.n.prune() %>% sample_data() %>% data.frame() %>% .$tree.ht

perf.n.clean <- phy.clean(perf.n)
lulu.n.clean <- phy.clean(lulu.n)

saveRDS(perf.n.clean, here(out.path, 'perf.n.clean.rds'))
saveRDS(lulu.n.clean, here(out.path, 'lulu.n.clean.rds'))

# Get session info ####
session.path <- here('output', 'sessions')
dir.create(session.path, recursive = T)
sessionInfo() %>% saveRDS(here(session.path, 'clean.sesh.rds'))
