# Obtain metrics derived from ALS point clouds for each sampled tree at each sampled height ####

library(gtools)
library(phyloseq)
library(future)
library(lidR)
library(tidyverse)
library(magrittr)
library(here)
# Consult "code/functions.R" for descriptions of function applications, inputs, and outputs
source(here('code', 'functions.R'))

# Set number of cores, minimum tree heights, and number of trees to sample
plan(multisession, workers = 8L)

# Make an output directory
out.path <- here('data', 'compile')
dir.create(out.path, recursive = T)

# Read in tree coordinates and sampled heights
hts <- readRDS(here(out.path, 'study.trees.rds'))
trees <- hts %>% distinct(tree, x, y)

# Clips out 10 m sections from each of the .laz files
cat <- readLAScatalog(here('data', 'lidar', 'raw'),
                      filter = '-keep_user_data 0 -keep_intensity_above 0',
                      chunk_size = 0,
                      chunk_buffer = 25)
opt_laz_compression(cat) <- T
opt_output_files(cat) <- here(out.path, 'clips.10m', '{ID}_10m')

clip_circle(cat, xcenter = trees[, 2], ycenter = trees[, 3], radius = 10)

clips.10m <- list.files(here(out.path, 'clips.10m'), full.names = T) %>%
  sort() %>% mixedsort()

crown <- study.crown(clips.10m, trees, hts, 8)

saveRDS(crown, here(out.path, 'study.crown.rds'))

# Get session info ####
session.path <- here('output', 'sessions')
dir.create(session.path, recursive = T)
sessionInfo() %>% saveRDS(here(session.path, 'study.sesh.rds'))
