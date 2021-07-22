# Identify trees from across the HJA Forest ####
# and obtain point cloud-derived metrics from randomly sampled trees at each randomly sampled height

library(gtools)
library(future)
library(lidR)
library(tidyverse)
library(magrittr)
library(here)
# Consult "code/functions.R" for descriptions of function applications, inputs, and outputs
source(here('code', 'functions.R'))

# Set number of cores, minimum tree heights, and number of trees to sample ####
# Change the number of workers as needed to avoid crashes/speed computation.
plan(multisession, workers = 6L)

# Make an output directory ####
out.path <- here('data', 'compile')
dir.create(out.path, recursive = T)

# Load the catalog of raw .laz files, dropping data with GPS time errors and intensity values < 0 ####
cat <- readLAScatalog(here('data', 'lidar', 'raw'),
                      filter = '-keep_user_data 0 -keep_intensity_above 0',
                      chunk_size = 0,
                      chunk_buffer = 25)

# Set random seed, and randomly sample point cloud clips with radius 50 m ####
# Smaller radii can likely be used, but large radii are more likely to "trap" complete tree crowns.
set.seed(666)
centers <- cat@data %>% rownames_to_column(var = 'file') %>% group_by(file) %>%
  summarize(random.plot(Min.X, Max.X, Min.Y, Max.Y)) %>% select(x.center, y.center)

opt_laz_compression(cat) <- T
opt_output_files(cat) <- here(out.path, 'clips.50m', '{ID}_50m')

clip_circle(cat, xcenter = centers$x.center, ycenter = centers$y.center, radius = 50)

# Normalize the elevation of vegetation points by the elevation of ground points to facilitate tree detection ####
clips.50m <- readLAScatalog(here(out.path, 'clips.50m'),
                            filter = '-keep_user_data 0 -keep_intensity_above 0',
                            chunk_size = 0,
                            chunk_buffer = 25)

opt_laz_compression(clips.50m) <- T
opt_output_files(clips.50m) <- here(out.path, 'normalized', '{ID}_50m_norm')
opt_stop_early(clips.50m) <- F

normalize_height(clips.50m, knnidw())

norm <- list.files(here(out.path, 'normalized'), full.names = T) %>%
  sort() %>% mixedsort()

# Identify trees in each 50 m, normalized point cloud clip ####
trees <- tree.samp(norm, 4)

# Saving and loading the tree data frame here allows you to skip the previous compute-intensive steps!
# Helpful for tweaking parameters!
saveRDS(trees, here(out.path, 'random.trees.rds'))
trees <- readRDS(here(out.path, 'random.trees.rds'))

# Create 10m-radius clips around each randomly selected tree ####
opt_output_files(cat) <- here(out.path, 'clips.10m', '{ID}_10m')
# This step takes a very long time. Set cores appropriately.
clip_circle(cat, xcenter = trees[, 1], ycenter = trees[, 2], radius = 10)

clips.10m <- list.files(here(out.path, 'clips.10m'), full.names = T) %>%
  sort() %>% mixedsort()

# Obtain crown variables at random heights for each randomly selected tree ####
crown <- random.crown(clips.10m, trees, 6)
# As indicated in random.crown(), I think I should just use the exact point coordinates
# identified either by the tree-finding algorithm or manually. This might needlessly drop
# trees from the analysis at the below filter step.
crown %<>% filter(top.diff < 0.5 & top.diff > -0.5)
crown$tree %<>% as.factor()

saveRDS(crown, file.path(out.path, 'random.crown.rds'))

# Get session info ####
session.path <- here('output', 'sessions')
dir.create(session.path, recursive = T)
sessionInfo() %>% saveRDS(here(session.path, 'random.sesh.rds'))
