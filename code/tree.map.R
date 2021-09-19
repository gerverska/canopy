# Create Fig. 1 from the manuscript, which:
# (a) shows a map of HJA showing the approximate location of each tree,
# (b) the relationship between crown closure and height in crown for climbed and un-climbed trees,
# (c) how each crown variable is calculated.

# Load packages ####
library(patchwork)
library(cowplot)
library(ggthemes)
library(ggnewscale)
library(sf)
library(terra)
library(tidyverse)
library(magrittr)
library(here)
# Consult "code/functions.R" for descriptions of function applications, inputs, and outputs
source(here('code', 'functions.R'))

# Read in phyloseq list
study.crowns <- readRDS(here('data', 'compile', 'study.crown.rds'))
study.crowns$tree %<>% str_replace('DT_NEIGHBOR', 'DT NEIGHBOR')
tree.tab <- study.crowns %>% dplyr::select(tree, x, y) %>% unique()
rand.crowns <- readRDS(here('data', 'compile', 'trad.random.crown.rds'))

# Read in crown variable schematic
crown.vars <- png::readPNG(here('data', 'media', 'crown.vars.png'), native = T)

# Set out.path
figure.out <- here('output', 'figures')
dir.create(figure.out, recursive = T)

# Height/closure profiles for each tree ####
study.crowns %<>% dplyr::select(-x, -y, -z)
min.zmed <- study.crowns$zmed %>% min()
rand.crowns %<>% dplyr::select(-exp.top, -top.diff) %>%
  dplyr::filter(zmed >= min.zmed)
# all.crowns <- rbind(study.crowns, rand.crowns)
# 
# all.loess <- loess.index(all.crowns)
# study.loess <- all.loess[1:71, ]

tree.plot <- ggplot(mapping = aes(x = height, y = closure)) +
  geom_line(data = rand.crowns, aes(group = tree), color = 'gray80', alpha = 0.5, show.legend = F) +
  # geom_smooth(data = all.crowns, aes(x = height, y = closure), color = 'red', formula = y ~ x, method = 'loess', se = F) +
  geom_point(data = study.crowns, mapping = aes(fill = tree), size = 2.5, shape = 21) +
  scale_y_continuous(breaks = c(0, 10, 20, 30, 40, 50)) +
  scale_fill_colorblind() +
  coord_equal() +
  ylab('Crown closure index\n') +
  xlab('\nHeight (m)') +
  labs(fill = 'Tree') +
  theme_cowplot() +
  theme(legend.title = element_text(size = 7, face = 'bold'),
        legend.text = element_text(size = 7),
        axis.title.x = element_text(face = 'bold', size = 7),
        axis.text.x = element_text(size = 7),
        axis.title.y = element_text(face = 'bold', size = 7),
        axis.text.y = element_text(size = 7))

# HJA raster map ####
# Load and trim the appropriate rasters
rasters <- list.files(here('data', 'rasters'),
                      full.names = T, recursive = T,
                      pattern = '^hillshade.bare|^elevation')

rasters %<>% lapply(rast)
rasters <- c(rasters[[1]], rasters[[2]]) %>%
  terra::aggregate(fact = 5)

# Make a box approximating the location Andrews
hja.ext <- ext(rasters)
hja.ext <- hja.ext@ptr$vector
hja.point <- c(median(hja.ext[1:2]), median(hja.ext[3:4]))

# Make the raster plot of HJA
# terra::as.data.frame currently gives many errors, but the output is still fine...
rasters.df <- terra::as.data.frame(rasters, xy = T)

# set seed to control point jitter
set.seed(123)
hja.map <- ggplot(mapping = aes(x = x, y = y)) +
  geom_raster(data = rasters.df, aes(fill = hillshade.bare), alpha = 0.5, show.legend = F) +
  scale_fill_gradient(low = 'black', high = 'white') +
  new_scale_fill() +
  geom_raster(data = rasters.df, aes(fill = elevation), alpha = 0.8) +
  scale_fill_gradient(low = 'grey20', high = 'white', name = 'Elevation (m)\n',
                      breaks = c(500, 750, 1000, 1250, 1500, 1750)) +
  new_scale_fill() +
  geom_jitter(data = tree.tab, aes(x = x, y = y,
                                   fill = tree),
              alpha = 0.95, size = 1.5, shape = 21,
              width = 150, height = 150, show.legend = F) +
  # geom_point(data = tree.tab, aes(x = x, y = y,
  #                                  fill = tree),
  #             alpha = 0.95, size = 1.5, shape = 21,
  #             show.legend = T) +
  scale_fill_colorblind(name = 'Tree') +
  scale_x_continuous(n.breaks = 8) +
  scale_y_continuous(n.breaks = 6) +
  coord_fixed() +
  xlab('\nEasting (m)') +
  ylab('Northing (m)\n') +
  theme_cowplot() +
  theme(axis.text.x = element_text(size = 7, angle = 45, hjust = 1, vjust = 1),
        axis.title.x = element_text(size = 7, face = 'bold'),
        axis.text.y = element_text(size = 7, angle = 45, hjust = 1, vjust = 1),
        axis.title.y = element_text(size = 7, face = 'bold'),
        legend.key.height = unit(1, 'line'),
        legend.key.width = unit(0.5, 'line'),
        legend.position = 'right',
        legend.direction = 'vertical',
        legend.title = element_text(face = 'bold', size = 7),
        legend.text = element_text(size = 7),
        panel.border = element_rect(color = 'black'))

# PNW (sensu lato) shapefile map ####
us <- read_sf(here('data', 'shapefiles', 'cb_2018_us_state_5m', 'cb_2018_us_state_5m.shp'))
us %<>% dplyr::select(st.prov = STUSPS, geometry)

ca <- read_sf(here('data', 'shapefiles', 'gpr_000b11a_e', 'gpr_000b11a_e.shp'))
ca %<>% dplyr::select(st.prov = PREABBR, geometry)

na <- bind_rows(us, ca)
pnw <- c('B.C.', 'Alta.', 'WA', 'ID', 'MT', 'WA', 'OR', 'CA', 'NV')
pnw <- filter(na, st.prov %in% pnw)

pnw %<>% st_transform(26910)

pnw.map <- ggplot(data = pnw) +
  geom_sf(color = 'black', fill = 'white', size = 0.1) +
  geom_point(x = hja.point[1], y = hja.point[2], color = 'red', size = 0.5,
             shape = 3, stroke = 0.5) +
  theme_map()

map.plot <- hja.map + inset_element(pnw.map, left = 0.13, bottom = 0.55, right = 0.45, top = 1, align_to = 'full')

tree.map <-  map.plot / (tree.plot + crown.vars) +
  plot_annotation(tag_levels = list(c('A', '', 'B', 'C'))) & 
  theme(plot.tag = element_text(size = 10, face = 'bold'))

ggsave(here(figure.out, 'fig.1.tiff'), tree.map, units = 'mm', width = 190, height = 190,
       dpi = 500, compression = 'lzw')

# Get session info ####
session.path <- here('output', 'sessions')
dir.create(session.path, recursive = T)
sessionInfo() %>% saveRDS(here(session.path, 'tree.map.sesh.rds'))
