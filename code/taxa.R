# Create Figure S3 from the manuscript, which shows:
# a species accumulation curve for each tree

# Create Figure 2 from the manuscript, which shows:
# (a) a stacked bar plot, with each tree getting its own bar,
# (b) the species abundance distribution of OTUs across all eight trees, and
# (c) the occupancy-abundance relationship for each OTU across all trees.

library(scales)
library(cowplot)
library(patchwork)
library(ggthemes)
library(ggtext)
library(gt)
library(BiodiversityR)
library(phyloseq)
library(tidyverse)
library(magrittr)
library(here)
# Consult "code/functions.R" for descriptions of function application, inputs, and outputs
source(here('code', 'functions.R'))

# Read in phyloseq list
perf.n <- readRDS(here('data', 'compile', 'perf.n.clean.rds'))

# Set out-path
figure.out <- here('output', 'figures')
dir.create(figure.out, recursive = T)
table.out <- here('output', 'tables')
dir.create(table.out, recursive = T)

# Soecies abundance distribution and occupancy-abundance object creation ####
# Creates warnings associated with merging
soa.tree <- sad.occ.abund(perf.n$counts$all, var = 'tree')

# Species abundance distribution plot ####
sad.tree <- ggplot(soa.tree$sad, aes(x = samples, y = otus)) +
  geom_bar(stat = 'identity', width = 0.5, color = 'black', fill = 'black') +
  scale_x_continuous(n.breaks = 8) +
  scale_fill_viridis_c() +
  scale_y_continuous(n.breaks = 8) +
  # xlab('\nOTU occupancy class\n(# trees occupied)') +
  xlab('\nOTU occupancy class\n(# trees occupied)') +
  ylab('Proportion of OTUs in\noccupancy class\n') +
  theme_cowplot() +
  theme(axis.title.x = element_text(face = 'bold', size = 8),
        axis.title.y = element_text(face = 'bold', size = 8, margin = margin(l = 25)),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6))

# Occupancy-abundance plot ####
set.seed(666)
occ.abund.tree <- ggplot(soa.tree$otus, aes(x = pa, y = ra, group = pa)) +
  # geom_jitter(size = 1, shape = 21, width = 0.1, alpha = 0.4, fill = 'black') +
  # geom_violin(fill = 'white', draw_quantiles = c(0.25, 0.5, 0.75)) +
  ggbeeswarm::geom_quasirandom(width = 0.2, size = 1, color = 'black', alpha = 0.5) +
  scale_x_continuous(n.breaks = 8) +
  scale_fill_viridis_c() +
  coord_trans(y = 'log10') +
  scale_y_continuous(breaks = trans_breaks('log10', function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  xlab('\nOTU occupancy class\n(# trees occupied)') +
  ylab('Mean relative abundance\n') +
  theme_cowplot() +
  theme(axis.title.x = element_text(face = 'bold', size = 8),
        axis.title.y = element_text(face = 'bold', size = 8, margin = margin(l = 25)),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6))

# Create a stacked bar plot of the top 95% OTU relative abundances, with samples merged at the tree level ####
# Samples can be merged for any categorical variable, and any top percentile cap
stack <- rank.n.stack(perf.n$counts$all, 'tree', 0.95)

stack.plot <- ggplot(stack, aes(x = Sample, y = Abundance, fill = combo)) +
  geom_bar(stat = 'identity', width = 0.75, color = NA, show.legend = T) +
  ylab('Relative abundance') +
  xlab('') +
  ylim(c(0, 1)) +
  labs(fill = 'Genus | Order') +
  scale_fill_manual(values = cols, na.value = 'black') +
  theme_cowplot() +
  theme(axis.title.y = element_text(face = 'bold', size = 8),
        axis.text.x = element_text(size = 8, angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 6),
        legend.title = element_text(size = 7, face = 'bold'),
        legend.key.size = unit(0.5, units = 'lines'),
        legend.text = element_markdown(size = 7))

stack.plot / (sad.tree + occ.abund.tree) +
  plot_annotation(tag_levels = list(c('A', 'B', 'C'))) &
  theme(plot.tag = element_text(size = 10, face = 'bold'))
ggsave(here(figure.out, 'fig.2.tiff'), units = 'mm', width = 190, height = 175,
       dpi = 500, compression = 'lzw')

# Macrofungal survey ####
# The authors found this table intriguing but difficult to incorporate into the manuscript.
# However, it will still be retained on the GitHub repository
macro <- subset_taxa(perf.n$relative$all, Family == 'Cudoniaceae' |
                      Genus %in% c('Geopyxis', 'Tricharina', 'Wilcoxina', 'Pseudoplectania',
                                   'Nemania', 'Corticium', 'Ganoderma', 'Perenniporia', 'Trichaptum',
                                   'Lenzites', 'Trametes', 'Coprinellus', 'Stereum', 'Calyptella', 'Xylaria'))

macro.otu <- macro@otu_table %>% data.frame()

macro.samples <- macro.otu
macro.samples[macro.samples > 0] <- 1
macro.samples <- colSums(macro.samples)

macro.max <- sapply(macro.otu, max) %>% round(digits = 4)

macro.tax <- macro@tax_table %>% data.frame()
macro.tax %<>% rownames_to_column(var = 'OTU')
macro.tax$OTU %<>% str_replace('otu', 'OTU')
macro.tax$`Maximum relative abundance` <- macro.max
macro.tax$`Number of samples` <- macro.samples
macro.tax$Kingdom <- NULL
macro.tax$Genus_species <- NULL
macro.tax$Family %<>% str_replace('Hymenochaetales_fam_Incertae_sedis', 'Incertae sedis')
macro.tax %<>% arrange(desc(`Maximum relative abundance`)) %>%
  arrange(desc(`Number of samples`)) %>%
  arrange(Order, Family, Genus)
macro.tax %>%
gt(groupname_col = c('Phylum')) %>%
  tab_header(title = 'Putatively endophytic macrofungal OTUs, applying the definition of Thiers and Halling 2018.') %>%
  cols_align(align = 'center') %>%
  tab_style(
    style = cell_text(style = 'italic'),
    locations = cells_body(
      columns = vars(Genus, Species))) %>%
  fmt_missing(columns = everything(), missing_text = '') %>%
  tab_source_note(source_note = 'Sorted first by maximum relative abundance, by occurrence in samples, and then by order, family, and genus.') %>%
  opt_table_font(font = google_font('Crimson Text')) %>%
  gtsave(filename = here(table.out, 'macrofungi.png'))

# OTU accumulation curve ####
otu.tab <- perf.n$counts$all@otu_table %>% data.frame()
sam.data <- perf.n$counts$all@sam_data %>% data.frame()

otu.accum <- accumcomp(otu.tab, y = sam.data, factor = 'tree',
                        method = 'exact', conditioned = FALSE, plotit = FALSE) %>% 
  accumcomp.long(ci = NA)

ggplot(otu.accum, aes(x = Sites, y = Richness, ymax = UPR, ymin = LWR, color = Grouping)) +
  geom_line(size = 1, alpha = 0.75) +
  geom_ribbon(aes(fill = Grouping), size = 0.25, alpha = 0.15, show.legend = F) +
  scale_color_colorblind() +
  scale_fill_colorblind() +
  xlab('\nSamples') +
  ylab('Estimated richness') +
  labs(color = 'Tree') +
  theme_cowplot() +
  theme(axis.title.x = element_text(face = 'bold', size = 8),
        axis.title.y = element_text(face = 'bold', size = 8),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        legend.title = element_text(size = 7, face = 'bold'),
        legend.key.size = unit(0.5, units = 'lines'),
        legend.text = element_markdown(size = 7))
ggsave(here(figure.out, 'fig.s3.tiff'), units = 'mm', width = 140,
       dpi = 500, compression = 'lzw')

# Get session info ####
session.path <- here('output', 'sessions')
dir.create(session.path, recursive = T)
sessionInfo() %>% saveRDS(here(session.path, 'taxa.sesh.rds'))
