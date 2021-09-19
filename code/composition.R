# Creates Tables 1, which: ####
# shows the results of PERMANOVA testing whether compositions differ among
# trees, needle age classes, exposure groups, and along crown variables more than expected by chance.
# Non-age variables are nested in age, and crown variables underwent marginal (type II) testing.

# Creates Figure 4, which: ####
# shows the results of db-RDA where compositional variation is constrained onto 
# height, depth, and crown closure, separated by needle age class.

# Creates Figure 5, which: ####
# shows the results of partial Mantel tests of uncorrelated community structures for pairs of
# needle age class transitions, separated by exposure group.

# Creates Figure 6, which: ####
# shows the relative abundances of Nothophaeocryptopus gaeumannii (a) and Rhabdocline parkeri
# (b) across age classes and split by exposure group.

# Creates Table S2, which: ####
# shows the results of PERMDISP2 tests of homogeneity of variance among trees, needle age classes, and exposure groups.

# Creates Tables S3-S5, which: ####
# show the results of PERMANOVAs testing whether compositions within individual trees vary along height
# and/or closure gradients. This involves a marginal test (S3), and tests only examining variation along
# height (S4) or closure (S5).

# Creates Tables S6-S7, which: ####
# subsets the dataset by needle age class and shows the results of PERMANOVAs testing whether
# compositions differ along crown variables (S6, marginal/type II testing) or between exposure groups (S7)
# more than expected by chance.

# Creates Tables S8, which: ####
# identifies indicator taxa for open and closed exposure groups

# Creates Figure S5, which: #####
# shows unconstrained NMDS ordinations for all retained samples, with points colored by
# (a) tree of origin and
# (b) needle age class.

library(patchwork)
library(ggthemes)
library(cowplot)
library(gt)
library(labdsv)
library(phyloseq)
library(vegan)
library(tidyverse)
library(magrittr)
library(here)
# Consult "code/functions.R" for descriptions of function applications, inputs, and outputs
source(here('code', 'functions.R'))

# Set out-paths ####
figure.out <- here('output', 'figures')
dir.create(figure.out, recursive = T)
table.out <- here('output', 'tables')
dir.create(table.out, recursive = T)

# Read in files ####
# Run "clean.R" prior to running this script
perf.n <- readRDS(here('data', 'compile', 'perf.n.clean.rds'))

counts <- perf.n$counts
relative <- perf.n$relative
logged <- perf.n$logged
tree.log <- perf.n$tree.log

# Test of no compositional differences among trees, nested within needle age class ####
npa.age.tree.terms.all <- permanova(logged$all, form = 'dist ~ age/tree', by = 'terms', n.perm = 999)
npa.age.tree.terms.all$row.group <- 'Age class and tree'
npa.age.tree.terms.all$Tree <- NULL
npa.age.tree.terms.all$Age <- NULL

# Marginal test of no compositional difference along crown variables, nested within needle age ####
npa.age.crown.margin.all <- permanova(logged$all, form = 'dist ~ age/(height + depth + closure)', by = 'margin', n.perm = 999)
npa.age.crown.margin.all$row.group <- 'Age class and crown variables'
npa.age.crown.margin.all$Tree <- NULL
npa.age.crown.margin.all$Age <- NULL

# Test of no compositional difference between exposure groups, nested in needle age class ####
npa.group.all <- permanova(logged$all, form = 'dist ~ age/group', by = 'terms', n.perm = 999)
npa.group.all$row.group <- 'Age class and exposure group'
npa.group.all$Tree <- NULL
npa.group.all$Age <- NULL

# Combine all 'all' PERMANOVAs
rbind(npa.age.tree.terms.all, npa.age.crown.margin.all, npa.group.all) %>%
  gt(groupname_col = 'row.group') %>%
  tab_header(title = 'Table 1. PERMANOVA results displaying the amount of compositional variation accounted for by age class; tree; the marginal variation of height, depth, and crown closure; and exposure group.') %>%
  cols_label(R2 = html('R<sup>2</sup>')) %>%
  tab_style(
    style = cell_text(style = 'italic'),
    locations = cells_body(
      columns = vars(Term),
      rows = grepl('^[[:lower:]]', Term))) %>%
  fmt_missing(everything(), missing_text = '') %>%
  tab_source_note(source_note = 'P values are the results of permutation tests with
                  pseudo F-ratios and 999 iterations each.') %>%
  opt_table_font(font = google_font('Crimson Text')) %>%
  cols_align(align = 'center') %>% gtsave(here(table.out, 'table.1.png'))

# dbRDA constrained on crown variables, applied to each age class separately ####
db.cols <- c('#d5a6bd', '#674ea7', '#1155cc', '#000000', '#666666')

# Vectors for crown variables are extracted from the "biplot" output of the dbRDA
# OTU relative abundance vectors are fit with the "envfit" function
a1.db <- dbrda.it(logged$a1, relative$a1, form = 'otu.tab ~ height + depth + closure',
                  vectors = T, min.r = 0.3, max.p = 0.05) +
  geom_point(aes(fill = tree, shape = group),
             size = 2,
             alpha = 0.7) +
  scale_shape_manual(values = 21:22) +
  guides(fill = guide_legend(override.aes = list(shape = c(21)))) +
  scale_fill_colorblind() +
  scale_color_manual(values = db.cols) +
  theme(legend.position = 'none')
a1.db$layers <- rev(a1.db$layers)

a2.db <- dbrda.it(logged$a2, relative$a2, form = 'otu.tab ~ height + depth + closure',
                  vectors = T, min.r = 0.3, max.p = 0.05) +
  geom_point(aes(fill = tree, shape = group),
             size = 2,
             alpha = 0.7) +
  scale_shape_manual(values = 21:22) +
  guides(fill = guide_legend(override.aes = list(shape = c(21)))) +
  scale_fill_colorblind() +
  scale_color_manual(values = db.cols, labels = c('Closure', 'Depth', 'Height',
                                               expression(italic('N. gaeumannii')),
                                               expression(italic('R. parkeri')))) +
  labs(fill = 'Tree', shape = 'Exposure', color = 'Vectors')
a2.db$layers <- rev(a2.db$layers)

a3.db <- dbrda.it(logged$a3, relative$a3, form = 'otu.tab ~ height + depth + closure',
                  vectors = T, min.r = 0.3, max.p = 0.05) +
  geom_point(aes(fill = tree, shape = group),
             size = 2,
             alpha = 0.7) +
  scale_shape_manual(values = 21:22) +
  guides(fill = guide_legend(override.aes = list(shape = c(21)))) +
  scale_fill_colorblind() +
  scale_color_manual(values = db.cols) +
  theme(legend.position = 'none')
a3.db$layers <- rev(a3.db$layers)

a4.db <- dbrda.it(logged$a4, relative$a4, form = 'otu.tab ~ height + depth + closure',
                  vectors = T, min.r = 0.3, max.p = 0.05) +
  geom_point(mapping = aes(fill = tree, shape = group),
             size = 2,
             alpha = 0.7) +
  scale_shape_manual(values = 21:22) +
  guides(fill = guide_legend(override.aes = list(shape = c(21)))) +
  scale_fill_colorblind() +
  scale_color_manual(values = db.cols) +
  theme(legend.position = 'none')
a4.db$layers <- rev(a4.db$layers)

# Combine all the plots
ages.db <- (a1.db | a2.db) / (a3.db | a4.db) + plot_annotation(tag_levels = list(c('A\n\nA1',
                                                         'B\n\nA2',
                                                         'C\n\nA3',
                                                         'D\n\nA4'))) +
  plot_layout(guides = 'collect') &
  theme(plot.tag = element_text(size = 10, face = 'bold'))
ggsave(here(figure.out, 'fig.4.tiff'), ages.db, units = 'mm', width = 190, height = 140,
       dpi = 500, compression = 'lzw')

# Mantel tests between ages ####
# Normal, unbootstrapped tests, accounting for xyz positions of each sampled height
mant.p <- rbind(
  time.mantel(logged$a1, logged$a2, 'Open', spatial = T),
  time.mantel(logged$a2, logged$a3, 'Open', spatial = T),
  time.mantel(logged$a3, logged$a4, 'Open', spatial = T),
  time.mantel(logged$a1, logged$a2, 'Closed', spatial = T),
  time.mantel(logged$a2, logged$a3, 'Closed', spatial = T),
  time.mantel(logged$a3, logged$a4, 'Closed', spatial = T)
)

# Bootstrapped tests to obtain mean and confidence intervals for each Mantel correlation
mant.boot <- rbind(
  boot.it.mantel(logged$a1, logged$a2, group = 'Open'),
  boot.it.mantel(logged$a2, logged$a3, group = 'Open'),
  boot.it.mantel(logged$a3, logged$a4, group = 'Open'),
  boot.it.mantel(logged$a1, logged$a2, group = 'Closed'),
  boot.it.mantel(logged$a2, logged$a3, group = 'Closed'),
  boot.it.mantel(logged$a3, logged$a4, group = 'Closed')
)

# Combine both normal and bootstrapped tests
# Bootstrapped tests do not partial out spatial variation
time <- left_join(mant.boot, mant.p, by = c('Group', 'Transition'))
# Identify significant correlations
time$sig <- ifelse(time$p < 0.05, '*', '')

time.plot <- ggplot(time, aes(x = Transition, y = mean)) +
  geom_pointrange(mapping = aes(ymin = lci, ymax = uci, color = Group),
                  size = 0.35,
                  position = position_dodge(width = 0.3)) +
  geom_text(mapping = aes(label = sig, group = Group, y = mean - 0.0075),
            color = 'white',
            fontface = 'bold',
            size = 4,
            position = position_dodge(width = 0.3)) +
  ylab('Spearman ranked\ncorrelation\n') +
  xlab('\nTransition') +
  scale_color_colorblind() +
  theme_cowplot() +
  theme(axis.title.x = element_text(size = 7, face = 'bold'),
        axis.text.x = element_text(size = 7),
        axis.title.y = element_text(size = 7, face = 'bold'),
        axis.text.y = element_text(size = 7),
        legend.title = element_text(size = 7, face = 'bold'),
        legend.text = element_text(size = 7))
time.plot
ggsave(here(figure.out, 'fig.5.tiff'), time.plot, units = 'mm', width = 90,
       dpi = 500, compression = 'lzw')

# Comparing the relative abundance of NOGA between exposure groups at each age ####
# Nothophaeocryptopus gaeumannii (OTU.1)
ng.ra <- relative$all@otu_table %>% data.frame() %>% .$otu.1
ng.sam.data <- relative$all@sam_data %>% data.frame()
ng.sam.data$ng.ra <- ng.ra

set.seed(666)
ng.ra.plot <- ggplot(ng.sam.data, aes(x = group, y = ng.ra, fill = tree)) +
  facet_grid(cols = vars(age)) +
  geom_violin(fill = 'white') +
  ggbeeswarm::geom_quasirandom(width = 0.4, shape = 21, size = 1, show.legend = F) +
  scale_fill_colorblind() +
  xlab('\nExposure group') +
  ylab(expression(atop(bolditalic(N.~gaeumannii), bold(relative~abundance)))) +
  labs(fill = 'Tree') +
  theme_cowplot() +
  theme(panel.spacing = unit(1, 'lines'),
        strip.text = element_text(size = 7),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        legend.title = element_text(face = 'bold'))

# Rhabdocline parkeri (OTU.2, OTU.6)
# Combine the relative abundances of the two indicator taxa we identified
otu.2.ra <- relative$all@otu_table %>% data.frame() %>% .$otu.2
otu.6.ra <- relative$all@otu_table %>% data.frame() %>% .$otu.6
rp.ra <- otu.2.ra + otu.6.ra
rp.sam.data <- relative$all@sam_data %>% data.frame()
rp.sam.data$rp.ra <- rp.ra

set.seed(666)
rp.ra.plot <- ggplot(rp.sam.data, aes(x = group, y = rp.ra, fill = tree)) +
  facet_grid(cols = vars(age)) +
  geom_violin(fill = 'white') +
  ggbeeswarm::geom_quasirandom(width = 0.4, shape = 21, size = 1) +
  scale_fill_colorblind() +
  xlab('\nExposure group') +
  ylab(expression(atop(bolditalic(R.~parkeri), bold(relative~abundance)))) +
  labs(fill = 'Tree') +
  theme_cowplot() +
  theme(panel.spacing = unit(1, 'lines'),
        strip.text = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_text(face = 'bold', size = 7),
        axis.text.x = element_text(size = 7, angle = 45, hjust = 1, vjust = 1),
        axis.title.y = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        legend.title = element_text(face = 'bold', size = 7),
        legend.text = element_text(size = 7))
ng.ra.plot / rp.ra.plot + plot_layout(guides = 'collect') + plot_annotation(tag_levels = list(c('A', 'B'))) &
  theme(plot.tag = element_text(size = 10, face = 'bold'))
ggsave(here(figure.out, 'fig.6.tiff'), units = 'mm', width = 140,
       dpi = 500, compression = 'lzw')

# Tests of equal variance among trees, age classes, exposure groups ####
pd.tree.all <- permdisp(logged$all, test = 'tree', n.perm = 999, bias = T)
pd.tree.all$Age <- NULL

pd.age.all <- permdisp(logged$all, test = 'age', n.perm = 999)
pd.age.all$Age <- NULL

pd.group.all <- permdisp(logged$all, test = 'group', n.perm = 999, bias = T)
pd.group.all$Age <- NULL

rbind(pd.tree.all, pd.age.all, pd.group.all) %>%
  gt(rowname_col = 'row') %>%
  tab_header(title = 'Table S2. Results of PERMDISP2 tests of homogeneity of variance for tree, age and group variables.') %>%
  tab_style(
    style = cell_text(style = 'italic'),
    locations = cells_body(
      columns = vars(Term),
      rows = grepl('^[[:lower:]]', Term))) %>%
  opt_table_font(font = google_font('Crimson Text')) %>%
  fmt_missing(everything(), missing_text = '') %>%
  cols_align(align = 'center') %>%
  tab_source_note(source_note = 'P values are the results of permutation tests with pseudo F-ratios and 999 iterations each.') %>%
  gtsave(here(table.out, 'table.s2.png'))

# Vertical stratification within individual trees ####
# All trees together, nesting crown variables in age
# Split by tree, nesting crown variables within age, A1 removed
# Implemented similarly to Harrison et al. 2016

# Marginal test of crown variables
npa.age.crown.trees <- lapply(tree.log[2:9], permanova, form = 'dist ~ age/(height + closure)',
                              by = 'margin', n.perm = 999) %>% bind_rows()
npa.age.crown.trees$Age <- NULL
npa.age.crown.trees %>%
  gt(rowname_col = 'row', groupname_col = 'Tree') %>%
  tab_header(title = 'Table S3. PERMANOVA results (separated by tree) displaying the marginal compositional variation accounted for by crown variables.') %>%
  fmt_missing(everything(), missing_text = '') %>%
  cols_align(align = 'center') %>%
  opt_table_font(font = google_font('Crimson Text')) %>%
  tab_style(
    style = cell_text(style = 'italic'),
    locations = cells_body(
      columns = vars(Term),
      rows = grepl('^[[:lower:]]', Term))) %>%
  tab_source_note(source_note = 'P values are the results of permutation tests with pseudo F-ratios and 999 iterations each.') %>%
  gtsave(here(table.out, 'table.s3.png'))

# Height
npa.age.ht.trees <- lapply(tree.log[2:9], permanova, form = 'dist ~ age/height',
                           by = 'terms', n.perm = 999) %>% bind_rows()
npa.age.ht.trees$Age <- NULL
npa.age.ht.trees %>%
  gt(groupname_col = 'Tree') %>%
  tab_header(title = 'Table S4. PERMANOVA results (separated by tree) displaying the amount of compositional variation accounted for by height.') %>%
  fmt_missing(everything(), missing_text = '') %>%
  cols_align(align = 'center') %>%
  opt_table_font(font = google_font('Crimson Text')) %>%
  tab_style(
    style = cell_text(style = 'italic'),
    locations = cells_body(
      columns = vars(Term),
      rows = grepl('^[[:lower:]]', Term))) %>%
  tab_source_note(source_note = 'P values are the results of permutation tests with pseudo F-ratios and 999 iterations each.') %>%
  gtsave(here(table.out, 'table.s4.png'))

# Closure
npa.age.closure.trees <- lapply(tree.log[2:9], permanova, form = 'dist ~ age/closure',
                                by = 'terms', n.perm = 999) %>% bind_rows()
npa.age.closure.trees$Age <- NULL
npa.age.closure.trees %>%
  gt(groupname_col = 'Tree') %>%
  tab_header(title = 'Table S5. PERMANOVA results (separated by tree) displaying the amount of compositional variation accounted for by crown closure.') %>%
  fmt_missing(everything(), missing_text = '') %>%
  cols_align(align = 'center') %>%
  opt_table_font(font = google_font('Crimson Text')) %>%
  tab_style(
    style = cell_text(style = 'italic'),
    locations = cells_body(
      columns = vars(Term),
      rows = grepl('^[[:lower:]]', Term))) %>%
  tab_source_note(source_note = 'P values are the results of permutation tests with pseudo F-ratios and 999 iterations each.') %>%
  gtsave(here(table.out, 'table.s5.png'))

# Marginal tests of no compositional difference along crown variables, applied to each age class separately ####
npa.crown.margin.ages <- lapply(logged[2:5], permanova, form = 'dist ~ height + depth + closure', by = 'margin', n.perm = 999) %>%
  bind_rows()
npa.crown.margin.ages$Tree <- NULL

npa.crown.margin.ages %>%
  gt(groupname_col = c('Age')) %>%
  tab_header(title = 'Table S6. PERMANOVA results (separated by age class) displaying the amount of marginal compositional variation accounted for by crown variables.') %>%
  cols_label(R2 = html('R<sup>2</sup>')) %>%
  tab_style(
    style = cell_text(style = 'italic'),
    locations = cells_body(
      columns = vars(Term),
      rows = grepl('^[[:lower:]]', Term))) %>%
  fmt_missing(everything(), missing_text = '') %>%
  tab_source_note(source_note = 'P values are the results of permutation tests with
                  pseudo F-ratios and 999 iterations.') %>%
  opt_table_font(font = google_font('Crimson Text')) %>%
  cols_align(align = 'center') %>% gtsave(here(table.out, 'table.s6.png'))

# Tests of no compositional differences between exposure groups, applied to each age class separately ####
npa.group.ages <- lapply(logged[2:5], permanova, form = 'dist ~ group', by = 'terms', n.perm = 999) %>% bind_rows()
npa.group.ages$Tree <- NULL

npa.group.ages %>%
  gt(groupname_col = c('Age')) %>%
  tab_header(title = 'Table S7. PERMANOVA results (separated by age class) displaying the amount of compositional variation accounted for by exposure group.') %>%
  cols_label(R2 = html('R<sup>2</sup>')) %>%
  tab_style(
    style = cell_text(style = 'italic'),
    locations = cells_body(
      columns = vars(Term),
      rows = grepl('^[[:lower:]]', Term))) %>%
  fmt_missing(everything(), missing_text = '') %>%
  tab_source_note(source_note = 'P values are the results of permutation tests with pseudo F-ratios and 999 iterations.') %>%
  opt_table_font(font = google_font('Crimson Text')) %>%
  cols_align(align = 'center') %>% gtsave(here(table.out, 'table.s7.png'))

# Indicator species analysis ####
isa.group.all <- indicate(counts$all, test = 'group')
isa.group.all %<>% arrange(desc(Group), desc(IV))
isa.group.all$IV %<>% round(digits = 3)
isa.group.all$P %<>% round(digits = 3)
isa.group.all$Taxon %<>% str_replace('Nothophaeocryptopus', 'N.')
isa.group.all$Age <- NULL
isa.group.all %>%
  gt() %>%
  tab_header(title = 'Table S8. Results of indicator species analysis for open and closed exposure groups.') %>%
  cols_align(align = 'center') %>%
  opt_table_font(font = google_font('Crimson Text')) %>%
  tab_style(
    style = cell_text(style = 'italic'),
    locations = cells_body(
      columns = vars(Taxon),
      rows = !grepl('sp\\.', Taxon))) %>%
  tab_source_note(source_note = 'P values were obtained after 999 permutations and corrected
                    for the false discovery rate using the Benjamini-Hochberg method.') %>%
  gtsave(here(table.out, 'table.s8.png'))

# NMDS ordination of all needle ages ####
all <- nmds.mc.par(phy.in = logged$all, n.cores = 6, trymax = 75, perm = 99)

all.tree <- nmds.plot(all, note = F, 'NMDS1', 'NMDS2') +
  geom_point(shape = 21, size = 1.5, aes(fill = tree)) +
  scale_fill_colorblind() +
  labs(fill = 'Tree') +
  xlab('') +
  coord_fixed()
all.age <- nmds.plot(all, note = T, 'NMDS1', 'NMDS2') +
  geom_point(shape = 21, size = 1.5, aes(fill = age)) +
  scale_fill_colorblind() +
  labs(fill = 'Age') +
  coord_fixed()

all.nmds <- all.tree / all.age + plot_annotation(tag_levels = list(c('A', 'B'))) &
  theme(plot.tag = element_text(size = 10, face = 'bold'))
ggsave(here(figure.out, 'fig.s5.tiff'), all.nmds, units = 'mm', width = 140, height = 140,
       dpi = 500, compression = 'lzw')

# Get session info ####
session.path <- here('output', 'sessions')
dir.create(session.path, recursive = T)
sessionInfo() %>% saveRDS(here(session.path, 'comp.sesh.rds'))
