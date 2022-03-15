# Creates Table 1, which: ####
# shows the results of PERMANOVA testing whether compositions differ among
# trees, needle age classes, and along crown variables more than expected by chance.
# Non-age variables are nested in age, and crown variables underwent marginal (type II) testing.

# Creates Fig. 3, which: ####
# shows the results of db-RDA where compositional variation is constrained onto 
# crown closure, separated by needle age class.

# Creates Fig. 5, which: ####
# shows the results of partial Mantel tests of uncorrelated community structures for pairs of
# needle age class transitions, separated by exposure group.

# Creates Fig. 6, which: ####
# shows the relative abundances of indicator OTUs for Nothophaeocryptopus gaeumannii (A) and
# Rhabdocline parkeri (B) across age classes and split by exposure group.

# Creates Table S1, which: ####
# shows the results of PERMDISP2 tests of homogeneity of variance among trees and needle age classes.

# Creates Table S2, which: ####
# subsets the dataset by needle age class and shows the results of PERMANOVAs testing whether
# compositions differ with crown closure more than expected by chance.

# Creates Table S3, which: ####
# identifies indicator taxa for open and closed exposure groups

# Creates Fig. S2, which: #####
# shows unconstrained NMDS ordinations for all retained samples, with points colored by
# (A) tree of origin and
# (B) needle age class.

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
source(here('code', 'functions.R')
       )

# Set out-paths ####
figure.out <- here('output', 'figures')
dir.create(figure.out, recursive = T)
table.out <- here('output', 'tables')
dir.create(table.out, recursive = T)

# Read in files ####
# Run "clean.R" prior to running this script
perf.n <- readRDS(here('data', 'compile', 'perf.n.clean.rds')
                  )
counts <- perf.n$counts
relative <- perf.n$relative
logged <- perf.n$logged

# Test of no compositional differences among trees and needle age class ####
npa.age.tree.terms.all <- permanova(logged$all, form = 'dist ~ tree*age', by = 'terms', n.perm = 999, block = 'tree')
npa.age.tree.terms.all$row.group <- 'Age class and tree'
npa.age.tree.terms.all$Tree <- NULL
npa.age.tree.terms.all$Age <- NULL

# Marginal test of no compositional difference along crown variables, nested within needle age (log-transformed relative abundance) ####
npa.age.crown.margin.all <- permanova(logged$all, form = 'dist ~ age/(height + closure)', by = 'margin', n.perm = 999, block = 'tree')
npa.age.crown.margin.all$row.group <- 'Age class and crown variables'
npa.age.crown.margin.all$Tree <- NULL
npa.age.crown.margin.all$Age <- NULL

# Marginal test of no compositional difference between exposure groups, nested within needle age (log-transformed relative abundance) ####
permanova(logged$all, form = 'dist ~ age/group', by = 'margin', n.perm = 999, block = 'tree')

# Marginal test of no compositional difference between exposure groups, nested within needle age (relative abundance) ####
permanova(relative$all, form = 'dist ~ age/group', by = 'margin', n.perm = 999, block = 'tree')

# Marginal test of no compositional difference along crown variables, nested within needle age (relative abundance) ####
permanova(relative$all, form = 'dist ~ age/(height + closure)', by = 'margin', n.perm = 999, block = 'tree')

# Combine all 'all' PERMANOVAs
rbind(npa.age.tree.terms.all, npa.age.crown.margin.all) %>%
  gt(groupname_col = 'row.group') %>%
  tab_header(title = 'Table 1. PERMANOVA results displaying the amount of compositional variation accounted for by needle age class and tree, and the marginal variation of height and crown closure.',
             ) %>%
  cols_label(R2 = html('R<sup>2</sup>')) %>%
  tab_options(
    heading.title.font.size = px(18),
    heading.align = 'left',
    heading.border.bottom.color = 'white',
    table.border.top.color = 'white',
    table.border.bottom.color = 'white',
    table_body.border.top.color = 'white',
    table_body.border.bottom.color = 'white',
    row_group.border.top.color = 'white',
    row_group.border.bottom.color = 'white',
    table_body.hlines.color = 'white'
    ) %>%
  tab_style(
    style = cell_text(weight = 'bold'),
    locations = cells_row_groups(
      groups = everything()
    )
  ) %>%
  tab_style(
    style = cell_text(style = 'italic'),
    locations = cells_body(
      columns = 'Term',
      rows = grepl('^[[:lower:]]', Term)
    )
  ) %>%
  tab_style(
    style = cell_borders(
      sides = 'top',
      color = '#000000',
      weight = px(1.5),
      style = 'solid'),
    locations = cells_body(
      columns = everything(),
      rows = Term == 'Total')
    ) %>%
  tab_style(
    style = cell_text(weight = 'bold'),
    locations = cells_column_labels(everything())) %>%
  tab_style(
    style = cell_text(color = 'white'),
    locations = cells_column_labels(columns = 'Term')) %>%
  fmt_missing(everything(), missing_text = '') %>%
  tab_source_note(source_note = 'P values are the results of permutation tests with
                  pseudo F-ratios and 999 iterations each, setting "tree" as a stratum.') %>%
  opt_table_font(font = google_font('Crimson Text')) %>%
  cols_align(align = 'center') %>%
  cols_width(everything()~px(100)) %>% 
  gtsave(here(table.out, 'table.1.png'))

# dbRDA constrained on crown variables, applied to each age class separately ####
db.cols <- c("#5a934b",
              "#9756a8",
              "#620000",
              "#d18b14",
              "#0080ff",
              "#bb0e03")
# Vectors for crown variables are extracted from the "biplot" output of the dbRDA
# OTU vectors are fit with the "envfit" function
# A1
a1.db <- dbrda.it(logged$a1, relative$a1, form = 'otu.tab ~ height + closure',
                  vectors = T, min.r = 0, max.p = 0.05) +
  geom_point(aes(fill = tree, shape = group),
             size = 2,
             alpha = 0.7) +
  scale_shape_manual(values = 21:22) +
  guides(fill = guide_legend(override.aes = list(shape = c(21)
                                                 )
                             )
         ) +
  scale_fill_colorblind() +
  scale_color_manual(values = db.cols) +
  theme(legend.position = 'none')
a1.db$layers <- rev(a1.db$layers)
a1.db

# A2
a2.db <- dbrda.it(logged$a2, relative$a2, form = 'otu.tab ~ height + closure',
                  vectors = T, min.r = 0, max.p = 0.05) +
  geom_point(aes(fill = tree, shape = group),
             size = 2,
             alpha = 0.7) +
  scale_shape_manual(values = 21:22) +
  guides(fill = guide_legend(override.aes = list(shape = c(21)
                                                 )
                             )
         ) +
  scale_fill_colorblind() +
  scale_color_manual(values = db.cols, labels = c('Closure',
                                                  'Height',
                                               expression(paste('OTU.1 - ', italic('N. gaeumannii')
                                                                )
                                                          ),
                                               expression(paste('OTU.2 - ', italic('R. parkeri')
                                                                )
                                                          ),
                                               expression(paste('OTU.3 - ', italic('R. parkeri')
                                                                )
                                                          ),
                                               expression(paste('OTU.6 - ', italic('R. parkeri')
                                                                )
                                                          )
                                               )
                     ) +
  labs(fill = 'Tree', shape = 'Exposure', color = 'Vectors')
a2.db$layers <- rev(a2.db$layers)
a2.db

# A3
a3.db <- dbrda.it(logged$a3, relative$a3, form = 'otu.tab ~ height + closure',
                  vectors = T, min.r = 0, max.p = 0.05) +
  geom_point(aes(fill = tree, shape = group),
             size = 2,
             alpha = 0.7) +
  scale_shape_manual(values = 21:22) +
  guides(fill = guide_legend(override.aes = list(shape = c(21)
                                                 )
                             )
         ) +
  scale_fill_colorblind() +
  scale_color_manual(values = db.cols) +
  theme(legend.position = 'none')
a3.db$layers <- rev(a3.db$layers)
a3.db

# A4
a4.db <- dbrda.it(logged$a4, relative$a4, form = 'otu.tab ~ height + closure',
                  vectors = T, min.r = 0, max.p = 0.05) +
  geom_point(mapping = aes(fill = tree, shape = group),
             size = 2,
             alpha = 0.7) +
  scale_shape_manual(values = 21:22) +
  guides(fill = guide_legend(override.aes = list(shape = c(21)
                                                 )
                             )
         ) +
  scale_fill_colorblind() +
  scale_color_manual(values = db.cols) +
  theme(legend.position = 'none')
a4.db$layers <- rev(a4.db$layers)
a4.db

# Combine all the plots
ages.db <- (a1.db | a2.db) / (a3.db | a4.db) + plot_annotation(tag_levels = list(c('A\n\nA1',
                                                                                   'B\n\nA2',
                                                                                   'C\n\nA3',
                                                                                   'D\n\nA4')
                                                                                 )
                                                               ) +
  plot_layout(guides = 'collect') &
  theme(plot.tag = element_text(size = 10, face = 'bold')
        )
ages.db
ggsave(here(figure.out, 'fig.3.tiff'), ages.db, units = 'mm', width = 190, height = 140,
       dpi = 300, compression = 'lzw')

# Relative abundance dbRDA, constrained on closure, for the third age class
a3.db <- dbrda.it(relative$a3, relative$a3, form = 'otu.tab ~ height + closure',
                  vectors = T, min.r = 0, max.p = 0.05) +
  geom_point(aes(fill = tree, shape = group),
             size = 2,
             alpha = 0.7) +
  scale_shape_manual(values = 21:22) +
  guides(fill = guide_legend(override.aes = list(shape = c(21)
  )
  )
  ) +
  scale_fill_colorblind() +
  scale_color_manual(values = db.cols, labels = c('Closure',
                                                  'Height',
                                                  expression(paste('OTU.1 - ', italic('N. gaeumannii')
                                                  )
                                                  ),
                                                  expression(paste('OTU.2 - ', italic('R. parkeri')
                                                  )
                                                  ),
                                                  expression(paste('OTU.3 - ', italic('R. parkeri')
                                                  )
                                                  ),
                                                  expression(paste('OTU.6 - ', italic('R. parkeri')
                                                  )
                                                  )
  )
  ) +
  labs(fill = 'Tree', shape = 'Exposure', color = 'Vectors')
a3.db$layers <- rev(a3.db$layers)
a3.db

# Mantel tests between ages ####
# Normal, un-bootstrapped tests, accounting for xyz positions of each sampled height
mant.p <- rbind(
  time.mantel(logged$a1, logged$a2, 'Open'),
  time.mantel(logged$a2, logged$a3, 'Open'),
  time.mantel(logged$a3, logged$a4, 'Open'),
  time.mantel(logged$a1, logged$a2, 'Closed'),
  time.mantel(logged$a2, logged$a3, 'Closed'),
  time.mantel(logged$a3, logged$a4, 'Closed')
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
time <- left_join(mant.boot, mant.p, by = c('Group', 'Transition')
                  )
# Identify significant correlations
sig <- time %>% filter(p < 0.05)

# Plot sample and bootstrapped correlations
time.plot <- ggplot(time, aes(x = Transition, y = mean)
                    ) +
  # geom_point(aes(y = r),
  #            size = 3,
  #            shape = 10,
  #            stroke = 0.3,
  #            color = 'red') +
  geom_pointrange(aes(ymin = lci, ymax = uci),
                  size = 0.5,
                  shape = 3,
                  stroke = 0.5,
                  position = position_dodge(width = 0.3)
                  ) +
  geom_point(data = sig,
             size = 2,
             shape = 8,
             stroke = 0.5,
             color = 'blue') +
  facet_grid(rows = vars(Group)
             ) +
  ylab('Spearman ranked\ncorrelation\n') +
  xlab('\nTransition') +
  theme_cowplot() +
  theme(axis.title.x = element_text(size = 7, face = 'bold'),
        axis.text.x = element_text(size = 7),
        axis.title.y = element_text(size = 7, face = 'bold'),
        axis.text.y = element_text(size = 7),
        legend.title = element_text(size = 7, face = 'bold'),
        legend.text = element_text(size = 7),
        strip.text.y = element_text(size = 7, face = 'bold'),
        plot.background = element_rect(fill = 'white', color = 'white')
        )
time.plot
ggsave(here(figure.out, 'fig.5.tiff'), time.plot, units = 'mm', width = 90,
       dpi = 300, compression = 'lzw')

# Comparing the relative abundance of NOGA between exposure groups at each age ####
# Nothophaeocryptopus gaeumannii (OTU.1)
ng.ra <- relative$all@otu_table %>% data.frame() %>% .$otu.1
ng.sam.data <- relative$all@sam_data %>% data.frame()
ng.sam.data$ng.ra <- ng.ra

ng.ra.a1 <- ng.sam.data %>% filter(age == 'A1') %>% arrange(tree.ht)
ng.ra.a2 <- ng.sam.data %>% filter(age == 'A2') %>% arrange(tree.ht)
ng.ra.a3 <- ng.sam.data %>% filter(age == 'A3') %>% arrange(tree.ht)
ng.ra.a4 <- ng.sam.data %>% filter(age == 'A4') %>% arrange(tree.ht)

wilcox.test(ng.ra.a1$ng.ra, ng.ra.a2$ng.ra, paired = T)
wilcox.test(ng.ra.a2$ng.ra, ng.ra.a3$ng.ra, paired = T)
wilcox.test(ng.ra.a3$ng.ra, ng.ra.a4$ng.ra, paired = T)

wilcox.test(ng.ra.a3[ng.ra.a3$group == 'Open', ]$ng.ra, ng.ra.a3[ng.ra.a3$group == 'Closed', ]$ng.ra, exact = F)

set.seed(666)
ng.ra.plot <- ggplot(ng.sam.data, aes(x = group, y = ng.ra, fill = tree)
                     ) +
  facet_grid(cols = vars(age)
             ) +
  geom_violin(fill = 'white') +
  ggbeeswarm::geom_quasirandom(width = 0.4, shape = 21, size = 1, show.legend = F) +
  scale_fill_colorblind() +
  xlab('\nExposure group') +
  ylab(expression(atop(bolditalic(N.~gaeumannii), bold(relative~abundance)
                       )
                  )
       ) +
  labs(fill = 'Tree') +
  theme_cowplot() +
  theme(panel.spacing = unit(1, 'lines'),
        strip.text = element_text(size = 7),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        legend.title = element_text(face = 'bold')
        )
ng.ra.plot

# Rhabdocline parkeri (OTU.2, OTU.6)
# Combine the relative abundances of the two indicator taxa we identified
otu.ra <- relative$all@otu_table %>% data.frame()
otu.2.ra <- otu.ra$otu.2
otu.3.ra <- otu.ra$otu.3
otu.6.ra <- otu.ra$otu.6
rp.ra <- otu.2.ra + otu.3.ra + otu.6.ra
# rp.ra <- otu.2.ra + otu.3.ra
rp.sam.data <- relative$all@sam_data %>% data.frame()
rp.sam.data$rp.ra <- rp.ra

rp.ra.a1 <- rp.sam.data %>% filter(age == 'A1') %>% arrange(tree.ht)
rp.ra.a2 <- rp.sam.data %>% filter(age == 'A2') %>% arrange(tree.ht)
rp.ra.a3 <- rp.sam.data %>% filter(age == 'A3') %>% arrange(tree.ht)
rp.ra.a4 <- rp.sam.data %>% filter(age == 'A4') %>% arrange(tree.ht)

wilcox.test(rp.ra.a1$rp.ra, rp.ra.a2$rp.ra, paired = T)
wilcox.test(rp.ra.a2$rp.ra, rp.ra.a3$rp.ra, paired = T)
wilcox.test(rp.ra.a3$rp.ra, rp.ra.a4$rp.ra, paired = T)

wilcox.test(rp.ra.a2[rp.ra.a2$group == 'Open', ]$rp.ra, rp.ra.a2[rp.ra.a2$group == 'Closed', ]$rp.ra, exact = T)
wilcox.test(rp.ra.a3[rp.ra.a3$group == 'Open', ]$rp.ra, rp.ra.a3[rp.ra.a3$group == 'Closed', ]$rp.ra, exact = T)

set.seed(666)
rp.ra.plot <- ggplot(rp.sam.data, aes(x = group, y = rp.ra, fill = tree)
                     ) +
  facet_grid(cols = vars(age)
  ) +
  geom_violin(fill = 'white') +
  ggbeeswarm::geom_quasirandom(width = 0.4, shape = 21, size = 1) +
  scale_fill_colorblind() +
  xlab('\nExposure group') +
  ylab(expression(atop(bolditalic(R.~parkeri), bold(relative~abundance)
                       )
                  )
       ) +
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
        legend.text = element_text(size = 7)
        )
rp.ra.plot
ng.ra.plot / rp.ra.plot + plot_layout(guides = 'collect') + plot_annotation(tag_levels = list(c('A', 'B')
                                                                                              )
                                                                            ) &
  theme(plot.tag = element_text(size = 10, face = 'bold')
        )
ggsave(here(figure.out, 'fig.6.tiff'), units = 'mm', width = 140,
       dpi = 300, compression = 'lzw')

# Tests of equal variance among trees, age classes, exposure groups ####
pd.tree.all <- permdisp(logged$all, test = 'tree', n.perm = 999, bias = T, block = 'tree')
pd.tree.all$test <- 'test.1'

pd.age.all <- permdisp(logged$all, test = 'age', n.perm = 999, block = 'tree')
pd.age.all$test <- 'test.2'

permdisp(logged$all, test = 'group', n.perm = 999, bias = T, block = 'tree')

rbind(pd.tree.all, pd.age.all) %>%
  gt(groupname_col = 'test') %>%
  tab_header(title = 'Table S1. Results of PERMDISP2 tests of homogeneity of variance for tree and needle age class.') %>%
  tab_options(
    heading.title.font.size = px(18),
    heading.align = 'left',
    heading.border.bottom.color = 'white',
    table.border.top.color = 'white',
    table.border.bottom.color = 'white',
    table_body.border.top.color = 'white',
    table_body.border.bottom.color = 'white',
    row_group.border.top.color = 'white',
    row_group.border.bottom.color = 'white',
    table_body.hlines.color = 'white'
  ) %>%
  tab_style(
    style = cell_text(style = 'italic'),
    locations = cells_body(
      columns = 'Term',
      rows = grepl('^[[:lower:]]', Term)
      )
    ) %>%
  tab_style(
    style = cell_borders(
      sides = 'top',
      color = 'black',
      weight = px(1.5),
      style = 'solid'),
    locations = cells_body(
      columns = everything(),
      rows = Term == 'Total')
    ) %>%
  tab_style(
    style = cell_text(color = 'white'),
    locations = cells_column_labels(columns = 'Term')
    ) %>%
  tab_style(
    style = cell_text(color = 'white'),
    locations = cells_row_groups(everything()
                                 )
  ) %>%
  tab_style(
    style = cell_text(weight = 'bold'),
    locations = cells_column_labels(everything()
                                    )
    ) %>%
  opt_table_font(font = google_font('Crimson Text')
                 ) %>%
  fmt_missing(everything(), missing_text = '') %>%
  cols_align(align = 'center') %>%
  cols_width(everything()~px(75)
             ) %>% 
  tab_source_note(source_note = 'P values are the results of permutation tests with pseudo F-ratios and 999 iterations each, setting "tree" as a stratum.') %>%
  gtsave(here(table.out, 'table.s1.png')
         )

# Marginal tests of no compositional difference along crown variables, applied to each age class separately (log-transformed relative abundance) ####
npa.crown.margin.ages <- lapply(logged[2:5], permanova, form = 'dist ~ height + closure', by = 'margin', n.perm = 999, block = 'tree') %>%
  bind_rows()
npa.crown.margin.ages$Tree <- NULL

npa.crown.margin.ages %>%
  gt(groupname_col = c('Age')) %>%
  tab_header(title = 'Table S2. PERMANOVA results displaying the amount of compositional variation accounted for by closure in each needle age class.') %>%
  cols_label(R2 = html('R<sup>2</sup>')
             ) %>%
  tab_options(
    heading.title.font.size = px(18),
    heading.align = 'left',
    heading.border.bottom.color = 'white',
    table.border.top.color = 'white',
    table.border.bottom.color = 'white',
    table_body.border.top.color = 'white',
    table_body.border.bottom.color = 'white',
    row_group.border.top.color = 'white',
    row_group.border.bottom.color = 'white',
    table_body.hlines.color = 'white'
  ) %>%
  tab_style(
    style = cell_text(weight = 'bold'),
    locations = cells_row_groups(
      groups = everything()
    )
  ) %>%
  tab_style(
    style = cell_text(style = 'italic'),
    locations = cells_body(
      columns = 'Term',
      rows = grepl('^[[:lower:]]', Term)
      )
    ) %>%
  tab_style(
    style = cell_text(color = 'white'),
    locations = cells_column_labels(columns = 'Term')
  ) %>%
  tab_style(
    style = cell_text(weight = 'bold'),
    locations = cells_column_labels(everything()
                                    )
    ) %>%
  tab_style(
    style = cell_borders(
      sides = 'top',
      color = 'black',
      weight = px(1.5),
      style = 'solid'),
    locations = cells_body(
      columns = everything(),
      rows = Term == 'Total')
    ) %>%
  fmt_missing(everything(), missing_text = '') %>%
  tab_source_note(source_note = 'P values are the results of permutation tests with
                  pseudo F-ratios and 999 iterations, setting "tree" as a stratum.') %>%
  opt_table_font(font = google_font('Crimson Text')) %>%
  cols_align(align = 'center') %>%
  cols_width(everything()~px(100)) %>% 
  gtsave(here(table.out, 'table.s2.png'))

# Tests of no compositional difference between exposure groups, applied to each needle age class separately (log-transformed relative abundance) ####
lapply(logged[2:5], permanova, form = 'dist ~ group', by = 'terms', n.perm = 999, block = 'tree') %>%
  bind_rows()

# Tests of no compositional difference between exposure groups, applied to each needle age class separately (relative abundance) ####
lapply(logged[2:5], permanova, form = 'dist ~ group', by = 'terms', n.perm = 999, block = 'tree') %>%
  bind_rows()

# Tests of no compositional difference along crown closure, applied to each needle age class separately (relative abundance) ####
lapply(relative[2:5], permanova, form = 'dist ~ height + closure', by = 'margin', n.perm = 999, block = 'tree') %>%
  bind_rows()

# Indicator species analysis ####
isa.group.all <- indicate(counts$all, test = 'group')
isa.group.all %<>% arrange(desc(Group), desc(IV)
                           )
isa.group.all$IV %<>% round(digits = 3)
isa.group.all$P %<>% round(digits = 3)
isa.group.all$Taxon %<>% str_replace('Nothophaeocryptopus', 'N.')
isa.group.all$Age <- NULL
isa.group.all %>%
  gt() %>%
  tab_header(title = 'Table S3. Results of indicator species analysis for open and closed exposure groups.') %>%
  tab_options(
    heading.title.font.size = px(18),
    heading.align = 'left',
    heading.border.bottom.color = 'white',
    table.border.top.color = 'white',
    table.border.bottom.color = 'white',
    table_body.border.top.color = 'white',
    table_body.border.bottom.color = 'white',
    row_group.border.top.color = 'white',
    row_group.border.bottom.color = 'white',
    column_labels.border.bottom.color = 'black',
    table_body.hlines.color = 'white'
  ) %>%
  tab_style(
    style = cell_text(style = 'italic'),
    locations = cells_body(
      columns = 'Taxon',
      rows = !grepl('sp\\.', Taxon)
      )
    ) %>%
  tab_style(
    style = cell_text(weight = 'bold'),
    locations = cells_column_labels(everything()
    )
  ) %>%
  cols_align(align = 'center') %>%
  cols_width(everything()~px(150)) %>% 
  opt_table_font(font = google_font('Crimson Text')
  ) %>%
  tab_source_note(source_note = 'P values were obtained after 999 permutations and corrected
                    for the false discovery rate using the Benjamini-Hochberg method.') %>%
  gtsave(here(table.out, 'table.s3.png'))

# NMDS ordination of all needle ages ####
all <- nmds.mc.par(phy.in = logged$all, n.cores = 6, trymax = 75, perm = 999)

age.cols <- c('#0172e2', '#fd8a1a', '#dcbaf0', '#b3003b')

all.tree <- nmds.plot(all, note = F, 'NMDS1', 'NMDS2') +
  geom_point(shape = 21, size = 1.5, aes(fill = tree)
             ) +
  scale_fill_colorblind() +
  labs(fill = 'Tree') +
  xlab('') +
  coord_fixed()
all.age <- nmds.plot(all, note = T, 'NMDS1', 'NMDS2') +
  geom_point(shape = 21, size = 1.5, aes(fill = age)
             ) +
  scale_fill_manual(values = age.cols) +
  labs(fill = 'Age') +
  coord_fixed()

all.nmds <- all.tree / all.age + plot_annotation(tag_levels = list(c('A', 'B')
                                                                   )
                                                 ) &
  theme(plot.tag = element_text(size = 10, face = 'bold')
        )
all.nmds
ggsave(here(figure.out, 'fig.s2.png'), all.nmds, units = 'mm', width = 140, height = 140,
       dpi = 300)

# Get session info ####
session.path <- here('output', 'sessions')
dir.create(session.path, recursive = T)
sessionInfo() %>% saveRDS(here(session.path, 'comp.sesh.rds')
                          )
