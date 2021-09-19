# Creates Table S1, which: ####
# reports the results of ANOVAs tests for whether diversity differs among tree, needle age class, and exposure group,
# differ more than expected by chance.

# Creates Figure 2, which: ####
# (a) shows how OTU richness varies among needle age classes,
# (b) models the relationship between estimated richness and closure 
# (c) how the Shannon diversity index differs among exposure groups, and
# (d) how the Shannon diversity index differs among trees.

# Creates Figure S4, which: ####
# models the relationship between the estimated Shannon index and closure

# New errors ####

# Warning message:
#   In is.na(x) : is.na() applied to non-(list or vector) of type 'expression'
# It's tied to the labels for the lme plots. Low priority, but does need addressing...

# Load packages ####
library(patchwork)
library(cowplot)
library(ggthemes)
library(gt)
library(iNEXT)
library(phyloseq)
library(tidyverse)
library(piecewiseSEM)
library(nlme)
library(car)
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
perf.n <- readRDS(here('data', 'compile', 'perf.n.clean.rds'))

# Image no longer needed
# ages <- png::readPNG(here('data', 'media', 'shaw.age.class.png'), native = T)

# Estimate richness and diversity ####
perf.n.div <- inext.div(perf.n$counts$all, cores = 6)

rich <- perf.n.div %>% filter(metric == 'q0')
shannon <- perf.n.div %>% filter(metric == 'q1')

# ANOVA ####
# Tree and age ####
perf.n.aov.tree.age <- perf.n.div %>%
  select(sample, age, tree, group, metric, estimate) %>%
  group_by(metric) %>%
  group_map(~ anova.crown(.x, form = 'estimate ~ age*tree'), .keep = T) %>%
  bind_rows()

# Group and age ####
perf.n.aov.group.age <- perf.n.div %>%
  select(sample, age, tree, group, metric, estimate) %>%
  group_by(metric) %>%
  group_map(~ anova.crown(.x, form = 'estimate ~ age*group'), .keep = T) %>%
  bind_rows()

rbind(perf.n.aov.tree.age, perf.n.aov.group.age) %>%
  gt(groupname_col = c('metric')) %>%
  tab_header(title = 'Table S1. ANOVA results testing variation in alpha diversity in and among trees, exposure groups, and needle ages.') %>%
  cols_align(align = 'center') %>%
  tab_style(
    style = cell_text(style = 'italic'),
    locations = cells_body(
      columns = vars(Term),
      rows = grepl('^[[:lower:]]', Term))) %>%
  fmt_missing(columns = everything(), missing_text = '') %>%
  tab_source_note(source_note = 'Estimates of OTU richness and the Shannon index of diversity were interpolated or extrapolated to a depth of 1000 reads.') %>%
  opt_table_font(font = google_font('Crimson Text')) %>%
  gtsave(filename = here(table.out, 'table.s1.png'))

# Tukey ####
# Age ####
perf.n.tukey.age <- perf.n.div %>%
  select(sample, age, tree, group, metric, estimate) %>%
  group_by(metric) %>%
  group_map(~ tukey(.x, form = 'age'), .keep = T) %>%
  bind_rows()
perf.n.tukey.age$estimate <- ifelse(perf.n.tukey.age$metric == 'q0', 22, 8.5)

# Tree ####
perf.n.tukey.tree <- perf.n.div %>%
  select(sample, age, tree, group, metric, estimate) %>%
  group_by(metric) %>%
  group_map(~ tukey(.x, form = 'tree'), .keep = T) %>%
  bind_rows()
perf.n.tukey.tree$estimate <- ifelse(perf.n.tukey.tree$metric == 'q0', 21.5, 3.1)

# Bootstrapped means and confidence intervals ####
# Age ####
perf.n.boot.age <- perf.n.div %>%
  select(age, tree, group, metric, estimate) %>%
  group_by(metric, age) %>%
  group_map(~ mean.boot(.x, .y, var = 'age')) %>%
  bind_rows %>%
  group_by(metric, factor) %>%
  summarize(mean = mean(estimate),
            lci = coxed::bca(estimate)[1],
            uci = coxed::bca(estimate)[2])
perf.n.boot.age$age <- perf.n.boot.age$factor$age
perf.n.boot.age$factor <- NULL

# Tree ####
perf.n.boot.tree <- perf.n.div %>%
  select(age, tree, group, metric, estimate) %>%
  group_by(metric, tree) %>%
  group_map(~ mean.boot(.x, .y, var = 'tree')) %>%
  bind_rows %>%
  group_by(metric, factor) %>%
  summarize(mean = mean(estimate),
            lci = coxed::bca(estimate)[1],
            uci = coxed::bca(estimate)[2])
perf.n.boot.tree$tree <- perf.n.boot.tree$factor$tree
perf.n.boot.tree$factor <- NULL

# Group ####
perf.n.boot.group <- perf.n.div %>%
  select(age, tree, group, metric, estimate) %>%
  group_by(metric, group) %>%
  group_map(~ mean.boot(.x, .y, var = 'group')) %>%
  bind_rows %>%
  group_by(metric, factor) %>%
  summarize(mean = mean(estimate),
            lci = coxed::bca(estimate)[1],
            uci = coxed::bca(estimate)[2])
perf.n.boot.group$group <- perf.n.boot.group$factor$group
perf.n.boot.group$factor <- NULL

# Linear mixed-effect modeling ####
# Richness vs crown variables ####
# Full model
rich.lme.full <- lme(estimate ~ age + height + depth + closure,
                     random = ~ 1|tree,
                     data = rich,
                     method = 'ML')

# Reduced - depth and height removed. Closure + age was the most parsimonious model.
rich.lme.red <- lme(estimate ~ age + closure,
                    random = ~ 1|tree,
                    data = rich,
                    method = 'ML')

# Random intercept only (for comparing overall model significance)
rich.lme.base <- lme(estimate ~ 1,
                     random = ~1|tree,
                     data = rich,
                     method = 'ML')
# Check residuals
plot(rich.lme.red)
# Overview
summary(rich.lme.red)
# Type II ANOVA
car::Anova(rich.lme.red)
# Compare full and reduced
anova(rich.lme.full, rich.lme.red)
# Compare reduced and base model to obtain overall model significance
rich.lme.p <- anova(rich.lme.red, rich.lme.base)$`p-value`[[2]] %>% round(digits = 3)
rich.lme.p[rich.lme.p < 0.001] <- '< 0.001'
# Obtain an R2 for the model
rich.lme.r2 <- rsquared(rich.lme.red)$Marginal %>% round(digits = 3)

# Diversity vs crown variables ####
# Full model
shannon.lme.full <- lme(log2(estimate) ~ height + depth + closure,
                        random = ~1|tree,
                        data = shannon,
                        method = 'ML')

# Reduced model - any crown variable suffices, but closure left for simplicity
# Age does not contribute to interpretation.
shannon.lme.red <- lme(log2(estimate) ~ closure,
                       random = ~1|tree,
                       data = shannon,
                       method = 'ML')
# Random intercept only (for comparing overall model significance)
shannon.lme.base <- lme(log2(estimate) ~ 1,
                        random = ~1|tree,
                        data = shannon,
                        method = 'ML')
# Check residuals
plot(shannon.lme.red)
# Overview
summary(shannon.lme.red)
# Compare full and reduced
anova(shannon.lme.full, shannon.lme.red)
# Compare reduced and base model to obtain overall model significance
shannon.lme.p <- anova(shannon.lme.red, shannon.lme.base)$`p-value`[[2]] %>% round(digits = 3)
shannon.lme.p[shannon.lme.p < 0.001] <- '< 0.001'
# Obtain an R2 for the model
shannon.lme.r2 <- rsquared(shannon.lme.red)$Marginal %>% round(digits = 3)

# Plots ####
# Richness, age, pairwise ####
age.q0.pair <- ggplot(filter(perf.n.div, metric == 'q0'),
                      aes(x = age, y = estimate)) +
  geom_violin(fill = 'white', color = 'black', draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_text(aes(label = Letters),
            filter(perf.n.tukey.age, metric == 'q0'),
            fontface = 'bold', size = 2.5) +
  geom_pointrange(aes(y = mean, ymin = lci, ymax = uci),
                  filter(perf.n.boot.age, metric == 'q0'),
                  size = 0.25, color = 'red') +
  xlab('') +
  ylab('Estimated richness\n') +
  labs(fill = 'Age') +
  scale_fill_colorblind() +
  theme_cowplot() +
  theme(axis.text.x = element_text(face = 'bold', size = 7),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(face = 'bold', size = 7),
        axis.text.y = element_text(size = 7))

# Richness vs closure plot ####
rich$lme.fixed <- rich.lme.red$fitted %>% data.frame() %>% .$fixed
rich.lab <- substitute(P~p*','~R^2~"="~r2, list(p = rich.lme.p, r2 = rich.lme.r2)) %>% as.expression()
age.cols <- c('#ff0000', '#ff0099', '#ff00ff', '#0000ff')

rich.closure.lme <- ggplot(rich, aes(x = closure, y = estimate)) +
  geom_point(size = 1.5, shape = 21, aes(fill = tree), show.legend = F) +
  geom_line(aes(y = lme.fixed, color = age), size = 1) +
  xlab('\nCrown closure index') +
  ylab('Estimated richness\n') +
  labs(color = 'Age', fill = 'Tree') +
  annotate('text', x = 10, y = 20, label = rich.lab, parse = T, size = 2.5) +
  scale_fill_colorblind() +
  scale_color_manual(values = age.cols) +
  theme_cowplot() +
  theme(legend.title = element_text(size = 7, face = 'bold'),
        legend.text = element_text(size = 7),
        axis.title.x = element_text(face = 'bold', size = 7),
        axis.text.x = element_text(size = 7),
        axis.title.y = element_text(face = 'bold', size = 7),
        axis.text.y = element_text(size = 7))

# Shannon, tree, pairwise ####
tree.q1.pair <- ggplot(filter(perf.n.div, metric == 'q1'),
                       aes(x = tree, y = log2(estimate))) +
  geom_violin(aes(fill = tree),
              color = 'white', draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_text(aes(label = Letters, y = estimate),
            filter(perf.n.tukey.tree, metric == 'q1'),
            fontface = 'bold', size = 2.5) +
  geom_pointrange(aes(y = log2(mean), ymin = log2(lci), ymax = log2(uci)),
                  filter(perf.n.boot.tree, metric == 'q1'),
                  size = 0.25, color = 'white') +
  xlab('') +
  ylab('Estimated Shannon index\n') +
  labs(fill = 'Tree') +
  scale_fill_colorblind() +
  theme_cowplot() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(face = 'bold', size = 7),
        axis.text.y = element_text(size = 7),
        legend.title = element_text(face = 'bold', size = 7),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.75, 'line'))

# Shannon, group ####
group.q1.pair <- ggplot(filter(perf.n.div, metric == 'q1')) +
  geom_violin(aes(x = group, y = log2(estimate)),
              fill = 'white', color = 'black',
              draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_pointrange(aes(x = group, y = log2(mean),
                      ymin = log2(lci), ymax = log2(uci)),
                  filter(perf.n.boot.group, metric == 'q1'),
                  size = 0.25, color = 'red') +
  # annotate('text', x = 'Closed', y = 2.5,
  #          label = '*', hjust = -6, size = 7.5, fontface = 2) +
  xlab('') +
  ylab('Estimated Shannon index\n') +
  labs(fill = 'Group') +
  scale_fill_colorblind() +
  theme_cowplot() +
  theme(axis.text.x = element_text(face = 'bold', size = 7),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(face = 'bold', size = 7),
        axis.text.y = element_text(size = 7))

# Figure output
(age.q0.pair + rich.closure.lme) / (tree.q1.pair + group.q1.pair)  +
  plot_annotation(tag_levels = list(c('A', 'B', 'C', 'D'))) &
  theme(plot.tag = element_text(size = 10, face = 'bold'))
ggsave(here(figure.out, 'fig.3.tiff'), units = 'mm', width = 190, height = 120,
       dpi = 500, compression = 'lzw')

# Diversity vs closure plot ####
shannon$lme.fixed <- shannon.lme.red$fitted %>% data.frame() %>% .$fixed
shannon.lab <- substitute(P~p*','~R^2~"="~r2, list(p = shannon.lme.p, r2 = shannon.lme.r2)) %>% as.expression()

shannon.closure.lme.plot <- ggplot(shannon, aes(x = closure, y = estimate)) +
  geom_point(size = 1.5, shape = 21, aes(fill = tree)) +
  geom_line(aes(y = lme.fixed), size = 1.5) +
  xlab('\nCrown closure index') +
  ylab('Estimated Shannon index\n') +
  labs(color = 'Age', fill = 'Tree') +
  annotate('text', x = 10, y = 8, label = shannon.lab, parse = T, size = 2.5) +
  scale_fill_colorblind() +
  theme_cowplot() +
  theme(legend.title = element_text(size = 7, face = 'bold'),
        legend.text = element_text(size = 7),
        axis.title.x = element_text(face = 'bold', size = 7),
        axis.text.x = element_text(size = 7),
        axis.title.y = element_text(face = 'bold', size = 7),
        axis.text.y = element_text(size = 7))

# Figure output
# rich.closure.lme.plot / shannon.closure.lme.plot +
#   plot_annotation(tag_levels = list(c('(a)', '(b)'))) &
#   theme(plot.tag = element_text(size = 10, face = 'bold'))
shannon.closure.lme.plot
ggsave(here(figure.out, 'fig.s4.tiff'), units = 'mm', width = 140,
       dpi = 500, compression = 'lzw')

# Get session info ####
session.path <- here('output', 'sessions')
dir.create(session.path, recursive = T)
sessionInfo() %>% saveRDS(here(session.path, 'div.sesh.rds'))
