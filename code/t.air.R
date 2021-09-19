# Create Figure S1 from the manuscript, which shows:
# (a) mean air temperature vs height,
# (b) seasonality of air temperature vs height,
# (c) mean daily air temperature range vs height
# (d) air temperature continentality vs height

library(patchwork)
library(ggthemes)
library(cowplot)
library(tidyverse)
library(magrittr)
library(here)

# Consult "code/functions.R" for descriptions of function applications, inputs, and outputs
source(here('code', 'functions.R'))

# Set out-paths ####
figure.out <- here('output', 'figures')
dir.create(figure.out, recursive = T)

# Read in files ####
t.air <- readRDS(here('data', 'compile', 'filled.pheno.t.air.rds')) %>%
  filter(year != 2011 & year != 2018) %>%
  filter(tree != 'PC17'| height != 21.5) # The logger for this height and tree appears to have malfunctioned.

# Match tree colors ####
tree.cols <- c('#000000', '#0072B2', '#D55E00', '#CC79A7')
line.size <- 1.5
point.size <- 3

# Calculate statistics ####
# Mean daily range
diff <- t.air %>%
  group_by(tree, height, ymd) %>%
  summarize(max.t.air = max(t.air, na.rm = T),
            min.t.air = min(t.air, na.rm = T),
            diff.t.air = max.t.air - min.t.air) %>%
  group_by(tree, height) %>% 
  summarize(mean.diff = mean(diff.t.air))

diff.plot <- ggplot(diff, aes(x = height, y = mean.diff, color = tree)) +
  geom_line(size = line.size) +
  geom_point(size = point.size) +
  scale_color_manual(values = tree.cols) +
  ylab(expression(bold(Mean~daily~air~temp.~range~(paste(degree, C))))) +
  xlab("Height (m)") +
  labs(color = 'Tree') +
  theme_cowplot() +
  theme(axis.title.x = element_text(face = 'bold', size = 7),
        axis.text.x = element_text(size = 7),
        axis.title.y = element_text(face = 'bold', size = 7),
        axis.text.y = element_text(size = 7),
        legend.position = 'none')
diff.plot

# Mean air temperature and seasonality (standard deviation of air temperature) ####
mean.sd.t.air <- t.air %>%
  group_by(tree, height) %>% 
  summarize(mean = mean(t.air, na.rm = T),
            sd = sd(t.air, na.rm = T))

# Mean air temperature (all time points averaged)
mean.plot <- ggplot(mean.sd.t.air, aes(x = height, y = mean, color = tree)) +
  geom_line(size = line.size) +
  geom_point(size = point.size) +
  scale_color_manual(values = tree.cols) +
  ylab(expression(atop(bold(Mean~air~temp.~(paste(degree, C))), ''))) +
  xlab("\nHeight (m)") +
  labs(color = 'Tree') +
  theme_cowplot() +
  theme(axis.title.x = element_text(face = 'bold', size = 7),
        axis.text.x = element_text(size = 7),
        axis.title.y = element_text(face = 'bold', size = 7),
        axis.text.y = element_text(size = 7),
        legend.position = 'none')
mean.plot

# Seasonality sensu Oita et al. 2021 (i.e., the standard deviation)
sd.plot <- ggplot(mean.sd.t.air, aes(x = height, y = sd, color = tree)) +
  geom_line(size = line.size) +
  geom_point(size = point.size) +
  scale_color_manual(values = tree.cols) +
  # ylab(expression(atop(bold(Seasonality~(paste('SD, ', degree, C))), ''))) +
  ylab(expression(atop(bold(Seasonality), bold((paste('standard deviation, ', degree, C)))))) +
  xlab("\nHeight (m)") +
  labs(color = 'Tree') +
  theme_cowplot() +
  theme(axis.title.x = element_text(face = 'bold', size = 7),
        axis.text.x = element_text(size = 7),
        axis.title.y = element_text(face = 'bold', size = 7),
        axis.text.y = element_text(size = 7),
        legend.position = 'none')
sd.plot

# Continentality (Mean of warmest month, subtracting the mean of the coldest month) ####
july.dec <- t.air %>%
  group_by(tree, height, month) %>% 
  summarize(mean = mean(t.air, na.rm = T)) %>% 
  filter(month == 7 | month == 12)

july <- july.dec %>% filter(month == 7)
dec <- july.dec %>% filter(month == 12)

july$cont <- july$mean - dec$mean

cont.plot <- ggplot(july, aes(x = height, y = cont, color = tree)) +
  geom_line(size = line.size) +
  geom_point(size = point.size) +
  scale_color_manual(values = tree.cols) +
  ylab(expression(atop(bold(Continentality), bold((paste('July - Dec. mean air temp., ', degree, C)))))) +
  xlab("\nHeight (m)") +
  labs(color = 'Tree') +
  theme_cowplot() +
  theme(axis.title.x = element_text(face = 'bold', size = 7),
        axis.text.x = element_text(size = 7),
        axis.title.y = element_text(face = 'bold', size = 7),
        axis.text.y = element_text(size = 7),
        legend.title = element_text(face = 'bold', size = 7),
        legend.text = element_text(size = 7))
cont.plot

# Combine all plots and output ####
t.air.plots <- (mean.plot | sd.plot) / (diff.plot | cont.plot) +
  plot_annotation(tag_levels = list(c('A', 'B', 'C', 'D'))) &
  theme(plot.tag = element_text(size = 10, face = 'bold'))
t.air.plots
ggsave(here(figure.out, 'fig.s1.tiff'), units = 'mm', width = 190, height = 140,
       dpi = 500, compression = 'lzw')

# Get session info ####
session.path <- here('output', 'sessions')
dir.create(session.path, recursive = T)
sessionInfo() %>% saveRDS(here(session.path, 't.air.sesh.rds'))
