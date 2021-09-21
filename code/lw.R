# Create Fig. S2 from the manuscript, which shows:
# (A) mean leaf wetness vs height,
# (B) seasonality of leaf wetness vs height,
# (C) mean daily leaf wetness range vs height
# (D) leaf wetness continentality vs height

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
lw <- readRDS(here('data', 'compile', 'discovery.lw.rds')) %>% 
  filter(year == 2017)

# Match tree colors ####
tree.cols <- '#56B4E9'
line.size <- 1.5
point.size <- 3

# Mean daily leaf wetness range ####
diff <- lw %>%
  group_by(tree, height, ymd) %>%
  summarize(max.lw = max(lw, na.rm = T),
            min.lw = min(lw, na.rm = T),
            diff.lw = max.lw - min.lw) %>%
  group_by(tree, height) %>% 
  summarize(mean.diff = mean(diff.lw))

diff.plot <- ggplot(diff, aes(x = height, y = mean.diff, color = tree)) +
  geom_line(size = line.size) +
  geom_point(size = point.size) +
  scale_color_manual(values = tree.cols) +
  ylab(expression(atop(bold(Mean~daily), bold(leaf~wetness~range~(mV))))) +
  xlab("Height (m)") +
  labs(color = 'Tree') +
  theme_cowplot() +
  theme(axis.title.x = element_text(face = 'bold', size = 7),
        axis.text.x = element_text(size = 7),
        axis.title.y = element_text(face = 'bold', size = 7),
        axis.text.y = element_text(size = 7),
        legend.position = 'none')
diff.plot

# Mean leaf wetness and seasonality (standard deviation of leaf wetness) ####
mean.sd.lw <- lw %>%
  group_by(tree, height) %>% 
  summarize(mean = mean(lw, na.rm = T),
            sd = sd(lw, na.rm = T))

# Mean leaf wetness (all time points averaged) ####
mean.plot <- ggplot(mean.sd.lw, aes(x = height, y = mean, color = tree)) +
  geom_line(size = line.size) +
  geom_point(size = point.size) +
  scale_color_manual(values = tree.cols) +
  ylab(expression(atop(bold(Mean~leaf~wetness~(mV)), ''))) +
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
sd.plot <- ggplot(mean.sd.lw, aes(x = height, y = sd, color = tree)) +
  geom_line(size = line.size) +
  geom_point(size = point.size) +
  scale_color_manual(values = tree.cols) +
  ylab(expression(atop(bold(Seasonality), bold((paste('standard deviation, ', mV)))))) +
  xlab("\nHeight (m)") +
  labs(color = 'Tree') +
  theme_cowplot() +
  theme(axis.title.x = element_text(face = 'bold', size = 7),
        axis.text.x = element_text(size = 7),
        axis.title.y = element_text(face = 'bold', size = 7),
        axis.text.y = element_text(size = 7),
        legend.position = 'none')
sd.plot

# Continentality (Mean of wettest month, subtracting the mean of the driest month) ####
mar.aug <- lw %>%
  group_by(tree, height, month) %>% 
  summarize(mean = mean(lw, na.rm = T)) %>% 
  filter(month == 3 | month == 8)

mar <- mar.aug %>% filter(month == 3)
aug <- mar.aug %>% filter(month == 8)

mar$cont <- mar$mean - aug$mean

cont.plot <- ggplot(mar, aes(x = height, y = mean, color = tree)) +
  geom_line(size = line.size) +
  geom_point(size = point.size) +
  scale_color_manual(values = tree.cols) +
  ylab(expression(atop(bold(Continentality), bold(('Mar. - Aug. mean leaf wetness, mV'))))) +
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
lw.plots <- (mean.plot | sd.plot) / (diff.plot | cont.plot) +
  plot_annotation(tag_levels = list(c('A', 'B', 'C', 'D'))) &
  theme(plot.tag = element_text(size = 10, face = 'bold'))
lw.plots
ggsave(here(figure.out, 'fig.s2.tiff'), units = 'mm', height = 140, width = 190,
       dpi = 300, compression = 'lzw')

# Get session info ####
session.path <- here('output', 'sessions')
dir.create(session.path, recursive = T)
sessionInfo() %>% saveRDS(here(session.path, 'lw.sesh.rds'))
