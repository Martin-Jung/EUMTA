# # # # # # # # # # # # # # # #
# This script makes some summary statistics and visualizations
#
# Author:
# Martin Jung (IIASA) - 2024
# # # # # # # # # # # # # # # #

library(ggplot2)
library(ggsci)
library(scales)
library(tidyverse)
library(assertthat)
library(purrr)

# Path to the outputs
path_output <- "export/"
dir.create(path_output, showWarnings = FALSE)

# ------------- #
# Load and combine all csv files
ll <- list.files(path_output, full.names = TRUE)
ll <- ll[has_extension(ll, "csv")]
ss <- ll |> map_dfr(read.csv)
ss$group <- factor(ss$art,
                   levels = c("Art_All", "Art_12","Art_17"),
                   labels = c("All", "Art 12", "Art 17"))

# Quick plot# Quick plot# Quick plot
g <- ggplot(ss, aes(x = year, y = y,
                    ymin = y - sd,
                    ymax = y + sd,
                    group = option, color = option)) +
  theme_bw(base_size = 20) +
  geom_line(linewidth = 1.5, alpha = .6) +
  # geom_errorbar(position = position_dodge(width = .25)) +
  ggsci::scale_color_d3() +
  theme(legend.position = "bottom") +
  facet_wrap(variable~group) +
  theme(axis.text.x.bottom = element_text(angle = 45, hjust = 1)) +
  labs(x = "", y = "Mean target achievement (MTA)")
g
