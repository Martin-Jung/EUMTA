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
ll <- ll[grep("MTA", ll)]
ss <- ll |> map_dfr(read.csv)
ss$group <- factor(ss$art,
                   levels = c("Art_All", "Art_12","Art_17"),
                   labels = c("All", "Art 12", "Art 17"))

# Quick plot
g1 <- ggplot(ss |> dplyr::filter(ctype=="combined"),
             aes(x = year, y = y,
                    ymin = y - sd,
                    ymax = y + sd,
                    group = option, color = option)) +
  theme_bw(base_size = 20) +
  geom_line(linewidth = 1.5, alpha = .6) +
  # geom_errorbar(position = position_dodge(width = .25)) +
  ggsci::scale_color_d3() +
  scale_y_continuous(limits = c(0,1)) +
  theme(legend.position = "bottom") +
  facet_wrap(variable~group) +
  theme(axis.text.x.bottom = element_text(angle = 45, hjust = 1)) +
  labs(x = "", y = "Mean target achievement (MTA)")
g1

# By scale
g2 <- ggplot(ss, aes(x = year, y = y,
                     ymin = y - sd,
                     ymax = y + sd,
                     group = option, color = option)) +
  theme_bw(base_size = 20) +
  geom_line() +
  # geom_errorbar(position = position_dodge(width = .25)) +
  ggsci::scale_color_d3() +
  scale_y_continuous(limits = c(0,1)) +
  theme(legend.position = "bottom") +
  facet_grid(variable~scale+group) +
  theme(axis.text.x.bottom = element_text(angle = 45, hjust = 1)) +
  labs(x = "", y = "Mean target achievement (MTA)")
g2

# ------------------ #
#### Area increase per year against MTA ####

# Load areas
sc_area <- read.csv("export/Overall_areastatistics.csv")

# Format MTA for overall
mta <- ss |> dplyr::filter(option == "loglinear",
                           category == "All",
                           scale == "MS") |>
  dplyr::select(year,ctype,y)

# Combine both
comb <- dplyr::full_join(sc_area, mta) |> dplyr::filter(ctype != "combined")

g <- ggplot(comb, aes(x = fullprop, y = y ) ) +
  theme_light(base_size = 16) +
  # coord_flip() +
  geom_line(aes(x = fullprop, y = y, group = ctype, color = ctype),
            inherit.aes = FALSE) +
    geom_point(size = 2) +
  scale_x_continuous(limits = c(0,1)) +
  scale_y_continuous(limits = c(0,1)) +
  guides(color = guide_legend(title = "")) +
  geom_abline(slope = 1,intercept = 0, linetype = "dashed", colour = "black") +
  labs(x = "EU land-area conserved (%)", y = "Mean target achievement (MTA)")
g
ggsave(g, filename = "img/Area_with_MTA.png",width = 8, height = 8)
