suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(purrr))

neurons <- readRDS("out/neurons.rds")
markers <- read.csv("out/markers.csv")

ggsave(
  "out/figure_4d.png",
  DimPlot(
    neurons,
    label = TRUE, pt.size = 1
  ) + FeaturePlot(
    neurons,
    features = c(
      "Elavl4",
      "Tubb3",
      "Uchl1",
      "Chat",
      "Nos1",
      "percent.td"
    ),
    cols = c("lightgrey", "red")
  ),
  units = "px",
  width = 4096,
  height = 2160
)

ggsave(
  "out/figure_6a.png",
  DimPlot(
    neurons,
    label = TRUE, pt.size = 1
  ) + (DoHeatmap(
    neurons,
    features = subset(markers, p_val_adj < 0.05)[["gene"]],
  ) + theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )),
  units = "px",
  width = 4096,
  height = 2160
)

ggsave(
  "out/figure_6b.png",
  FeaturePlot(
    neurons,
    features = c("Grp", "Calcb", "Sst"),
    cols = c("lightgrey", "blue")
  ),
  units = "px",
  width = 4096,
  height = 2160
)
