# Load required libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)

# Define color map
cluster_colors <- c(
  "dI1" = "orange", "dI2" = "gold", "dI3" = "tomato1", "dI4" = "goldenrod3", "dI5" = "olivedrab2",
  "dI6" = "darkolivegreen", "V0" = "lightslateblue", "V1" = "maroon1", "V2a" = "violet",
  "V2b" = "violetred", "V3" = "lightpink1", "MNs" = "green3"
)

# Marker genes
genes <- c(
  "gad2.L", "slc17a7.L", "en1.L", "foxp2.L", "tal1.S", "gata3.L", "gata2.L",
  "wt1.L", "dmrt3.L", "evx1.L", "evx2.S", "lhx3.L", "vsx2.S", "shox2.S",
  "nkx2-2.L", "sim1.S", "slit2.S", "slc18a3.S"
)

# Load processed objects (INs + MNs only)
neu38 <- readRDS("INs_and_MNs_st38.rds")
neu47 <- readRDS("INs_and_MNs_st47.rds")
neu54 <- readRDS("INs_and_MNs_st54.rds")
neu57 <- readRDS("INs_and_MNs_st57.rds")

# UMAP if missing
for (obj in list(neu38, neu47, neu54, neu57)) {
  if (is.null(Embeddings(obj, "umap"))) {
    obj <- RunUMAP(obj, dims = 2:50, verbose = FALSE)
  }
}

# Set levels and cluster identities
for (obj in list(neu38, neu47, neu54, neu57)) {
  levels(obj) <- rev(c("V1", "V2b", "dI6", "V0", "V2a", "V3", "MNs"))
}

# DimPlot and DotPlot for each stage
DimPlot(neu38, label = FALSE, pt.size = 1, cols = cluster_colors) + NoAxes()
DotPlot(neu38, features = genes, cols = c("white", "black")) + RotatedAxis()

# Merge objects with stage labels
neu38$stage <- "38"
neu47$stage <- "47"
neu54$stage <- "54"
neu57$stage <- "57"
neu <- merge(neu38, c(neu47, neu54, neu57))
neu <- JoinLayers(neu)

saveRDS(neu, "INs_MNs_merged_all_stages.rds")

# --- Proportion by cluster and stage ---

make_prop_df <- function(obj, stage) {
  df <- as.data.frame(table(Idents(obj)))
  colnames(df) <- c("cluster", "count")
  df$prop <- df$count / sum(df$count)
  df$stage <- stage
  return(df)
}

df38 <- make_prop_df(neu38, "Stage 38")
df47 <- make_prop_df(neu47, "Stage 47")
df54 <- make_prop_df(neu54, "Stage 54")
df57 <- make_prop_df(neu57, "Stage 57")
df_all <- bind_rows(df38, df47, df54, df57)

ggplot(df_all, aes(x = cluster, y = prop, fill = cluster)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = cluster_colors) +
  facet_wrap(~stage) +
  labs(title = "Cluster Proportions by Stage", x = "Cluster", y = "Proportion") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# --- Proportion of MNs vs INs ---

df_class <- df_all %>%
  mutate(class = ifelse(cluster == "MNs", "MNs", "INs")) %>%
  group_by(stage, class) %>%
  summarise(prop = sum(prop)) %>%
  ungroup()

ggplot(df_class, aes(x = stage, y = prop, fill = class)) +
  geom_col() +
  scale_fill_manual(values = c("MNs" = "green3", "INs" = "dodgerblue")) +
  labs(title = "MNs vs INs Proportions", x = "Stage", y = "Proportion") +
  coord_cartesian(ylim = c(0, 0.45)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

