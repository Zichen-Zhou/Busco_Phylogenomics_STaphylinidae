# Internal labels taxa50/75/90 correspond to matrix50/75/90 in the final manuscript.

# 5.2 RF-distance topology comparison

library(ape)
library(tidyverse)
library(ggrepel)

# Read topology files
tree_files <- list.files(
  "treefiles_all",
  pattern = "\\.(treefile|nwk)$",
  full.names = TRUE
)

# Parse method and dataset metadata
parse_metadata <- function(filepath) {
  
  fname <- basename(filepath)
  
  is_4group <- fname == "PMSF_FLT_taxa50_4group.treefile"
  is_sse_s <- fname == "PMSF_FLT_taxa50_SSE-S.treefile"
  is_constrain <- is_4group | is_sse_s
  
  model <- str_extract(fname, "LG|GHOST|EXEHO|PMSF|Dayhoff6|ASTRAL")
  
  dataset <- str_extract(
    fname,
    "ALL_taxa90|ALL_taxa75|ALL_taxa50|SRH_taxa90|SRH_taxa75|SRH_taxa50|FLT_taxa75|FLT_taxa50"
  )
  
  is_len70 <- str_detect(fname, "_len70")
  label_core <- paste(model, dataset, sep = "_")
  
  label_plot <- case_when(
    is_4group ~ "PMSF_FLT_taxa50_4group",
    is_sse_s ~ "PMSF_FLT_taxa50_SSE-S",
    is_len70 ~ paste0(label_core, "_len70"),
    TRUE ~ label_core
  )
  
  tibble(
    file = filepath,
    label_id = fname,
    label_plot = label_plot,
    model = model,
    dataset = dataset,
    is_len70 = is_len70,
    is_4group = is_4group,
    is_sse_s = is_sse_s,
    is_constrain = is_constrain
  )
}

metadata <- purrr::map_dfr(tree_files, parse_metadata)

# Read trees
trees <- purrr::map(metadata$file, read.tree)

# Restrict all trees to the shared tip set
common_tips <- purrr::reduce(
  purrr::map(trees, ~ .x$tip.label),
  intersect
)

pruned_trees <- purrr::map(
  trees,
  ~ drop.tip(.x, setdiff(.x$tip.label, common_tips))
)

# Unroot trees
unrooted_trees <- purrr::map(pruned_trees, unroot)

# Compute pairwise Robinson-Foulds distances
n <- length(unrooted_trees)

rf_matrix <- matrix(
  0,
  nrow = n,
  ncol = n,
  dimnames = list(metadata$label_id, metadata$label_id)
)

for (i in 1:(n - 1)) {
  for (j in (i + 1):n) {
    rf <- dist.topo(unrooted_trees[[i]], unrooted_trees[[j]])
    rf_matrix[i, j] <- rf
    rf_matrix[j, i] <- rf
  }
}

# Perform classical multidimensional scaling
mds_res <- cmdscale(as.dist(rf_matrix), k = 2)

mds_df <- as.data.frame(mds_res) %>%
  tibble::rownames_to_column("label_id") %>%
  left_join(metadata, by = "label_id")

# Plot settings
model_shapes <- c(
  "LG" = 21,
  "GHOST" = 22,
  "EXEHO" = 23,
  "PMSF" = 24,
  "Dayhoff6" = 25,
  "ASTRAL" = 8
)

dataset_levels <- c(
  "ALL_taxa50",
  "SRH_taxa50",
  "FLT_taxa50",
  "ALL_taxa75",
  "SRH_taxa75",
  "FLT_taxa75",
  "ALL_taxa90",
  "SRH_taxa90"
)

dataset_colors <- c(
  "ALL_taxa50" = "#ffafaf",
  "ALL_taxa75" = "#e3c378",
  "ALL_taxa90" = "#8499ce",
  "SRH_taxa50" = "#f76060",
  "SRH_taxa75" = "#f9c261",
  "SRH_taxa90" = "#535497",
  "FLT_taxa50" = "#a71200",
  "FLT_taxa75" = "#a77b40"
)

base_size <- 5

# Define plotting groups
mds_df <- mds_df %>%
  mutate(
    dataset = factor(dataset, levels = dataset_levels),
    style_flag = ifelse(
      is_len70 | grepl("^FLT", dataset),
      "len70",
      "len0"
    ),
    coverage = factor(
      ifelse(style_flag == "len0", "all sequences", "length≥70%"),
      levels = c("all sequences", "length≥70%")
    )
  )

# Split data by display style
df_len0 <- mds_df %>%
  filter(style_flag == "len0", !is_constrain)

df_len70_other <- mds_df %>%
  filter(style_flag == "len70", model != "ASTRAL", !is_constrain)

df_len70_astral <- mds_df %>%
  filter(style_flag == "len70", model == "ASTRAL", !is_constrain)

df_4group <- mds_df %>%
  filter(is_4group)

df_sse_s <- mds_df %>%
  filter(is_sse_s)

df_labels <- mds_df

# Coverage legend
df_cov_legend <- data.frame(
  V1 = rep(mds_df$V1[1], 2),
  V2 = rep(mds_df$V2[1], 2),
  coverage = factor(
    c("all sequences", "length≥70%"),
    levels = c("all sequences", "length≥70%")
  )
)

# Plot with labels
p_label <- ggplot() +
  
  geom_point(
    data = df_cov_legend,
    aes(x = V1, y = V2, alpha = coverage),
    inherit.aes = FALSE,
    size = 0,
    color = "black",
    show.legend = TRUE
  ) +
  
  geom_point(
    data = df_len0,
    aes(V1, V2, shape = model, color = dataset, fill = dataset),
    size = base_size,
    stroke = 0.5
  ) +
  
  geom_point(
    data = df_len70_other,
    aes(V1, V2, shape = model, color = dataset),
    fill = NA,
    size = base_size,
    stroke = 0.5
  ) +
  
  geom_point(
    data = df_len70_astral,
    aes(V1, V2, shape = model, color = dataset),
    fill = NA,
    size = base_size * 0.7,
    stroke = 0.5
  ) +
  
  geom_point(
    data = df_4group,
    aes(V1, V2),
    shape = 11,
    fill = NA,
    color = "black",
    size = base_size,
    stroke = 1.2
  ) +
  
  geom_point(
    data = df_sse_s,
    aes(V1, V2),
    shape = 11,
    fill = NA,
    color = "#A9D08E",
    size = base_size,
    stroke = 1.2
  ) +
  
  geom_text_repel(
    data = df_labels,
    aes(V1, V2, label = label_plot, color = dataset),
    box.padding = 0.25,
    point.padding = 0.2,
    min.segment.length = 0,
    max.overlaps = Inf,
    force = 1,
    size = 2.5,
    segment.color = "gray40",
    direction = "y",
    show.legend = FALSE
  ) +
  
  scale_shape_manual(
    values = model_shapes,
    name = "Method"
  ) +
  
  scale_color_manual(
    values = dataset_colors,
    name = "Dataset",
    drop = FALSE
  ) +
  
  scale_fill_manual(
    values = dataset_colors,
    guide = "none",
    drop = FALSE
  ) +
  
  scale_alpha_manual(
    values = c("all sequences" = 1, "length≥70%" = 1),
    name = "Length Coverage",
    breaks = c("all sequences", "length≥70%")
  ) +
  
  theme_minimal() +
  
  coord_fixed() +
  
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  ) +
  
  labs(
    x = "MDS1",
    y = "MDS2",
    title = "MDS of RF Distances"
  ) +
  
  scale_x_continuous(
    breaks = c(-15, -10, -5, 0, 5, 10, 15),
    minor_breaks = NULL
  ) +
  
  scale_y_continuous(
    breaks = c(-15, -10, -5, 0, 5, 10, 15),
    limits = c(-15, 15),
    minor_breaks = NULL
  ) +
  
  guides(
    shape = guide_legend(order = 1),
    color = guide_legend(order = 2),
    alpha = guide_legend(
      order = 3,
      override.aes = list(
        shape = c(16, 1),
        size = c(3, 3),
        color = "black",
        fill = NA,
        alpha = 1,
        stroke = c(0.8, 1.2)
      )
    )
  )

p_label

# Plot with jitter
p_jitter <- ggplot() +
  
  geom_point(
    data = df_cov_legend,
    aes(x = V1, y = V2, alpha = coverage),
    inherit.aes = FALSE,
    size = 0,
    color = "black",
    show.legend = TRUE
  ) +
  
  geom_jitter(
    data = df_len0,
    aes(V1, V2, shape = model, color = dataset, fill = dataset),
    size = base_size,
    stroke = 0.9,
    width = 0.08,
    height = 0.08
  ) +
  
  geom_jitter(
    data = df_len70_other,
    aes(V1, V2, shape = model, color = dataset),
    fill = NA,
    size = base_size,
    stroke = 1.0,
    width = 0.08,
    height = 0.08
  ) +
  
  geom_jitter(
    data = df_len70_astral,
    aes(V1, V2, shape = model, color = dataset),
    fill = NA,
    size = base_size * 0.7,
    stroke = 1.0,
    width = 0.08,
    height = 0.08
  ) +
  
  geom_point(
    data = df_4group,
    aes(V1, V2),
    shape = 11,
    fill = NA,
    color = "black",
    size = base_size,
    stroke = 1.2
  ) +
  
  geom_point(
    data = df_sse_s,
    aes(V1, V2),
    shape = 11,
    fill = NA,
    color = "#A9D08E",
    size = base_size,
    stroke = 1.2
  ) +
  
  scale_shape_manual(
    values = model_shapes,
    name = "Method"
  ) +
  
  scale_color_manual(
    values = dataset_colors,
    name = "Dataset",
    drop = FALSE
  ) +
  
  scale_fill_manual(
    values = dataset_colors,
    guide = "none",
    drop = FALSE
  ) +
  
  scale_alpha_manual(
    values = c("all sequences" = 1, "length≥70%" = 1),
    name = "Length Coverage",
    breaks = c("all sequences", "length≥70%")
  ) +
  
  theme_minimal() +
  
  coord_fixed() +
  
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  ) +
  
  labs(
    x = "MDS1",
    y = "MDS2",
    title = "MDS of RF Distances"
  ) +
  
  scale_x_continuous(
    breaks = c(-15, -10, -5, 0, 5, 10, 15),
    minor_breaks = NULL
  ) +
  
  scale_y_continuous(
    breaks = c(-15, -10, -5, 0, 5, 10, 15),
    limits = c(-15, 15),
    minor_breaks = NULL
  ) +
  
  guides(
    shape = guide_legend(
      order = 1,
      override.aes = list(size = 3)
    ),
    color = guide_legend(
      order = 2,
      override.aes = list(size = 3)
    ),
    alpha = guide_legend(
      order = 3,
      override.aes = list(
        shape = c(16, 1),
        size = c(3, 3),
        color = "black",
        fill = NA,
        alpha = 1,
        stroke = c(0.8, 1.2)
      )
    )
  )

p_jitter


ggsave(
  "RF_MDS_label.pdf",
  p_label,
  width = 12,
  height = 8,
  dpi = 600
)

ggsave(
  "RF_MDS_jitter.pdf",
  p_jitter,
  width = 12,
  height = 8,
  dpi = 600
)
