# 4.1 Define locus sets for supermatrix construction

library(dplyr)
library(readr)

stats_file <- "data.eva_stat_combined.csv"
halfhalf_file <- "halfhalf.txt"
output_dir <- "loci_list"

total_samples <- 72
srh_threshold <- 0.05

dir.create(output_dir, showWarnings = FALSE)

df <- read_csv(stats_file, show_col_types = FALSE)
halfhalf <- readLines(halfhalf_file)

# "halfhalf" identifies loci containing comparable numbers of two orthologs.
df <- df %>%
  mutate(completeness = faa / total_samples)

# Internal dataset labels taxa50/75/90 correspond to
# matrix50/75/90 in the final manuscript.

ALL_taxa50 <- df %>%
  filter(completeness >= 0.50, !(Locus %in% halfhalf)) %>%
  pull(Locus)

ALL_taxa75 <- df %>%
  filter(completeness >= 0.75, !(Locus %in% halfhalf)) %>%
  pull(Locus)

ALL_taxa90 <- df %>%
  filter(completeness >= 0.90, !(Locus %in% halfhalf)) %>%
  pull(Locus)

SRH_taxa50 <- df %>%
  filter(symtest_pval >= srh_threshold,
         completeness >= 0.50,
         !(Locus %in% halfhalf)) %>%
  pull(Locus)

SRH_taxa75 <- df %>%
  filter(symtest_pval >= srh_threshold,
         completeness >= 0.75,
         !(Locus %in% halfhalf)) %>%
  pull(Locus)

SRH_taxa90 <- df %>%
  filter(symtest_pval >= srh_threshold,
         completeness >= 0.90,
         !(Locus %in% halfhalf)) %>%
  pull(Locus)

writeLines(ALL_taxa50, file.path(output_dir, "loci_list_ALL_taxa50.txt"))
writeLines(ALL_taxa75, file.path(output_dir, "loci_list_ALL_taxa75.txt"))
writeLines(ALL_taxa90, file.path(output_dir, "loci_list_ALL_taxa90.txt"))
writeLines(SRH_taxa50, file.path(output_dir, "loci_list_SRH_taxa50.txt"))
writeLines(SRH_taxa75, file.path(output_dir, "loci_list_SRH_taxa75.txt"))
writeLines(SRH_taxa90, file.path(output_dir, "loci_list_SRH_taxa90.txt"))
