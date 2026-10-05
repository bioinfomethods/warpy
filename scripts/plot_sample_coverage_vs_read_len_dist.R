#!/usr/bin/env Rscript
# Plot cumulative coverage vs sequence length threshold, per sample, with
# reference lines for: read length N50, and the length threshold at which
# coverage drops to 20X and 50X (if those levels are reached).
#
# For each sample, at every observed read length x, this computes the total
# bases (bp) contributed by all reads with length >= x:
#     bases(x) = sum(length * count) for length >= x
# This lets you answer: "how much coverage/yield comes from reads >= x bp?"
#
# If genome_size_bp is provided, bases are converted to coverage depth (X)
# by dividing by genome size; otherwise the y-axis is raw total bp (in which
# case the 20X/50X reference lines are skipped, since they require depth).
#
# Usage:
#   Rscript plot_sample_coverage_vs_read_len_dist.R <input_tsv> <output_pdf> <summary_tsv> [genome_size_bp]
#
# Example:
#   Rscript plot_sample_coverage_vs_read_len_dist.R read_length_dist.tsv depth_vs_read_len.pdf depth_vs_read_len.txt 3100000000

suppressMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
input_tsv   <- ifelse(length(args) >= 1, args[1], "")
output_pdf  <- ifelse(length(args) >= 2, args[2], "depth_vs_read_len.pdf")
summary_tsv <- ifelse(length(args) >= 3, args[3], "depth_vs_read_len.txt")
genome_size <- ifelse(length(args) >= 4, as.numeric(args[4]), 3088286401)

df <- read.delim(input_tsv, header = TRUE, stringsAsFactors = FALSE)

# guard against bad rows (e.g. duplicated header lines from concatenated
# TSVs), which would otherwise import as NA and corrupt the cumulative sum
df$length <- suppressWarnings(as.numeric(df$length))
df$count  <- suppressWarnings(as.numeric(df$count))
n_bad <- sum(is.na(df$length) | is.na(df$count))
if (n_bad > 0) {
  warning(n_bad, " row(s) had non-numeric length/count and were dropped ",
          "(check for duplicated header lines in the input TSV)")
  df <- df %>% filter(!is.na(length), !is.na(count))
}

# total bases per (sample, length) bin
# NOTE: must force double precision here. cumsum() on an integer-typed
# vector silently overflows past .Machine$integer.max (~2.1 billion) and
# returns NA for every value from that point on -- easy to hit since total
# bases per sample are typically tens of billions. sum() is overflow-safe
# internally, but cumsum() is not, hence as.numeric() below.
df <- df %>% mutate(bases = as.numeric(length) * as.numeric(count))

# for each sample, at each length x, sum bases for all reads with length >= x
# (reverse cumulative sum after sorting by length ascending)
df_cov <- df %>%
  arrange(sample_id, length) %>%
  group_by(sample_id) %>%
  mutate(cum_bases_ge = rev(cumsum(rev(as.numeric(bases))))) %>%
  ungroup()

if (!is.na(genome_size)) {
  df_cov <- df_cov %>% mutate(coverage = cum_bases_ge / genome_size)
  y_col   <- "coverage"
  y_lab   <- "Coverage (X)"
} else {
  df_cov <- df_cov %>% mutate(coverage = cum_bases_ge)
  y_col   <- "coverage"
  y_lab   <- "Cumulative bases (bp)"
}

# ---- per-sample reference metrics -----------------------------------------
# N50: the largest length x such that bases from reads >= x is at least half
# the sample's total bases (standard read-length N50 definition).
# 20X/50X length: the largest length x at which coverage is still >= that
# depth (i.e. how long reads need to be, at minimum, to retain that much
# coverage). NA if the sample never reaches that depth even at the shortest
# read length.
metrics <- df_cov %>%
  group_by(sample_id) %>%
  summarise(
    total_bases = max(cum_bases_ge),  # value at the smallest length = total
    n50_length  = max(length[cum_bases_ge >= total_bases / 2]),
    cov20_length = if (!is.na(genome_size) && any(coverage >= 20))
      max(length[coverage >= 20]) else NA_real_,
    cov50_length = if (!is.na(genome_size) && any(coverage >= 50))
      max(length[coverage >= 50]) else NA_real_,
    .groups = "drop"
  )

# cat("Per-sample reference metrics:\n")
# print(metrics)

max_length <- ceiling(max(metrics[3:5], na.rm = TRUE) / 1000) * 1000

# ---- build reference-line data for plotting --------------------------------
ref_lines <- bind_rows(
  metrics %>% transmute(sample_id, x = n50_length, metric = "N50"),
  metrics %>% filter(!is.na(cov20_length)) %>%
    transmute(sample_id, x = cov20_length, metric = "20X"),
  metrics %>% filter(!is.na(cov50_length)) %>%
    transmute(sample_id, x = cov50_length, metric = "50X")
)

metric_linetypes <- c("N50" = "dashed", "20X" = "dotted", "50X" = "dotdash")

# x-axis display limit only -- computed above from the FULL distribution so
# N50/20X/50X reference values stay correct even if the plotted range is
# truncated. Only affects what's drawn, not what's calculated.
plot_df <- df_cov
if (!is.na(max_length)) {
  plot_df <- plot_df %>% filter(length <= max_length)
  ref_lines <- ref_lines %>% filter(x <= max_length)
}

p <- ggplot(plot_df, aes(x = length, y = .data[[y_col]], color = sample_id)) +
  geom_step(linewidth = 0.6) +
  { if (nrow(ref_lines) > 0)
    geom_vline(
      data = ref_lines,
      aes(xintercept = x, color = sample_id, linetype = metric),
      linewidth = 0.4, alpha = 0.7
    )
  } +
  scale_linetype_manual(name = "Reference", values = metric_linetypes) +
  labs(
    x = "Minimum read length threshold, x (bp)",
    y = y_lab,
    color = "Sample ID",
    title = "Cumulative Coverage for Reads \u2265 x bp"
  ) +
  theme_bw()

ggsave(output_pdf, plot = p, width = 9, height = 6)
write_tsv(metrics, summary_tsv)
cat("Plot written to:", output_pdf, "\n")
if (!is.na(max_length)) {
  n_hidden <- nrow(bind_rows(
    metrics %>% transmute(sample_id, x = n50_length, metric = "N50"),
    metrics %>% filter(!is.na(cov20_length)) %>% transmute(sample_id, x = cov20_length, metric = "20X"),
    metrics %>% filter(!is.na(cov50_length)) %>% transmute(sample_id, x = cov50_length, metric = "50X")
  ) %>% filter(x > max_length))
  if (n_hidden > 0) {
    cat(n_hidden, "reference line(s) fall beyond max_length and are not shown on the plot ",
        "(see printed metrics table above for their values).\n")
  }
}
if (is.na(genome_size)) {
  cat("Note: no genome_size given, y-axis is raw bp yield, not depth (X).\n")
  cat("Pass genome size in bp as 3rd argument to convert to coverage depth ",
      "and enable the 20X/50X reference lines.\n", sep = "")
}
