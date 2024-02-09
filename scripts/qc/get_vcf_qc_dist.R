if (!requireNamespace("docopt", quietly = TRUE)) {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
  install.packages("docopt")
}

suppressPackageStartupMessages({
  library(VariantAnnotation)
  library(tools)
  library(data.table)
  library(ggplot2)
  library(scales)
  library(docopt)
})

convert_to_long_dt <- function(samples.attr.tab, attr.name) {
  samples.attr.dt <- as.data.table(samples.attr.tab)
  samples.attr.dt[, var_id := rownames(samples.attr.tab)]
  
  samples.attr.long.dt <- melt(samples.attr.dt, 
                               id.vars = "var_id", 
                               variable.name = "Sample", 
                               value.name = attr.name)
  
  return(samples.attr.long.dt)
}

read_vcf <- function(vcf_file_path) {
  samples.gt.long.dt <- convert_to_long_dt(readGeno(vcf_file_path, "GT"), "GT")
  samples.gq.long.dt <- convert_to_long_dt(readGeno(vcf_file_path, "GQ"), "GQ")

  samples.geno.long.dt <- samples.gt.long.dt[samples.gq.long.dt, on = .(var_id = var_id, Sample = Sample)]
  samples.geno.long.dt[is.na(GQ), GQ := -1]
  
  return(samples.geno.long.dt)
}

classify_genotype <- function(genotype) {
  allele_pair <- unlist(strsplit(genotype, "[|/]"))

  if (allele_pair[1] == allele_pair[2]) {
    if (allele_pair[1] == 0) {
      return("hom_ref")
    }
    
    return("hom_alt")
  }
  
  return("het")
}

filter_by_genotype <- function(samples.geno.long.dt, target_gt_str) {
  if (is.null(target_gt_str)) {
    return(samples.geno.long.dt)
  }
  
  target_gts <- unlist(strsplit(target_gt_str, ","))
  
  is_hom_ref <- FALSE
  is_hom_alt <- FALSE
  is_het <- FALSE
  
  for (target_gt in target_gts) {
    if (target_gt == "hom_ref") {
      is_hom_ref <- TRUE
      next
    }
    
    if (target_gt == "hom_alt") {
      is_hom_alt <- TRUE
      next
    }
    
    if (target_gt == "hom") {
      is_hom_ref <- TRUE
      is_hom_alt <- TRUE
      next
    }
    
    if (target_gt == "het") {
      is_het <- TRUE
      next
    }
  }
  
  if (!(is_hom_ref || is_hom_alt || is_het)) {
    return(samples.geno.long.dt)
  }
  
  samples.geno.long.dt[, gt.class := lapply(GT, classify_genotype)]
  
  target_gt.query <- c()
  
  if (is_hom_ref) {
    target_gt.query <- c(target_gt.query, "hom_ref")
  }
  
  if (is_hom_alt) {
    target_gt.query <- c(target_gt.query, "hom_alt")
  }
  
  if (is_het) {
    target_gt.query <- c(target_gt.query, "het")
  }
  
  samples.geno.long.dt <- samples.geno.long.dt[gt.class %in% target_gt.query]
  samples.geno.long.dt[, gt.class := NULL]
  
  return(samples.geno.long.dt)
}

get_sample_gq_dist <- function(samples.geno.long.dt) {
  max_gq <- max(samples.geno.long.dt$GQ, 99)
  
  num_of_placeholder_vars <- max_gq + 2
  
  placeholder.var_id <- paste0("var", 1:num_of_placeholder_vars)
  #placeholder.sample <- rep("placeholder", num_of_placeholder_vars)
  placeholder.gt <- rep("", num_of_placeholder_vars)
  placeholder.gq <- -1:max_gq
  
  inflated.samples.geno.long.dt <- samples.geno.long.dt
  samples.geno.long.dt <- NULL
  sample_ids <- unique(inflated.samples.geno.long.dt$Sample)
  
  for (sample_id in sample_ids) {
    geno.placeholder.dt <- data.table(var_id = placeholder.var_id,
                                      Sample = rep(sample_id, num_of_placeholder_vars),
                                      GT = placeholder.gt,
                                      GQ = placeholder.gq)
    
    inflated.samples.geno.long.dt <- rbind(inflated.samples.geno.long.dt, 
                                           geno.placeholder.dt)
  }
  
  samples.gq.dist.long.dt <- inflated.samples.geno.long.dt[, .N, .(Sample, GQ)]
  setnames(samples.gq.dist.long.dt, "N", "Count")
  samples.gq.dist.long.dt[, Count := Count - 1]
  setkey(samples.gq.dist.long.dt, Sample, GQ)
  
  return(samples.gq.dist.long.dt)
}

normalise_gq_dist <- function(samples.gq.dist.long.dt) {
  samples.gq.norm.dist.long.dt <- copy(samples.gq.dist.long.dt)
  samples.gq.norm.dist.long.dt[, Count := Count / sum(Count), Sample]
  
  return(samples.gq.norm.dist.long.dt)
}

plot_gq_dist <- function(samples.gq.dist.long.dt, is_plot_cum = FALSE, is_norm = FALSE, is_histo = FALSE) {
    gq.endpoint <- ceiling(max(samples.gq.dist.long.dt[, "GQ"]) / 10) * 100
    
    if (is_plot_cum) {
        samples.gq.dist.long.dt[, cum_Count := ave(Count, Sample, FUN = cumsum)]
        p <- ggplot(samples.gq.dist.long.dt, aes(x = GQ, y = cum_Count, color = Sample)) + 
            geom_line()
        
        if (is_norm) {
            plot.title <- "Normalised cumulative GQ distribution of sample variants"
            y_axis.title <- "Normalised cumulative frequency"
        }
        else {
            plot.title <- "Cumulative GQ distribution of sample variants"
            y_axis.title <- "Cumulative frequency"  
        }
    }
    else {
        if (is_histo) {
            p <- ggplot(samples.gq.dist.long.dt, aes(x = GQ, y = Count, fill = Sample)) + 
                geom_col(position = "dodge")
        }
        else {
            p <- ggplot(samples.gq.dist.long.dt, aes(x = GQ, y = Count, color = Sample)) + 
                geom_line()
        }
        
        if (is_norm) {
            plot.title <- "Normalised GQ distribution of sample variants"
            y_axis.title <- "Normalised frequency"
        }
        else {
            plot.title <- "GQ distribution of sample variants"
            y_axis.title <- "No. of variants"
        }
    }
    
    p <- p +
        scale_x_continuous(breaks = seq(0, gq.endpoint, 10)) +
        labs(title = plot.title, x = "GQ", y = y_axis.title) +
        theme_bw() +
        theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
            axis.title.x = element_text(size = 12, face = "bold"),
            axis.title.y = element_text(size = 12, face = "bold"),
            axis.text.x = element_text(size = 12, face = "bold"),
            axis.text.y = element_text(size = 12, face = "bold"),
            legend.position = "top",
            legend.title = element_blank(),
            legend.text = element_text(size = 9, face = "bold"))
    
    if (is_norm) {
        p <- p + scale_y_continuous()
    }
    else {
        p <- p + scale_y_continuous(labels = label_number_si())
    }
    
    return(p)
}

plot_gq_density <- function(samples.var.gq.dt) {
    work.samples.var.gq.dt <- copy(samples.var.gq.dt)
    work.samples.var.gq.dt[, var.id := 1:nrow(work.samples.var.gq.dt)]
    
    samples.var.gq.long.dt <- melt(work.samples.var.gq.dt, 
                                   id.vars = "var.id", 
                                   variable.name = "Sample", 
                                   value.name = "GQ")
    
    p <- ggplot(samples.var.gq.long.dt, aes(x = GQ, color = Sample)) +
        geom_density() + 
        labs(title = "Density Plot of GQ", x = "GQ") +
        scale_y_continuous() +
        labs(title = "GQ density of sample variants", x = "GQ", y = "No. of variants") +
        theme_bw() +
        theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
              axis.title.x = element_text(size = 12, face = "bold"),
              axis.title.y = element_text(size = 12, face = "bold"),
              axis.text.x = element_text(size = 12, face = "bold"),
              axis.text.y = element_text(size = 12, face = "bold"),
              legend.position = "top",
              legend.title = element_blank(),
              legend.text = element_text(size = 9, face = "bold"))
    
    return(p)
}

export_plot_to_png <- function(output.plot, output.png_file_path) {
    ggsave(filename = output.png_file_path, 
           plot = output.plot, 
           device = "png", 
           width = 6, 
           height = 4, 
           units = "in", 
           dpi = 300)
}

doc <- "
    Usage: get_vcf_qc_dist.sh --vcf=</path/to/VCF/file> --out=</path/to/output/directory> [--geno=<hom|hom_ref|hom_alt|het>]

    Options:
        -h --help     Show this help message and exit
        -i --vcf=</path/to/VCF/file>      VCF file path
        -o --out=</path/to/output/directory>      Output directory path
        -g --geno=<hom|hom_ref|hom_alt|het>     Target genotype(s) separated by comma(s)
"

# Parse the command-line arguments using docopt
args <- docopt(doc)

# Access the parsed arguments
vcf_file_path <- args[["--vcf"]]
output_dir_path <- args[["--out"]]
target_genotypes <- args[["--geno"]]

pdf(NULL)

message(paste0("Calculating GQ distribution for '", vcf_file_path, "'"))
samples.geno.long.dt <- read_vcf(vcf_file_path)

target_genotypes.filename_prefix <- ""

if (!is.null(target_genotypes)) {
  variant_count.unfiltered <- nrow(samples.geno.long.dt)
  samples.geno.long.dt <- filter_by_genotype(samples.geno.long.dt, target_genotypes)
  variant_count.filtered <- nrow(samples.geno.long.dt)
  
  if (variant_count.filtered < variant_count.unfiltered) {
    target_genotypes.filename_prefix <- paste0(".", gsub(",", "-", target_genotypes))
  }
}

samples.gq.dist.long.dt <- get_sample_gq_dist(samples.geno.long.dt)
samples.norm.gq.dist.long.dt <- normalise_gq_dist(samples.gq.dist.long.dt)

#GQ TSV file
output_file_name.prefix <- paste0(file_path_sans_ext(basename(vcf_file_path)), 
                                  target_genotypes.filename_prefix)
output_tsv_file_path <- file.path(output_dir_path, 
                                  paste0(output_file_name.prefix, 
                                         ".gq_dist.tsv"))
write.table(samples.gq.dist.long.dt, file = output_tsv_file_path, sep = "\t", quote = FALSE, row.names = FALSE)
message(paste0("GQ distribution exported to '", output_tsv_file_path, "'"))

#GQ distribution plot
samples.gq.dist.plot <- plot_gq_dist(samples.gq.dist.long.dt)
output_png_file_path <- file.path(output_dir_path, 
                                  paste0(output_file_name.prefix, 
                                         ".gq_dist.png"))
export_plot_to_png(samples.gq.dist.plot, output_png_file_path)
message(paste0("GQ distribution plot exported to '", output_png_file_path, "'"))

#Normalised GQ distribution plot
samples.norm.gq.dist.plot <- plot_gq_dist(samples.norm.gq.dist.long.dt, 
                                          is_norm = TRUE)
output_png_file_path <- file.path(output_dir_path, 
                                  paste0(output_file_name.prefix, 
                                         ".norm.gq_dist.png"))
export_plot_to_png(samples.norm.gq.dist.plot, output_png_file_path)
message(paste0("Normalised GQ distribution plot exported to '", 
               output_png_file_path, "'"))

#Cumulative GQ distribution plot
samples.gq.cum_dist.plot <- plot_gq_dist(samples.gq.dist.long.dt, 
                                         is_plot_cum = TRUE)
output_png_file_path <- file.path(output_dir_path, 
                                  paste0(output_file_name.prefix, 
                                         ".gq_cum_dist.png"))
export_plot_to_png(samples.gq.cum_dist.plot, output_png_file_path)
message(paste0("Cumulative GQ distribution plot exported to '", output_png_file_path, "'"))

#Normalised cumulative GQ distribution plot
samples.norm.gq.cum_dist.plot <- plot_gq_dist(samples.norm.gq.dist.long.dt, 
                                              is_plot_cum = TRUE, 
                                              is_norm = TRUE)
output_png_file_path <- file.path(output_dir_path, 
                                  paste0(output_file_name.prefix, 
                                         ".norm.gq_cum_dist.png"))
export_plot_to_png(samples.norm.gq.cum_dist.plot, output_png_file_path)
message(paste0("Normalised cumulative GQ distribution plot exported to '", 
               output_png_file_path, "'"))

#GQ density plot
# samples.gq.density.plot <- plot_gq_density(samples.var.gq.dt)
# output_png_file_path <- file.path(output_dir_path, 
#                                   paste0(output_file_name.prefix, 
#                                          ".gq_density.png"))
# export_plot_to_png(samples.gq.density.plot, output_png_file_path)
# message(paste0("GQ density plot exported to '", output_png_file_path, "'"))

#GQ distribution histogram
samples.gq.dist.histo.plot <- plot_gq_dist(samples.gq.dist.long.dt, 
                                           is_histo = TRUE)
output_png_file_path <- file.path(output_dir_path, 
                                  paste0(output_file_name.prefix, 
                                         ".gq_dist.histo.png"))
export_plot_to_png(samples.gq.dist.histo.plot, output_png_file_path)
message(paste0("Normalised GQ distribution histogram exported to '", 
               output_png_file_path, "'"))

#Normalised GQ distribution histogram
samples.norm.gq.dist.histo.plot <- plot_gq_dist(samples.norm.gq.dist.long.dt, 
                                                is_norm = TRUE,
                                                is_histo = TRUE)
output_png_file_path <- file.path(output_dir_path, 
                                  paste0(output_file_name.prefix, 
                                         ".norm.gq_dist.histo.png"))
export_plot_to_png(samples.norm.gq.dist.histo.plot, output_png_file_path)
message(paste0("Normalised GQ distribution histogram exported to '", 
               output_png_file_path, "'"))

message("Process completed")

invisible(dev.off())
