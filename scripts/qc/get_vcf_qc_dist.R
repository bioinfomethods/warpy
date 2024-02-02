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

read_vcf <- function(vcf_file_path) {
  samples.var.gq.dt <- as.data.table(readGeno(vcf_file_path, x = "GQ"))
  samples.ids <- colnames(samples.var.gq.dt)
  
  #Missing GQ (i.e. NA) set to -1
  samples.var.gq.dt[, (samples.ids) := lapply(.SD, function(x) replace(x, is.na(x), -1)), .SDcols = samples.ids]
  
  return(samples.var.gq.dt)
}

get_sample_gq_dist <- function(samples.var.gq.dt) {
    max_gq <- max(samples.var.gq.dt, 99)
    
    samples.ids <- colnames(samples.var.gq.dt)
    sort.samples.ids <- paste0("S", 1:length(samples.ids))
    colnames(samples.var.gq.dt) <- sort.samples.ids
    
    gq.placeholder.dt <- as.data.table(matrix(-1:max_gq, ncol = length(samples.ids), nrow = max_gq + 2))
    colnames(gq.placeholder.dt) <- sort.samples.ids
    inflated.samples.var.gq.dt <- rbind(samples.var.gq.dt, gq.placeholder.dt)
    samples.var.gq.dt <- NULL
    
    vcf.samples.gq.dist.dt <- data.table("GQ" = -1:max_gq)
    
    for (sample.id in sort.samples.ids) {
        sample.gq.dist.dt <- inflated.samples.var.gq.dt[, .N, by = sample.id]
        
        eval(parse(text = paste0("setkey(sample.gq.dist.dt, ", sample.id, ")")))
        sample.gq.dist.dt[, N := N - 1]
        
        eval(parse(text = paste0("vcf.samples.gq.dist.dt[, ", sample.id, " := sample.gq.dist.dt[, N]]")))
    }
    
    colnames(vcf.samples.gq.dist.dt) <- c("GQ", samples.ids)
    
    return(vcf.samples.gq.dist.dt)
}

normalise_gq_dist <- function(vcf.samples.gq.dist.dt) {
    vcf.samples.norm.gq.dist.dt <- copy(vcf.samples.gq.dist.dt)
    samples.ids <- colnames(vcf.samples.norm.gq.dist.dt)[-1]
    #print(samples.ids)
    vcf.samples.norm.gq.dist.dt[, (samples.ids) := lapply(.SD, function(x) x / sum(x)), .SDcols = samples.ids]
    
    return(vcf.samples.norm.gq.dist.dt)
}

plot_gq_dist <- function(vcf.samples.gq.dist.dt, is_plot_cum = FALSE, is_norm = FALSE) {
    samples.gq.dist.long.dt <- melt(vcf.samples.gq.dist.dt, 
                                    id.vars = "GQ", 
                                    variable.name = "Sample", 
                                    value.name = "Count")
    
    gq.endpoint <- ceiling(max(samples.gq.dist.long.dt[, "GQ"]) / 10) * 100
    
    if (is_plot_cum) {
        samples.gq.dist.long.dt[, cum_N := ave(Count, Sample, FUN = cumsum)]
        p <- ggplot(samples.gq.dist.long.dt, aes(x = GQ, y = cum_N, color = Sample))
        
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
        p <- ggplot(samples.gq.dist.long.dt, aes(x = GQ, y = Count, color = Sample))
        
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
        geom_line() +
        scale_x_continuous(breaks = seq(0, gq.endpoint, 10)) +
        #scale_y_continuous(labels = label_number_si()) +
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
    Usage: get_vcf_qc_dist.sh --vcf=</path/to/VCF/file> --out=</path/to/output/directory>

    Options:
        -h --help     Show this help message and exit.
        -i --vcf=</path/to/VCF/file> VCF file path.
        -o --out=</path/to/output/directory> Output directory path.
"

# Parse the command-line arguments using docopt
args <- docopt(doc)

# Access the parsed arguments
vcf_file_path <- args[["--vcf"]]
output_dir_path <- args[["--out"]]

pdf(NULL)

message(paste0("Calculating GQ distribution for '", vcf_file_path, "'"))

samples.var.gq.dt <- read_vcf(vcf_file_path)
vcf.samples.gq.dist.dt <- get_sample_gq_dist(samples.var.gq.dt)
vcf.samples.norm.gq.dist.dt <- normalise_gq_dist(vcf.samples.gq.dist.dt)

samples.gq.dist.plot <- plot_gq_dist(vcf.samples.gq.dist.dt)
samples.gq.cum_dist.plot <- plot_gq_dist(vcf.samples.gq.dist.dt, is_plot_cum = TRUE)
samples.norm.gq.dist.plot <- plot_gq_dist(vcf.samples.norm.gq.dist.dt, is_norm = TRUE)
samples.norm.gq.cum_dist.plot <- plot_gq_dist(vcf.samples.norm.gq.dist.dt, is_plot_cum = TRUE, is_norm = TRUE)
samples.gq.density.plot <- plot_gq_density(samples.var.gq.dt)

#GQ TSV file
output_file_name.prefix <- file_path_sans_ext(basename(vcf_file_path))
output_tsv_file_path <- file.path(output_dir_path, 
                                  paste0(output_file_name.prefix, 
                                         ".gq_dist.tsv"))
write.table(vcf.samples.gq.dist.dt, file = output_tsv_file_path, sep = "\t", quote = FALSE, row.names = FALSE)
message(paste0("GQ distribution exported to '", output_tsv_file_path, "'"))

#GQ distribution plot
gq.dist.output_png_file_path <- file.path(output_dir_path, 
                                          paste0(output_file_name.prefix, 
                                                 ".gq_dist.png"))
export_plot_to_png(samples.gq.dist.plot, gq.dist.output_png_file_path)
message(paste0("GQ distribution plot exported to '", gq.dist.output_png_file_path, "'"))

#Cumulative GQ distribution plot
gq.cum_dist.output_png_file_path <- file.path(output_dir_path, 
                                              paste0(output_file_name.prefix, 
                                                     ".gq_cum_dist.png"))
export_plot_to_png(samples.gq.cum_dist.plot, gq.cum_dist.output_png_file_path)
message(paste0("Cumulative GQ distribution plot exported to '", gq.cum_dist.output_png_file_path, "'"))

#Normalised GQ distribution plot
norm.gq.dist.output_png_file_path <- file.path(output_dir_path, 
                                               paste0(output_file_name.prefix, 
                                                      ".norm.gq_dist.png"))
export_plot_to_png(samples.norm.gq.dist.plot, norm.gq.dist.output_png_file_path)
message(paste0("Normalised GQ distribution plot exported to '", 
               norm.gq.dist.output_png_file_path, "'"))

#Normalised cumulative GQ distribution plot
norm.gq.cum_dist.output_png_file_path <- file.path(output_dir_path, 
                                                   paste0(output_file_name.prefix, 
                                                          ".norm.gq_cum_dist.png"))
export_plot_to_png(samples.norm.gq.cum_dist.plot, norm.gq.cum_dist.output_png_file_path)
message(paste0("Normalised cumulative GQ distribution plot exported to '", 
               norm.gq.cum_dist.output_png_file_path, "'"))

#GQ density plot
gq.density.output_png_file_path <- file.path(output_dir_path, 
                                             paste0(output_file_name.prefix, 
                                                    ".gq_density.png"))
export_plot_to_png(samples.gq.density.plot, gq.density.output_png_file_path)
message(paste0("GQ density plot exported to '", gq.dist.output_png_file_path, "'"))

message("Process completed")

invisible(dev.off())
