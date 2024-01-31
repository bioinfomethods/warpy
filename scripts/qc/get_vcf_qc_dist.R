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

get_sample_gq_dist <- function(vcf_file_path) {
    samples.var.gq.dt <- as.data.table(readGeno(vcf_file_path, x = "GQ"))
    samples.ids <- colnames(samples.var.gq.dt)
    
    #Missing GQ (i.e. NA) set to -1
    samples.var.gq.dt[, (samples.ids) := lapply(.SD, function(x) replace(x, is.na(x), -1)), .SDcols = samples.ids]
    
    max_gq <- max(samples.var.gq.dt, 99)
    
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

plot_gq_dist <- function(vcf.samples.gq.dist.dt) {
    samples.gq.dist.long.dt <- melt(vcf.samples.gq.dist.dt, 
                                    id.vars = "GQ", 
                                    variable.name = "Sample", 
                                    value.name = "Count")
    
    gq.endpoint <- ceiling(max(samples.gq.dist.long.dt[, "GQ"]) / 10) * 100
    
    p <- ggplot(samples.gq.dist.long.dt, aes(x = GQ, y = Count, color = Sample)) +
        geom_line() +
        scale_x_continuous(breaks = seq(0, gq.endpoint, 10)) +
        scale_y_continuous(labels = label_number_si()) +
        labs(title = "GQ distribution of sample variants", x = "GQ", y = "No. of variants") +
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

vcf.samples.gq.dist.dt <- get_sample_gq_dist(vcf_file_path)
samples.gq.dist.plot <- plot_gq_dist(vcf.samples.gq.dist.dt)

output_file_name.prefix <- file_path_sans_ext(basename(vcf_file_path))

output_tsv_file_path <- file.path(output_dir_path, 
                                  paste0(output_file_name.prefix, 
                                         ".gq_dist.tsv"))
write.table(vcf.samples.gq.dist.dt, file = output_tsv_file_path, sep = "\t", quote = FALSE, row.names = FALSE)

message(paste0("GQ distribution exported to '", output_tsv_file_path, "'"))

output_png_file_path <- file.path(output_dir_path, 
                                  paste0(output_file_name.prefix, 
                                         ".gq_dist.png"))
ggsave(filename = output_png_file_path, 
       plot = samples.gq.dist.plot, 
       device = "png", 
       width = 6, 
       height = 4, 
       units = "in", 
       dpi = 300)

message(paste0("GQ distribution plot exported to '", output_png_file_path, "'"))
message("Process completed")

invisible(dev.off())
