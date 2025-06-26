#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript plot_maf_report.R <input_maf_file> <sample_id>")
}

input_maf <- args[1]
sample_id <- args[2]
outdir <- paste0("plots_", sample_id)
dir.create(outdir, showWarnings = FALSE)

suppressMessages(library(maftools))

# Read MAF
maf <- tryCatch({
  read.maf(maf = input_maf, vc_nonSyn = NULL)
}, error = function(e) {
  cat("⚠️ read.maf() failed:", conditionMessage(e), "\n")
  quit(status = 0)
})

# Placeholder if empty
if (nrow(maf@data) == 0) {
  cat("⚠️ No mutations found. Creating placeholder plots.\n")
  pdf(file.path(outdir, "maf_summary_plot.pdf")); plot.new(); title("No mutations found"); dev.off()
  pdf(file.path(outdir, "oncoplot_top10.pdf")); plot.new(); title("No mutations found"); dev.off()
  pdf(file.path(outdir, "titv_plot.pdf")); plot.new(); title("No mutations found"); dev.off()
  pdf(file.path(outdir, "rainfall_plot.pdf")); plot.new(); title("No mutations found"); dev.off()
  quit(status = 0)
}

# Summary plot
pdf(file.path(outdir, "maf_summary_plot.pdf"), width = 10, height = 6)
plotmafSummary(maf, dashboard = TRUE, addStat = 'median')
dev.off()

# Oncoplot
pdf(file.path(outdir, "oncoplot_top10.pdf"), width = 8, height = 6)
oncoplot(maf, top = 10)
dev.off()

# TiTv
titv_res <- titv(maf, plot = FALSE)
pdf(file.path(outdir, "titv_plot.pdf"), width = 7, height = 5)
plotTiTv(titv_res)
dev.off()

# Rainfall plot
pdf(file.path(outdir, "rainfall_plot.pdf"), width = 10, height = 4)
rainfallPlot(maf)
dev.off()

cat("✅ Basic plots completed in:", outdir, "\n")
