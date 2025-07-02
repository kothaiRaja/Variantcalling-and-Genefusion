#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript plot_maf_report.R <input_maf_file> <sample_id>")
}

input_maf <- args[1]
sample_id <- args[2]
outdir <- "plots"
dir.create(outdir, showWarnings = FALSE)

suppressMessages(library(maftools))

# Read MAF with error handling
maf <- tryCatch({
  read.maf(maf = input_maf, vc_nonSyn = NULL)
}, error = function(e) {
  cat("⚠️ read.maf() failed:", conditionMessage(e), "\n")
  quit(status = 0)
})

# Define output filenames
plot_files <- c(
  paste0(sample_id, "_maf_summary.pdf"),
  paste0(sample_id, "_oncoplot.pdf"),
  paste0(sample_id, "_titv.pdf"),
  paste0(sample_id, "_rainfall.pdf")
)

# Placeholder for empty MAF
if (nrow(maf@data) == 0) {
  cat("⚠️ No mutations found. Creating placeholder plots.\n")
  for (f in plot_files) {
    pdf(file.path(outdir, f))
    plot.new()
    title(paste("No mutations found -", sample_id))
    dev.off()
  }
  quit(status = 0)
}

# Generate main plots
pdf(file.path(outdir, plot_files[1]), width = 10, height = 6)
plotmafSummary(maf, dashboard = TRUE, addStat = 'median')
dev.off()

pdf(file.path(outdir, plot_files[2]), width = 8, height = 6)
oncoplot(maf, top = 10)
dev.off()

pdf(file.path(outdir, plot_files[3]), width = 7, height = 5)
plotTiTv(titv(maf, plot = FALSE))
dev.off()

pdf(file.path(outdir, plot_files[4]), width = 10, height = 4)
rainfallPlot(maf)
dev.off()

# Lollipop plot for top mutated gene (safe)
top_gene <- tryCatch({
  getGeneSummary(maf)[1, "Hugo_Symbol"]
}, error = function(e) {
  NA
})

if (!is.na(top_gene)) {
  tryCatch({
    pdf(file.path(outdir, paste0(sample_id, "_lollipop_", top_gene, ".pdf")), width = 10, height = 4)
    lollipopPlot(maf, gene = top_gene)
    dev.off()
  }, error = function(e) {
    cat(sprintf("⚠️ Skipping lollipop plot for gene %s: %s\n", top_gene, conditionMessage(e)))
  })
}


# Save gene summary as TSV
tryCatch({
  write.table(getGeneSummary(maf),
              file = file.path(outdir, paste0(sample_id, "_gene_summary.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
}, error = function(e) {
  cat("⚠️ Failed to write gene summary:", conditionMessage(e), "\n")
})

# Save session info
writeLines(capture.output(sessionInfo()),
           con = file.path(outdir, paste0(sample_id, "_sessionInfo.txt")))

# Final report
cat("✅ Plots created in:", outdir, "\n")
cat("Files created:\n")
print(list.files(outdir, full.names = TRUE))
