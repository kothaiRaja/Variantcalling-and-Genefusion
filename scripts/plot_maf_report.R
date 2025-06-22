#!/usr/bin/env Rscript

# ─────────────────────────────────────────────────────
# 0. Setup and Argument Parsing
# ─────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript plot_maf_report.R <input_maf_file> <sample_id>")
}

input_maf <- args[1]
sample_id <- args[2]

cat("Starting MAF visualization for:", sample_id, "\n")
cat("Reading file:", input_maf, "\n")

# ─────────────────────────────────────────────────────
# 1. Load Library
# ─────────────────────────────────────────────────────
suppressMessages(library(maftools))


# ─────────────────────────────────────────────────────
# 2. Read MAF File (Allow all variant classes)
# ─────────────────────────────────────────────────────
# ─────────────────────────────────────────────────────
# 2. Read MAF File (Allow all variant classes)
# ─────────────────────────────────────────────────────
maf <- tryCatch({
  read.maf(maf = input_maf, vc_nonSyn = NULL)
}, error = function(e) {
  cat("⚠️ read.maf() failed:", conditionMessage(e), "\n")
  quit(status = 0)
})

# Check for mutations and create placeholders if needed
if (nrow(maf@data) == 0 || nrow(getGeneSummary(maf)) == 0) {
  cat("⚠️ No mutations found in MAF file. Creating empty placeholder plots.\n")
  outdir <- paste0("plots_", sample_id)
  dir.create(outdir, showWarnings = FALSE)

  # Create blank PDFs as placeholders
  pdf(file.path(outdir, "maf_summary_plot.pdf")); plot.new(); title("No mutations found"); dev.off()
  pdf(file.path(outdir, "oncoplot_top10.pdf")); plot.new(); title("No mutations found"); dev.off()
  pdf(file.path(outdir, "titv_plot.pdf")); plot.new(); title("No mutations found"); dev.off()
  pdf(file.path(outdir, "rainfall_plot.pdf")); plot.new(); title("No mutations found"); dev.off()
  pdf(file.path(outdir, "lollipop_plot_placeholder.pdf")); plot.new(); title("No mutations found"); dev.off()

  cat(" Placeholder plots generated in:", outdir, "\n")
  quit(status = 0)
}


# Create output directory
outdir <- paste0("plots_", sample_id)
dir.create(outdir, showWarnings = FALSE)

# ─────────────────────────────────────────────────────
# 3. Summary Dashboard Plot
# ─────────────────────────────────────────────────────
pdf(file.path(outdir, "maf_summary_plot.pdf"), width = 10, height = 6)
plotmafSummary(maf = maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, showBarcodes = TRUE)
dev.off()

# ─────────────────────────────────────────────────────
# 4. Oncoplot for Top 10 Genes
# ─────────────────────────────────────────────────────
pdf(file.path(outdir, "oncoplot_top10.pdf"), width = 8, height = 6)
oncoplot(maf = maf, top = 10, fontSize = 12, showTumorSampleBarcodes = TRUE)
dev.off()

# ─────────────────────────────────────────────────────
# 5. Transition/Transversion Plot
# ─────────────────────────────────────────────────────
titv_res <- titv(maf = maf, plot = FALSE)
pdf(file.path(outdir, "titv_plot.pdf"), width = 7, height = 5)
plotTiTv(res = titv_res)
dev.off()

# ─────────────────────────────────────────────────────
# 6. Rainfall Plot
# ─────────────────────────────────────────────────────
pdf(file.path(outdir, "rainfall_plot.pdf"), width = 10, height = 4)
rainfallPlot(maf = maf, detectChangePoints = TRUE)
dev.off()

# ─────────────────────────────────────────────────────
# 7. Lollipop Plot (Gene-Level)
# ─────────────────────────────────────────────────────
genes <- getGeneSummary(maf)$Hugo_Symbol[1:1]  # Top mutated gene

for (g in genes) {
  cat("Processing gene:", g, "\n")
  pdf(file.path(outdir, paste0("lollipop_plot_", g, ".pdf")), width = 10, height = 4)

  tryCatch({
    lollipopPlot(maf = maf, gene = g)
  }, error = function(e) {
    cat("⚠️ Skipping gene", g, "- no structure found or other issue:", conditionMessage(e), "\n")
    plot.new()
    title(paste("Lollipop plot not available for", g))
  })

  dev.off()
}
