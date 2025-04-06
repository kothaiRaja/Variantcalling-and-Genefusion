# Load the maftools library
library(maftools)

# Read the MAF file
maf <- read.maf(maf = "sample_1.maf")

# Print summary info to console
print(maf)

# ─────────────────────────────────────────────────────
# 1. Summary Dashboard Plot
# ─────────────────────────────────────────────────────
pdf("maf_summary_plot.pdf", width = 10, height = 6)
plotmafSummary(maf = maf,
               rmOutlier = TRUE,
               addStat = 'median',
               dashboard = TRUE,
               showBarcodes = TRUE)
dev.off()

# ─────────────────────────────────────────────────────
# 2. Oncoplot: Top 10 mutated genes
# ─────────────────────────────────────────────────────
pdf("oncoplot_top10.pdf", width = 8, height = 6)
oncoplot(maf = maf,
         top = 10,
         fontSize = 12,
         showTumorSampleBarcodes = TRUE)
dev.off()

# ─────────────────────────────────────────────────────
# 3. Ti/Tv plot (Transition/Transversion)
# ─────────────────────────────────────────────────────
titv_res <- titv(maf = maf, plot = FALSE)
pdf("titv_plot.pdf", width = 7, height = 5)
plotTiTv(res = titv_res)
dev.off()

# ─────────────────────────────────────────────────────
# 4. Rainfall plot: mutation clustering
# ─────────────────────────────────────────────────────
pdf("rainfall_plot.pdf", width = 10, height = 4)
rainfallPlot(maf = maf, detectChangePoints = TRUE)
dev.off()

# ─────────────────────────────────────────────────────
# 5. Lollipop plot (for a specific gene)
# ─────────────────────────────────────────────────────
pdf("lollipop_plot_BCR.pdf", width = 10, height = 4)
lollipopPlot(maf = maf, gene = "BCR")
dev.off()
