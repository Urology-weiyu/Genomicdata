# ============================================================
# Title:    Annotate SMCHD1 CUT&Tag Peaks
# Purpose:  Perform peak annotation and generate annotation pie chart
# ============================================================

setwd('/path/to/workdir')

# ---------- [Load libraries] ----------
library(ChIPseeker)
library(ggplot2)
library(dplyr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(biomaRt)

# ---------- [Read peak files] ----------
a <- list.files(pattern = "\\.narrowPeak$", full.names = TRUE)
peaks <- lapply(a, readPeakFile)
names(peaks) <- "SMCDH1"

# ---------- [Prepare promoter regions] ----------
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
promoter <- getPromoters(TxDb = txdb, upstream = 3000, downstream = 3000)

# ---------- [Tag matrix calculation] ----------
tagMatrixList <- lapply(peaks, getTagMatrix, windows = promoter)

# ---------- [Peak annotation] ----------
peakAnnoList <- lapply(peaks, annotatePeak,
                       TxDb = txdb,
                       tssRegion = c(-3000, 3000),
                       verbose = FALSE,
                       addFlankGeneInfo = TRUE,
                       flankDistance = 5000,
                       annoDb = "org.Hs.eg.db")

# ---------- [Export annotated peaks] ----------
for (i in 1:length(peakAnnoList)) {
  tmp <- as.data.frame(peakAnnoList[[i]])
  filen <- paste0(names(peaks)[i], "_annotated.bed")
  write.table(tmp, file = filen, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)
}

# ---------- [Add gene biotype information using biomaRt] ----------
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
results <- getBM(attributes = c('external_gene_name', 'gene_biotype'),
                 filters = 'external_gene_name',
                 values = tmp$SYMBOL,
                 mart = ensembl)
colnames(results)[1] <- 'SYMBOL'
tmp <- merge(tmp, results, by = 'SYMBOL')

write.table(tmp, file = filen, sep = "\t", quote = FALSE, row.names = FALSE)

# ---------- [Session Info Summary] ----------
# R version: 4.3.1 (2023-06-16)
# Platform:  x86_64-w64-mingw32/x64 (Windows 11)
# Key packages: biomaRt_2.58.2 TxDb.Hsapiens.UCSC.hg38.knownGene_3.18.0 ChIPseeker_1.38.0
