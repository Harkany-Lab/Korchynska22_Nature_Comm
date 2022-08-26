## Nat Neuroscie TH analysis from Hemberg lab (https://hemberg-lab.github.io/scRNA.seq.datasets/mouse/brain/)
## created 170521
library(tidyverse)
library(magrittr)
library(SingleCellExperiment)
library(scater)
library(Seurat)
sce_all <- read_rds("romanov.rds")
srt_all <- as.Seurat(sce_all, counts = "counts")
Idents(srt_all) <- "cell_type2"
print(levels(srt_all))

Idents(srtTh) %>% table()
srt_all %<>%
  SCTransform(return.only.var.genes = FALSE) %>%
  RunPCA(npcs = 30) %>%
  FindNeighbors(dims = 1:30) %>%
  RunUMAP(dims = 1:30)
srtTh <- srt_all %>% subset(ident = c("Dopamine 1", "Dopamine 2 (low VMAT2)", "Dopamine 3", "Dopamine 4"))
srt_th_dge <- FindAllMarkers(srtTh, min.cells.feature = 2, logfc.threshold = 0.05)

markers_neuro_final <-
  srt_th_dge %>%
  group_by(cluster) %>%
  filter(p_val < 0.01) %>%
  top_n(n = 15, wt = pct.1) %>%
  top_n(n = 10, wt = avg_log2FC) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(pct.1)) %>%
  dplyr::distinct(gene, .keep_all = TRUE) %>%
  dplyr::arrange(cluster) %>%
  .$gene



plHMTh <- DoHeatmap(srtTh, features = c('Gad1', 'Gad2', 'Slc32a1',
                                        markers_neuro_final,
                                        'Per2', 'Nmbr', 'Nts',
                                        'Pnoc', 'Penk', 'Pdyn'), size = 3)

plDotTh <- DotPlot(srtTh, features = c('Gad1', 'Gad2', 'Slc32a1',
                                       markers_neuro_final,
                                       'Per2', 'Nmbr', 'Nts',
                                       'Pnoc', 'Penk', 'Pdyn'))+ RotatedAxis()
cowplot::save_plot(filename = 'DotPlot_Th_170521_2017-natneuro.pdf',
                   plot = plDotTh,
                   base_height = 6,
                   base_asp = 2.5)
cowplot::save_plot(filename = 'HeatMap_Th_170521_2017-natneuro.pdf',
                   plot = plHMTh,
                   base_height = 12,
                   base_asp = 1.618)
