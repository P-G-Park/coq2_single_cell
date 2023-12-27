rm(list=ls())
pacman::p_load(tidyverse, readxl, Seurat, data.table, ggsci, ggpubr, Matrix, harmony)

setwd('..')

####  Define QC parameters  #######

Matrix_disease <- Read10X('./raw_data/81230608')
seurat_disease <- CreateSeuratObject(Matrix_han, min.cells = 3)
seurat_disease$disease <- 'disease'
seurat_disease$channel <- 'this_paper'

all <- merge(seurat_disease, list(control_1, control_2, control_3, control_4, control_5, control_6, control_7, control_8))
all$percent.mt <- PercentageFeatureSet(object = all, pattern = "^MT-")
all <- subset(all, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 50)

all <- NormalizeData(all)
all <- FindVariableFeatures(all)
all <- ScaleData(all)
all <- RunPCA(all)
all <- IntegrateLayers(
  object = all, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "integrated"
)
all[['RNA']] <- JoinLayers(all[["RNA"]])
all <- all %>% 
  RunUMAP(reduction = 'integrated', dims = 1:30) %>% 
  FindNeighbors(reduction = 'integrated',  dims = 1:30) %>% 
  FindClusters(resolution = 0.4)

# save(all, file = './raw_data/all_v1.RData')


##### 230905
load('./raw_data/all_v1.RData')

quick(all)
quickdot(all, feat = Marker)
all <-  all %>% RenameIdents(
  '0' = 'PT',
  '1' = 'Injured PT',
  '15' = 'Injured PT',
  '5' = 'TAL',
  '7' ='DCT',
  '2' = 'CD-P',
  '10' = 'CD-P',
  '8' = 'CD-I',
  '12' = 'IC',
  '9' = 'gEC',
  '3' = 'ptEC',
  '13' = 'Arterial EC',
  '11' = 'SMC',
  '14' = 'SMC',
  '18' = 'PERI',
  '19' = 'PODO',
  '6' = 'Mac',
  '4' = 'T',
  '16' = 'T',
  '17' = 'B'
)
Marker <- c('CUBN','VCAM1', 'UMOD', 
                'SLC12A1', 'SLC12A3', 'CALB1', 'AQP2',
                'SPINK1',
                'SLC26A4',
                'CD34', 'EHD3', 'PLVAP', 'CLDN5',
                'ACTA2', 'PDGFRB',
                'NPHS1', 'NPHS2',
                  'C1QC',  'S100A8','CD3D', 'CD79A',
                'PTPRC')

DimPlot(all, group.by = 'disease', reduction = "umap")+
  theme(plot.title = element_blank())

DimPlot(all, reduction = "umap", label = TRUE)+
    theme(plot.title = element_blank()) & NoLegend()

DotPlot(all, features = rev(Marker)) +
  theme(axis.text.x = element_text(hjust = 1,angle = 30))+
  coord_flip() & NoLegend()


degs <- tibble('Type' = character(),'Control' = integer(), 'Disease' = integer())
for(i in levels(Idents(all))){
  obj <- subset(all, idents = i)
  Idents(obj) <- obj$disease
  obj_deg <- FindMarkers(obj, ident.1 = 'disease')
  obj_deg_pos <- obj_deg %>% filter(avg_log2FC > 0.5, p_val_adj < 0.05)
  obj_deg_neg <- obj_deg %>% filter(avg_log2FC < -0.5, p_val_adj < 0.05)
  degs <- degs %>% 
    add_row(Type = i, Control = nrow(obj_deg_neg), Disease = nrow(obj_deg_pos))
}
degs %>% write.table("clipboard", sep = '\t') 



require(writexl)
DEGs <- list()
for(i in levels(Idents(all))){
  subset_i <- subset(all, idents = i)
  Idents(subset_i) <- subset_i$disease
  DEGs_i <- FindMarkers(subset_i, ident.1 = 'disease',  only.pos = TRUE) %>% 
    rownames_to_column('Gene')
  DEGs[[i]] <- DEGs_i
}
writexl::write_xlsx(DEGs, path = str_c('all_disease.xlsx'))

cell.prop <- table(Idents(all), all$disease) %>% 
  as.data.frame() %>% pivot_wider(names_from = Var2, values_from = Freq)

tibble(type = c('Control', 'Disease'), Parenchyme = c(6721, 2896), Immune = c(1064, 841)) %>% 
  pivot_longer(cols = c('Parenchyme', 'Immune')) %>% 
  ggplot(aes(x = type, fill = name, y = value)) +
  geom_bar(stat = 'identity', position = 'fill', color = 'black', width = .5, linewidth = 1) +
  scale_fill_manual(values = c('brown1', 'white')) +
  ggpubr::theme_pubr() +
  theme(legend.title = element_blank()) + xlab('') + ylab('Proportion')


podo <- subset(all, idents = 'PODO')
Idents(podo) <- podo$disease
podo_deg <- FindMarkers(podo, ident.1 = 'disease')

pt <- subset(all, idents = 'PT')
gec <- subset(all, idents = 'gEC')
ptec <- subset(all, idents = 'ptEC')
krm <- subset(all, idents = 'Mac')

coq <-  c('PDSS1', 'PDSS2', 'COQ2', 'COQ3', 'COQ4', 'COQ5', 'COQ6', 'COQ7',
          'COQ8A', 'COQ8B', 'COQ9', 'COQ10A', 'COQ10B', 'FDX1L', 'FDXR', 'ALDH3A1')
mt <- c('MT-ND2', 'MT-ND4', 'MT-ND5', 'MT-ND1', 'MT-ND6', 'MT-ND3', 'MT-CO1', 'MT-ATP6', 'MT-ATP8')

c1 <- c('NDUFA1', 'NDUFA2', 'NDUFA3', 'NDUFB1', 'NDUFB2')
c2 <- c('SDHA', 'SDHB', 'SDHC')
c3 <- c('UQCRB', 'UQCRQ', 'CYC1')
c4 <- c('COX5A', 'COX5B', 'COX6C', 'COX8A', 'COX7B')
c5 <- c('ATP5PD','ATP5PF','ATP5MG', 'ATP5MC2','ATP5MC3', 'ATP5ME', 'ATP5F1E')
etc <- c(c1, c2, c3, c4, c5)

glycolysis <- c('ALDOA', 'BPGM', 'ENO1', 'ENO2', 'GAPDH', 'GPI', 'HK1', 'HK2', 
                'PFKL', 'PFKM', 'PGAM1', 'PGK1', 'PKM', 'TPI1')


mt_fc <- FoldChange(all, ident.1 = 'disease', feat = mt, group.by = 'disease', subset.ident = 'PT')
mt_fc <- mt_fc[0]
for(i in levels(Idents(all))){
  asdf <- FoldChange(all, ident.1 = 'disease', feat = mt, group.by = 'disease', subset.ident = i)
  asdf <- asdf[1]
  colnames(asdf) <- i
  mt_fc <- cbind(mt_fc, asdf)
}


coq_fc <- FoldChange(all, ident.1 = 'disease', feat = coq, group.by = 'disease', subset.ident = 'PT')
coq_fc <- coq_fc[0]
for(i in levels(Idents(all))){
  asdf <- FoldChange(all, ident.1 = 'disease', feat = coq, group.by = 'disease', subset.ident = i)
  asdf <- asdf[1]
  colnames(asdf) <- i
  coq_fc <- cbind(coq_fc, asdf)
}


etc_fc <- FoldChange(all, ident.1 = 'disease', feat = etc, group.by = 'disease', subset.ident = 'PT')
etc_fc <- etc_fc[0]
for(i in levels(Idents(all))){
  asdf <- FoldChange(all, ident.1 = 'disease', feat = etc, group.by = 'disease', subset.ident = i)
  asdf <- asdf[1]
  colnames(asdf) <- i
  etc_fc <- cbind(etc_fc, asdf)
}

glyco_fc <- FoldChange(all, ident.1 = 'disease', feat = glycolysis, group.by = 'disease', subset.ident = 'PT')
glyco_fc <- glyco_fc[0]
for(i in levels(Idents(all))){
  asdf <- FoldChange(all, ident.1 = 'disease', feat = glycolysis, group.by = 'disease', subset.ident = i)
  asdf <- asdf[1]
  colnames(asdf) <- i
  glyco_fc <- cbind(glyco_fc, asdf)
}

ComplexHeatmap::Heatmap((coq_fc %>% dplyr::select(-14) %>% as.matrix()), 
                        cluster_rows = FALSE, cluster_columns = FALSE, name = 'Fold Change')

f1 = colorRamp2::colorRamp2(c(-3.5, 0, 3.5), c("blue", "#EEEEEE", "red"))
ComplexHeatmap::Heatmap((mt_fc %>% dplyr::select(-14) %>% as.matrix()), 
                        col = f1, cluster_rows = FALSE, cluster_columns = FALSE, name = 'Fold Change')
ComplexHeatmap::Heatmap((etc_fc %>% dplyr::select(-14) %>% as.matrix()), cluster_rows = FALSE, cluster_columns = FALSE, name = 'Fold Change')

f1 = colorRamp2::colorRamp2(c(-3.5, 0, 3.5), c("blue", "#EEEEEE", "red"))
ComplexHeatmap::Heatmap((glyco_fc %>% dplyr::select(-14) %>% as.matrix()), col = f1, 
                        cluster_rows = FALSE, cluster_columns = FALSE, name = 'Fold Change')


kid_dev <- c('PODXL', 'NPHS1', 'NPHS2', 'STAT1', 'PLCE1', 'SYNPO', 'COL6A2', 'FBLN1','FBN1', 
             'COL5A1')

kid_fc <- FoldChange(all, ident.1 = 'disease', feat = kid_dev, group.by = 'disease', subset.ident = 'PT')
kid_fc <- kid_fc[0]
for(i in levels(Idents(all))){
  asdf <- FoldChange(all, ident.1 = 'disease', feat = kid_dev, group.by = 'disease', subset.ident = i)
  asdf <- asdf[1]
  colnames(asdf) <- i
  kid_fc <- cbind(kid_fc, asdf)
}
ComplexHeatmap::Heatmap((kid_fc %>% select(-14) %>% as.matrix()), 
                        cluster_rows = FALSE, cluster_columns = FALSE,  col = f1,
                        name = 'Fold Change')

pacman::p_load(fgsea, msigdbr, ggsci)

podo_gene <- podo_deg$avg_log2FC
names(podo_gene) <- rownames(podo_deg)

h_pathway <- msigdbr(species = "human", category = "C5", subcategory = 'GO:BP')
h_list <- split(x = h_pathway$gene_symbol, f = h_pathway$gs_name)

gsea_H <- fgsea(pathways=h_list, stats=podo_gene, nperm=10000) %>% 
  as_tibble() %>% 
  mutate(pathway = str_replace(pathway,'GOBP_', ''),
         pathway = str_replace(pathway, '_', ' '),
         pathway = str_replace(pathway, '_', ' '),
         pathway = str_replace(pathway, '_', ' '),
         pathway = str_replace(pathway, '_', ' ')) %>% 
  arrange(-NES) %>% filter(!is.na(NES))

ggplot(gsea_H %>% slice(1:7, (nrow(gsea_H)-6):nrow(gsea_H)), aes(reorder(pathway, NES), NES)) +
  geom_col(fill = pal_nejm(alpha = 0.8)(2)[2]) +
  coord_flip()+ xlab('')+
  theme_minimal()
# donut chart

all$ident1 <- recode(Idents(all),
                             B = 'Imm',
                             Neutrophil = 'Imm',
                             'CD4+T' = 'Imm',
                             'KRM'= 'Imm',
                             'CD8+T' = 'Imm',
                             'gEC' = 'EC',
                             'ptEC' = 'EC',
                             'Arterial EC' = 'EC'
)

