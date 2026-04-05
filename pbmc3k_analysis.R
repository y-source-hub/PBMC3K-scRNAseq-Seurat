# ============================================
# PBMC3K 数据分析（文件已手动下载）
# ============================================

# 1. 加载包
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

# 2. 确保文件夹存在
dir.create("data", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

# 3. 检查文件是否已存在
if (file.exists("data/pbmc3k_filtered_gene_bc_matrices.tar.gz")) {
  message("✓ 数据文件已找到，开始解压...")
  
  # 4. 解压
  untar("data/pbmc3k_filtered_gene_bc_matrices.tar.gz", exdir = "data/")
  
  # 5. 读取数据
  pbmc.data <- Read10X(data.dir = "data/filtered_gene_bc_matrices/hg19/")
  message(paste("✓ 数据读取成功！维度：", paste(dim(pbmc.data), collapse = " x ")))
  
  # 6. 创建 Seurat 对象
  pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k",
                             min.cells = 3, min.features = 200)
  
  # 7. 质控
  pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
  pbmc <- subset(pbmc, subset = nFeature_RNA > 200 &
                   nFeature_RNA < 2500 &
                   percent.mt < 5)
  
  # 8. 标准化
  pbmc <- NormalizeData(pbmc)
  
  # 9. 找高变基因
  pbmc <- FindVariableFeatures(pbmc, nfeatures = 2000)
  
  # 10. 缩放
  pbmc <- ScaleData(pbmc)
  
  # 11. PCA
  pbmc <- RunPCA(pbmc, npcs = 50)
  
  # 12. 聚类
  pbmc <- FindNeighbors(pbmc, dims = 1:15)
  pbmc <- FindClusters(pbmc, resolution = 0.5)
  
  # 13. UMAP
  pbmc <- RunUMAP(pbmc, dims = 1:15)
  
  # 14. 细胞类型注释
  new.cluster.ids <- c(
    "Naive CD4+ T", "CD14+ Mono", "Naive CD4+ T",
    "Memory CD4+ T", "B cell", "CD8+ T",
    "FCGR3A+ Mono", "NK cell", "DC", "Platelet"
  )
  names(new.cluster.ids) <- levels(pbmc)
  pbmc <- RenameIdents(pbmc, new.cluster.ids)
  
  # 15. 画图
  umap_plot <- DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
  print(umap_plot)
  
  # 16. 保存结果
  saveRDS(pbmc, file = "results/pbmc3k_final.rds")
  ggsave("results/umap_pbmc3k_final.png", plot = umap_plot, width = 10, height = 8, dpi = 300)
  
  message("\n🎉 分析完成！结果保存在 results/ 文件夹")
  
} else {
  message("✗ 数据文件不存在！")
  message("请手动下载文件并放到 data/ 文件夹中：")
  message("https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz")
}

# 重新绘制并保存 UMAP 图
umap_plot <- DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave("results/umap_pbmc3k_final.png", plot = umap_plot, width = 10, height = 8, dpi = 300)

# 保存标记基因表
if (!file.exists("results/marker_genes.csv")) {
  pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  write.csv(pbmc.markers, file = "results/marker_genes.csv")
}

message("✓ 结果文件已保存到 results/ 文件夹")

#补充火山图
library(Seurat)
unique(Idents(pbmc))

# 1. 加载包
library(Seurat)
library(dplyr)
library(ggplot2)

# 2. 加载数据（如果还没加载）
pbmc <- readRDS("results/pbmc3k_final.rds")

# 3. 确认细胞类型名称
unique(Idents(pbmc))

# 4. 比较 CD14+ Mono 和 FCGR3A+ Mono
monocyte_markers <- FindMarkers(pbmc, 
                                ident.1 = "CD14+ Mono", 
                                ident.2 = "FCGR3A+ Mono",
                                min.pct = 0.25,
                                logfc.threshold = 0.25)

# 5. 查看结果
head(monocyte_markers)

# 6. 添加基因名列
monocyte_markers$gene <- rownames(monocyte_markers)

# 7. 标记显著性
monocyte_markers$sig <- ifelse(monocyte_markers$p_val_adj < 0.05 & 
                                 abs(monocyte_markers$avg_log2FC) > 0.5, 
                               "Significant", "Not Significant")

# 8. 画火山图
volcano_plot <- ggplot(monocyte_markers, 
                       aes(x = avg_log2FC, 
                           y = -log10(p_val_adj), 
                           color = sig)) +
  geom_point(alpha = 0.6, size = 1) +
  scale_color_manual(values = c("grey", "red")) +
  theme_minimal() +
  labs(title = "CD14+ Mono vs FCGR3A+ Mono",
       x = "Log2 Fold Change",
       y = "-log10(Adjusted P-value)") +
  theme(legend.position = "bottom")

# 9. 显示并保存
print(volcano_plot)
ggsave("results/volcano_plot.png", plot = volcano_plot, width = 8, height = 6, dpi = 300)

message("✅ 火山图已保存到 results/volcano_plot.png")


# ============================================
# PBMC3K 单细胞数据分析 - 补充分析
# 包括：热图、GO富集、小提琴图、特征图、山脊线图
# ============================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

# 加载数据
pbmc <- readRDS("results/pbmc3k_final.rds")

# 创建 results 文件夹
dir.create("results", showWarnings = FALSE)

# ============================================
# 1. 热图（展示 top 差异基因）
# ============================================

# 获取所有群的 top5 标记基因
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, 
                               min.pct = 0.25, 
                               logfc.threshold = 0.25)

top5_markers <- pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

# 绘制热图
heatmap_plot <- DoHeatmap(pbmc, features = top5_markers$gene) + 
  NoLegend() +
  labs(title = "Top 5 Marker Genes per Cell Type")

ggsave("results/heatmap.png", plot = heatmap_plot, width = 12, height = 10, dpi = 300)

# ============================================
# 2. GO 富集分析（需要安装包）
# ============================================

install.packages("BiocManager")

# 加载 BiocManager
library(BiocManager)

BiocManager::install("clusterProfiler", ask = FALSE)
BiocManager::install("org.Hs.eg.db", ask = FALSE)

install.packages("cachem")
install.packages("RSQLite")

install("org.Hs.eg.db")

BiocManager::install("GO.db")

# 加载包
library(Seurat)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)

# 加载数据
pbmc <- readRDS("results/pbmc3k_final.rds")

# 找差异基因（B细胞 vs CD8+ T细胞）
b_vs_t_markers <- FindMarkers(pbmc, 
                              ident.1 = "B cell", 
                              ident.2 = "CD8+ T",
                              min.pct = 0.25,
                              logfc.threshold = 0.5)

# 查看差异基因数量
sig_genes <- b_vs_t_markers %>% filter(p_val_adj < 0.05)
cat("显著差异基因数量:", nrow(sig_genes), "\n")

# 提取上调基因（在B细胞中高表达的基因）
up_genes <- b_vs_t_markers %>%
  filter(avg_log2FC > 1 & p_val_adj < 0.05) %>%
  rownames()

cat("用于分析的基因数量:", length(up_genes), "\n")

# 过滤掉核糖体蛋白基因
# 核糖体蛋白基因通常以 RP、RPL、RPS 开头
up_genes_filtered <- up_genes[!grepl("^RP[LS]|^RPL|^RPS", up_genes)]

cat("过滤后的基因数量:", length(up_genes_filtered), "\n")

# 转换基因名为 Entrez ID
entrez_ids <- bitr(up_genes, 
                   fromType = "SYMBOL", 
                   toType = "ENTREZID", 
                   OrgDb = org.Hs.eg.db)

# GO 富集分析
go_result <- enrichGO(gene = entrez_ids$ENTREZID,
                      OrgDb = org.Hs.eg.db,
                      ont = "BP",           # 生物学过程
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.10)

# 查看结果
cat("显著富集的通路数量:", nrow(go_result), "\n")
head(go_result, 10)

# 绘制并保存 dotplot
if (nrow(go_result) > 0) {
  # 使用更美观的绘图方式
  go_plot <- dotplot(go_result, showCategory = 15) + 
    ggtitle("GO Enrichment Analysis: B cells vs CD8+ T cells (Ribosomal genes filtered)") +
    theme(plot.title = element_text(hjust = 0.5, size = 10))
  
  print(go_plot)
  ggsave("results/go_enrichment.png", plot = go_plot, width = 12, height = 10, dpi = 300)
  
  # 保存完整结果表格
  write.csv(as.data.frame(go_result), "results/go_enrichment_results.csv")
  
  message("\n✅ 分析完成！结果已保存到 results/go_enrichment_improved.png")
} else {
  message("没有显著的富集结果，可以尝试放宽 pvalueCutoff 到 0.1")
}

# ============================================
# 3. 小提琴图（展示关键基因在各细胞类型的表达分布）
# ============================================

# 加载包
library(Seurat)
library(ggplot2)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)

# 加载数据
pbmc <- readRDS("results/pbmc3k_final.rds")

# 创建 results 文件夹
dir.create("results", showWarnings = FALSE)

# 正确的标签对应关系（之前有错误标注）

# 重新标注
correct_ids <- c(
  "CD14+ Mono",        # 原 cluster 0
  "Naive CD4+ T",      # 原 cluster 1
  "B cell",            # 原 cluster 2
  "CD8+ T",            # 原 cluster 3
  "FCGR3A+ Mono",      # 原 cluster 4
  "NK cell",           # 原 cluster 5
  "NK cell",           # 原 cluster 6（如果标记基因是 NK 相关）
  "DC"                 # 原 cluster 7
)

# 应用新标签
names(correct_ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, correct_ids)

# 验证
unique(Idents(pbmc))

# 重新找标记基因确认
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, 
                               min.pct = 0.25, 
                               logfc.threshold = 0.25)

# 查看新的 top 基因
top5 <- pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

print(top5[, c("cluster", "gene")], n = 50)

# 保存修正后的对象
saveRDS(pbmc, file = "results/pbmc3k_final.rds")

# ============================================
# 重新生成所有图片和表格（基于修正后的标签）
# ============================================

library(Seurat)
library(ggplot2)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)

# 加载修正后的数据
pbmc <- readRDS("results/pbmc3k_final.rds")

# 创建 results 文件夹
dir.create("results", showWarnings = FALSE)

# ============================================
# 重新生成标记基因表 (marker_genes.csv)
# ============================================

message("正在重新生成标记基因表...")

pbmc.markers <- FindAllMarkers(pbmc,
                               only.pos = TRUE,
                               min.pct = 0.25,
                               logfc.threshold = 0.25)

# 保存完整标记基因表
write.csv(pbmc.markers, file = "results/marker_genes.csv", row.names = FALSE)

# 同时保存每个细胞类型 top 10 标记基因（方便查看）
top10_markers <- pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)
write.csv(top10_markers, file = "results/marker_genes_top10.csv", row.names = FALSE)

message("✅ marker_genes.csv 已生成")

# ============================================
# UMAP 图
# ============================================

umap_plot <- DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave("results/umap_pbmc3k_final.png", plot = umap_plot, width = 10, height = 8, dpi = 300)

# ============================================
# 小提琴图（展示关键基因在各细胞类型的表达分布）
# ============================================

key_genes <- c("CD3D", "CD79A", "CD14", "FCGR3A", "NKG7", "CD8A")
vln_plot <- VlnPlot(pbmc, features = key_genes, ncol = 3, pt.size = 0)
ggsave("results/violin_plot.png", plot = vln_plot, width = 12, height = 8, dpi = 300)

# ============================================
# 特征图（UMAP 上标基因表达）
# ============================================

feature_plot <- FeaturePlot(pbmc, features = key_genes, ncol = 3)
ggsave("results/feature_plot.png", plot = feature_plot, width = 15, height = 10, dpi = 300)

# ============================================
# 山脊线图（基因表达密度分布）
# ============================================

ridge_plot <- RidgePlot(pbmc, features = key_genes, ncol = 3)
ggsave("results/ridge_plot.png", plot = ridge_plot, width = 12, height = 10, dpi = 300)

# ============================================
# 热图
# ============================================

top5_markers <- pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

heatmap_plot <- DoHeatmap(pbmc, features = top5_markers$gene) + NoLegend()
ggsave("results/heatmap.png", plot = heatmap_plot, width = 12, height = 10, dpi = 300)

# ============================================
# 火山图（单核细胞亚群比较）
# ============================================

monocyte_markers <- FindMarkers(pbmc, 
                                ident.1 = "CD14+ Mono", 
                                ident.2 = "FCGR3A+ Mono",
                                min.pct = 0.25,
                                logfc.threshold = 0.25)

monocyte_markers$gene <- rownames(monocyte_markers)
monocyte_markers$sig <- ifelse(monocyte_markers$p_val_adj < 0.05 & 
                                 abs(monocyte_markers$avg_log2FC) > 0.5, 
                               "Significant", "Not Significant")

volcano_plot <- ggplot(monocyte_markers, 
                       aes(x = avg_log2FC, 
                           y = -log10(p_val_adj), 
                           color = sig)) +
  geom_point(alpha = 0.6, size = 1) +
  scale_color_manual(values = c("grey", "red")) +
  theme_minimal() +
  labs(title = "CD14+ Mono vs FCGR3A+ Mono",
       x = "Log2 Fold Change",
       y = "-log10(Adjusted P-value)") +
  theme(legend.position = "bottom")

ggsave("results/volcano_plot.png", plot = volcano_plot, width = 8, height = 6, dpi = 300)

# 保存火山图的差异表达结果
write.csv(monocyte_markers, file = "results/monocyte_diff_expression.csv", row.names = FALSE)

# ============================================
# GO 富集分析（B细胞 vs CD8+ T细胞）
# ============================================

message("正在重新生成 GO 富集分析...")

# 加载必要的包
library(Seurat)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)

# 加载数据
pbmc <- readRDS("results/pbmc3k_final.rds")

# 找差异基因
b_vs_t_markers <- FindMarkers(pbmc, 
                              ident.1 = "B cell", 
                              ident.2 = "CD8+ T",
                              min.pct = 0.25,
                              logfc.threshold = 0.5)

# 提取上调基因
up_genes <- b_vs_t_markers %>%
  filter(avg_log2FC > 1 & p_val_adj < 0.05) %>%
  rownames()

# 转换基因名
entrez_ids <- bitr(up_genes, 
                   fromType = "SYMBOL", 
                   toType = "ENTREZID", 
                   OrgDb = org.Hs.eg.db)

# GO 分析
go_result <- enrichGO(gene = entrez_ids$ENTREZID,
                      OrgDb = org.Hs.eg.db,
                      ont = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.10)

# 打印显著富集的通路
cat("\n========== GO 富集分析结果 ==========\n")
cat(paste("显著富集的通路数量:", nrow(go_result), "\n\n"))

if (nrow(go_result) > 0) {
  # 创建结果数据框
  go_df <- as.data.frame(go_result)
  
  # 选择关键列显示
  go_display <- go_df[, c("ID", "Description", "p.adjust", "Count")]
  
  # 按 p 值排序
  go_display <- go_display[order(go_display$p.adjust), ]
  
  # 打印所有通路
  cat("所有显著富集的通路 (按 p 值排序):\n")
  cat("----------------------------------------\n")
  for (i in 1:nrow(go_display)) {
    cat(sprintf("%d. %s\n", i, go_display$Description[i]))
    cat(sprintf("   p.adjust = %.2e, 基因数 = %d\n", 
                go_display$p.adjust[i], go_display$Count[i]))
    cat("\n")
  }

  # 保存完整结果到文件
  write.csv(go_df, file = "results/go_enrichment_full_results.csv", row.names = FALSE)
  cat("\n✅ 完整结果已保存到 results/go_enrichment_full_results.csv\n")
  
} else {
  cat("没有发现显著富集的通路\n")
}

# 保存 GO 结果表格
if (nrow(go_result) > 0) {
  write.csv(as.data.frame(go_result), file = "results/go_enrichment_results.csv", row.names = FALSE)
  
  go_plot <- dotplot(go_result, showCategory = 15) + 
    ggtitle("GO Enrichment: B cells vs CD8+ T cells")
  ggsave("results/go_enrichment.png", plot = go_plot, width = 12, height = 10, dpi = 300)
  
  message("✅ go_enrichment_results.csv 已生成")
} else {
  message("⚠️ GO 富集没有显著结果，已跳过")
}

# ============================================
# 完成
# ============================================

message("\n✅ 所有分析完成！")
message("生成的文件：")
message("  图片: umap_pbmc3k_final.png, violin_plot.png, feature_plot.png, ridge_plot.png, heatmap.png, volcano_plot.png, go_enrichment.png")
message("  表格: marker_genes.csv, marker_genes_top10.csv, monocyte_diff_expression.csv, go_enrichment_results.csv")

