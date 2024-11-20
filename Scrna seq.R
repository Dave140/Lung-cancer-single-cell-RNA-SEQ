# Load Required Libraries
library(Seurat)        # For single-cell data analysis
library(SingleR)       # For cell annotation
library(DoubletFinder) # For doublet detection
library(dplyr)         # For data manipulation
library(ggplot2)       # For visualization
library(enrichR)       # For enrichment analysis
library(Harmony)       # For batch correction (optional)
library(CellTypist)    # For cell-type annotation (optional)

# Set Parameters
variable_genes_list <- c(2000, 1000, 500, 100, 50, 10) # Different numbers of top variable genes to test
pca_dims_list <- c(10, 25)                            # Different PCA dimensions to test
resolution_list <- c(0.4, 0.6, 0.8, 1.0)              # Different clustering resolutions to test

# Step 1: Data Loading
# Replace with your dataset path
sc_data <- Read10X(data.dir = "path_to_your_10X_data/")
seurat_obj <- CreateSeuratObject(counts = sc_data, project = "scRNAseq", min.cells = 3, min.features = 200)

# Step 2: Quality Control
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Step 3: Doublet Detection with DoubletFinder
# Estimate the expected number of doublets (adjust this based on your data)
seurat_obj <- doubletFinder_v3(seurat_obj, PCs = 1:20, pN = 0.25, pK = 0.05, nExp = 100, reuse.pANN = FALSE)

# Add doublet identification to metadata (e.g., "Doublet" column)
seurat_obj$doublet_status <- seurat_obj$DF.classifications_0.25_0.05_100

# Step 4: Remove Doublets
# Filter out doublets (doublet cells identified as "Doublet" in the metadata)
seurat_obj <- subset(seurat_obj, subset = doublet_status == "Singlet")

# Step 5: Re-normalization
# Re-normalize the data after removing doublets to ensure proper scaling
seurat_obj <- SCTransform(seurat_obj, verbose = FALSE)

# Step 6: Batch Correction (if needed)
# If you have multiple batches (e.g., different samples), you can use Harmony for batch correction
# seurat_obj <- RunHarmony(seurat_obj, group.by.vars = "batch")

# Step 7: Loop Through Different Variable Genes, PCA Dimensions, and Clustering Resolutions
# Initialize a list to store results
results <- list()

for (top_variable_genes in variable_genes_list) {
  for (pca_dims in pca_dims_list) {
    for (resolution in resolution_list) {
      cat("Processing with", top_variable_genes, "variable genes,", pca_dims, "PCA dimensions, and", resolution, "resolution...\n")
      
      # Normalization and Variable Feature Selection
      seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = top_variable_genes)
      
      # Scaling and PCA
      seurat_obj <- ScaleData(seurat_obj, features = VariableFeatures(seurat_obj))
      seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj), npcs = max(pca_dims_list))
      
      # Elbow Plot for PCA
      ElbowPlot(seurat_obj, ndims = max(pca_dims_list))
      
      # Clustering and UMAP
      seurat_obj <- FindNeighbors(seurat_obj, dims = 1:pca_dims)
      seurat_obj <- FindClusters(seurat_obj, resolution = resolution)
      seurat_obj <- RunUMAP(seurat_obj, dims = 1:pca_dims)
      
      # Identify Consistently Expressed Genes
      consistently_expressed_genes <- rownames(seurat_obj@assays$RNA@counts)[
        rowMeans(seurat_obj@assays$RNA@counts > 0) == 1
      ]
      
      # Save UMAP Plot for Visualization
      plot_title <- paste0("UMAP (Top ", top_variable_genes, " genes, ", pca_dims, " PCA dims, Res", resolution, ")")
      DimPlot(seurat_obj, reduction = "umap", label = TRUE) + ggtitle(plot_title)
      
      # Differential Expression for Clusters
      markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
      selected_markers <- markers %>%
        filter(gene %in% consistently_expressed_genes) %>%
        arrange(desc(avg_log2FC))
      
      # Save Results
      results[[paste0("VarGenes_", top_variable_genes, "_PCADims_", pca_dims, "_Res_", resolution)]] <- list(
        "UMAP_Plot" = plot_title,
        "Markers" = selected_markers
      )
      
      # Print Top Markers
      cat("Top markers for", top_variable_genes, "variable genes,", pca_dims, "PCA dimensions, and resolution", resolution, ":\n")
      print(head(selected_markers))
    }
  }
}

# Step 8: Correlation Test
# Example: Correlation between gene expression values of two specific genes
# Here I select two genes as an example, replace with genes of interest

gene_1 <- "TP53"  # Replace with gene of interest
gene_2 <- "KRAS"  # Replace with gene of interest

# Extract expression values for the two genes
gene_1_expr <- FetchData(seurat_obj, vars = gene_1)
gene_2_expr <- FetchData(seurat_obj, vars = gene_2)

# Pearson Correlation
pearson_corr <- cor(gene_1_expr, gene_2_expr, method = "pearson")

# Spearman Correlation
spearman_corr <- cor(gene_1_expr, gene_2_expr, method = "spearman")

# Print Correlation Results
cat("Pearson Correlation between", gene_1, "and", gene_2, ":", pearson_corr, "\n")
cat("Spearman Correlation between", gene_1, "and", gene_2, ":", spearman_corr, "\n")

# Optional: Correlation between cluster identities and gene expression
# Get cluster identities
seurat_obj$cluster <- seurat_obj$seurat_clusters

# Perform correlation for each cluster
cluster_gene_correlation <- list()
for (cluster_id in unique(seurat_obj$cluster)) {
  cluster_cells <- WhichCells(seurat_obj, idents = cluster_id)
  cluster_expr <- FetchData(seurat_obj, vars = c(gene_1, gene_2), cells = cluster_cells)
  pearson_corr_cluster <- cor(cluster_expr[, gene_1], cluster_expr[, gene_2], method = "pearson")
  spearman_corr_cluster <- cor(cluster_expr[, gene_1], cluster_expr[, gene_2], method = "spearman")
  
  cluster_gene_correlation[[paste("Cluster", cluster_id)]] <- list(
    pearson = pearson_corr_cluster,
    spearman = spearman_corr_cluster
  )
}

# Print correlation for each cluster
for (cluster_id in names(cluster_gene_correlation)) {
  cat(cluster_id, ": Pearson =", cluster_gene_correlation[[cluster_id]]$pearson, 
      ", Spearman =", cluster_gene_correlation[[cluster_id]]$spearman, "\n")
  
  # Step 9: Cell Annotation with SingleR
  # Use the best-performing configuration based on the consistency of markers
  ref <- celldex::HumanPrimaryCellAtlasData()
  pred <- SingleR(test = GetAssayData(seurat_obj, slot = "data"), ref = ref, labels = ref$label.main)
  seurat_obj$SingleR <- pred$labels
  
  # Step 10: Enrichment Analysis for Selected Marker Genes
  # Perform enrichment analysis for selected marker genes
  enrich_results <- enrichr(
    unique(unlist(lapply(results, function(x) x$Markers$gene))),
    databases = c("GO_Biological_Process", "KEGG_2021_Human")
  )
  print(enrich_results)
  
  # Step 11: Cell Cycle Scoring
  # Perform cell cycle scoring (optional)
  seurat_obj <- CellCycleScoring(seurat_obj, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = TRUE)
  DimPlot(seurat_obj, group.by = "Phase", reduction = "umap")
  
  # Step 12: Save Final Object and Results
  saveRDS(seurat_obj, file = "seurat_object_final.rds")
  write.csv(enrich_results, "enrichment_results.csv")
  