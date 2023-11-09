install.packages("dplyr")
options(repos = c(CRAN = "https://cloud.r-project.org"))
install.packages('biclust')
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("bioDist")

install.packages('dendextend')
library(dendextend)
install.packages('ape')
library(ape)

install.packages(c("htmlwidgets", "png", "base64enc"))

BiocManager::install('ComplexHeatmap')
BiocManager::install('InteractiveComplexHeatmap')
BiocManager::install('cola')

install.packages("Cairo")
library(Cairo)
library(cola)
library(ComplexHeatmap)
library(InteractiveComplexHeatmap)
install.packages('/Users/gautham/microbiome_analysis/d3heatmap_0.6.0.tar.gz', repos = NULL, type = 'source')
library(d3heatmap)
install.packages('jaccard_heatmap')
install.packages('dynamicTreeCut')
library(dynamicTreeCut)
install.packages("withr")
library(biclust)
library(dplyr)
install.packages("pheatmap")
install.packages("gridExtra")
library(pheatmap)
library(gridExtra)
install.packages("infotheo")
library(infotheo)
install.packages('proxy')
library(proxy)
install.packages("entropy")
library(entropy)

install.packages("Information")
library(Information)


df <- test_array_output
selected_columns <- df %>%
  select(contains("GCF"))

tcid_acc_cols <- df[c("TCID", "Acc")]
tc_list <- list()

for (i in 1:nrow(tcid_acc_cols)){
  
  tcid <- tcid_acc_cols[i, 1]
  acc <- tcid_acc_cols[i, 2]
  
  tc_acc <- paste(tcid, acc, sep = "-")
  tc_list[[i]] <- tc_acc
}



genome_matrix_df <- selected_columns %>% 
  mutate_all(~as.numeric(ifelse(. == '+', 1, ifelse(. == '-', 0, .))))


rownames(genome_matrix_df) <- tc_list
genome_matrix <- as.matrix(genome_matrix_df)

res <- biclust(x=genome_matrix, method=BCBimax(), minr=4, minc=4, number=17)

# this is using manhattan distance measurement as default for both rows and columns using kmeans clustering
#original_heatmap <- pheatmap(genome_matrix, cluster_rows = TRUE, cluster_cols = TRUE, scale='none', main = 'Original heatmap', kmeans_k = 17, clustering_distance_rows = 'manhattan', clustering_distance_cols = 'manhattan')
# temp <- proxy::dist(genome_matrix, method = 'Manhattan')
# hc <- hclust(temp, method = 'complete')
# plot(hc)
# dendrogram <- hclust(dist(genome_matrix, method = "manhattan"))
# plot(dendrogram)
# 
# row_clusters <- cutree(dendrogram, k=17)
# col_clusters <- cutree(dendrogram, k=17)
# 
# new_data <- genome_matrix[row_clusters, col_clusters]
# heatmap <- pheatmap(new_data, clustering_distance_rows = "manhattan",
#          clustering_distance_cols = "manhattan",
#          cluster_rows = TRUE, cluster_cols = TRUE,
#          main = "Hierarchical Clustering Heatmap (Manhattan Distance)")

original_heatmap <- pheatmap(genome_matrix, cluster_rows = TRUE, cluster_cols = TRUE, scale='none', main = 'Original heatmap', clustering_method = "complete", clustering_distance_rows = 'manhattan', clustering_distance_cols = 'manhattan')

row_dendogran <- original_heatmap$tree_row
row_groups <- cutree(row_dendogran, h = 1)
plot(row_dendogran)
rect.hclust(row_dendogran, k=17, border = 1:17)
protein_clusters = cutree(row_dendogran, k=17)

col_dendogram <- original_heatmap$tree_col
col_groups <- cutree(col_dendogram, h = 1)
plot(col_dendogram)
rect.hclust(col_dendogram, k=6, border = 1:17)
genome_clusters = cutree(col_dendogram, k=17)

# row_kmeans_clusters <- kmeans(genome_matrix, centers = 17)$cluster
# col_kmeans_clusters <- kmeans(t(genome_matrix), centers = 17)$cluster

# Create dataframes for row and column clusters
# row_clusters_df <- data.frame(Row = rownames(genome_matrix), RowCluster = row_kmeans_clusters)
# col_clusters_df <- data.frame(Column = colnames(genome_matrix), ColCluster = col_kmeans_clusters)

mi_matrix <- matrix(0, ncol(genome_matrix), ncol(genome_matrix))
# for (i in 1:(ncol(genome_matrix) - 1)) {
#   for (j in (i + 1):ncol(genome_matrix)) {
#     mi_value <- mutinformation(genome_matrix[, i], genome_matrix[, j])
#     mi_matrix[i, j] <- mi_value
#     mi_matrix[j, i] <- mi_value
#   }
# }
transpose_genome_matrix <- as.data.frame(t(genome_matrix_df))
                            
mi_matrix <- mutinformation(transpose_genome_matrix, method='emp') # big runtime
mi_matrix_genomes <- mutinformation(genome_matrix_df, method='emp')

# Convert mutual information to a distance measure (e.g., using 1 - MI)
distance_matrix_protein <- 1 - mi_matrix
distance_matrix_genome <- 1 - mi_matrix_genomes

dist_object_row <- as.dist(distance_matrix_protein)
dist_object_col <- as.dist(distance_matrix_genome)

custom_heatmap <- pheatmap(
  genome_matrix,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  scale = 'none',
  main = 'Using mutual information',
  clustering_method = "complete",
  clustering_distance_rows = dist_object_row,
  clustering_distance_cols = dist_object_col
)


d3heatmap(scale(genome_matrix), colors = 'RdYlBu', k_row = 17, k_col = 6)

test <- Heatmap(genome_matrix,
                   name = "test",
                   clustering_distance_rows = dist_object_row,
                   clustering_distance_columns = dist_object_col)
test <- draw(test)
htShiny(test)
# # test <- Heatmap(genome_matrix)
# # test <- draw(test)
# 
#htShiny(test)

genome_matrix.clust <- cbind(genome_matrix, cluster = cutree(custom_heatmap$tree_row, k = 17))

for (i in 1:17) {
  subset_data <- genome_matrix.clust[genome_matrix.clust[,"cluster"] == i, -ncol(genome_matrix.clust)]
  pheatmap(
    subset_data,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    scale = 'none',
    main = paste("Cluster", i),
    clustering_method = "complete",
    clustering_distance_rows = dist_object_row,
    clustering_distance_cols = dist_object_col
  )
  
}

my_tree <- as.phylo(custom_heatmap$tree_row)

write.tree(phy=my_tree, file="file_name.tree")  
cd <- cutreeDynamic(dendro = custom_heatmap$tree_row, minClusterSize = 20)

row_dendrogram <- custom_heatmap$tree_row
dend <- as.dendrogram(row_dendrogram)

for (i in seq_along(cd)) {
  # Color the branches of the dendrogram based on the current cluster assignment
  dend_colored <- color_branches(dend, k = cd[[i]])
  
  # Plot the dendrogram with colored clusters
  plot(dend_colored, main = paste("Cluster", i))
  break
}



row_dendogran <- custom_heatmap$tree_row
row_groups <- cutree(row_dendogran, h = 0.9)
plot(row_dendogran)
rect.hclust(row_dendogran, k=6, border = 1:17)
protein_clusters = cutree(row_dendogran, k=17)

col_dendogram <- custom_heatmap$tree_col
col_groups <- cutree(col_dendogram, h = 0.9)
plot(col_dendogram)
rect.hclust(col_dendogram, k=6, border = 1:17)
genome_clusters = cutree(col_dendogram, k=6)

# cluster_heatmaps <- list()
# 
# for (i in 1:max(protein_clusters)) {
#   for (j in 1:max(genome_clusters)) {
#     # Identify the rows and columns belonging to the current cluster
#     rows_in_cluster <- which(protein_clusters == i)
#     cols_in_cluster <- which(genome_clusters == j)
#     
#     # Subset the data for the current cluster
#     subset_data <- genome_matrix[rows_in_cluster, cols_in_cluster]
#     
#     # Create a heatmap for the current cluster
#     heatmap <- pheatmap(
#       subset_data,
#       scale = 'none',
#       main = paste("Row Cluster", i, " - ", "Column Cluster", j),
#       clustering_method = "complete"
#     )
#     
#     row_dendogran <- heatmap$tree_row
#     row_groups <- cutree(row_dendogran, h = 1)
#     plot(row_dendogran)
#     rect.hclust(row_dendogran, k=6, border = 1:17)
#     temp_p = cutree(row_dendogran, k=6)
#     
#     col_dendogram <- heatmap$tree_col
#     col_groups <- cutree(col_dendogram, h = 4)
#     plot(col_dendogram)
#     rect.hclust(col_dendogram, k=2, border = 1:17)
#     temp_g = cutree(col_dendogram, k=6)
#     # Store the heatmap in the list
#     cluster_heatmaps[[paste("RowCluster_", i, "_ColumnCluster_", j)]] <- heatmap
#   }
# }

cluster_heatmaps <- list()

# Iterate through the clusters and create individual heatmaps
for (i in 1:max(protein_clusters)) {
  for (j in 1:max(genome_clusters)) {
    # Subset the data for the current cluster
    subset_data <- genome_matrix[row_groups == i & col_groups == j, ]
    
    if (is.null(subset_data) || length(subset_data) == 0) {
      cat("subset_data is empty or contains no data.\n")
    } 
    else {
      heatmap <- pheatmap(
        subset_data,
        cluster_rows = FALSE,  # Avoid re-clustering within the cluster
        cluster_cols = FALSE,  # Avoid re-clustering within the cluster
        scale = 'none',
        main = paste("Row Cluster", i, " - ", "Column Cluster", j),
        clustering_method = "none"  # No clustering within the cluster
      )
      
      # Store the heatmap in the list
      cluster_heatmaps[[paste("RowCluster_", i, "_ColumnCluster_", j)]] <- heatmap
      
    }
    # Create a heatmap for the current cluster
      
    }
}

# test <- biclust(x=genome_matrix, method=BCBimax(), minr=4, minc=4, number=17, dist = list(row = dist_object_row, col = dist_object_col))
