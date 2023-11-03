install.packages("dplyr")
options(repos = c(CRAN = "https://cloud.r-project.org"))
install.packages('biclust')
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("bioDist")

install.packages("withr")
library(biclust)
library(dplyr)
install.packages("pheatmap")
install.packages("gridExtra")
library(pheatmap)
library(gridExtra)
install.packages("infotheo")
library(infotheo)

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
dendrogram <- hclust(dist(genome_matrix, method = "manhattan"))

row_clusters <- cutree(dendrogram, k=17)
col_clusters <- cutree(dendrogram, k=17)

new_data <- genome_matrix[row_clusters, col_clusters]
heatmap <- pheatmap(new_data, clustering_distance_rows = "manhattan",
         clustering_distance_cols = "manhattan",
         cluster_rows = TRUE, cluster_cols = TRUE,
         main = "Hierarchical Clustering Heatmap (Manhattan Distance)")

row_kmeans_clusters <- kmeans(genome_matrix, centers = 17)$cluster
col_kmeans_clusters <- kmeans(t(genome_matrix), centers = 17)$cluster

# Create dataframes for row and column clusters
row_clusters_df <- data.frame(Row = rownames(genome_matrix), RowCluster = row_kmeans_clusters)
col_clusters_df <- data.frame(Column = colnames(genome_matrix), ColCluster = col_kmeans_clusters)

mi_matrix <- matrix(0, ncol(genome_matrix), ncol(genome_matrix))
for (i in 1:(ncol(genome_matrix) - 1)) {
  for (j in (i + 1):ncol(genome_matrix)) {
    mi_value <- mutinformation(genome_matrix[, i], genome_matrix[, j])
    mi_matrix[i, j] <- mi_value
    mi_matrix[j, i] <- mi_value
  }
}

# Convert mutual information to a distance measure (e.g., using 1 - MI)
distance_matrix <- 1 - mi_matrix

dist_object <- as.dist(distance_matrix)
custom_heatmap <- pheatmap(genome_matrix, cluster_rows = TRUE, cluster_cols = TRUE, scale='none', main = 'Original heatmap', kmeans_k = 17, clustering_distance_rows = dist_object, clustering_distance_cols = dist_object)

