install.packages("dplyr")
install.packages("isa2")
options(repos = c(CRAN = "https://cloud.r-project.org"))
library(isa2)
install.packages('biclust')

install.packages("withr")
library(biclust)
library(dplyr)


df <- test_array_output
selected_columns <- df %>%
  select(contains("GCF"))

genome_matrix_df <- selected_columns %>% 
  mutate_all(~ifelse(. == '+', 1, ifelse(. == '-', 0, .)))

genome_matrix <- matrix(genome_matrix_df)
# genome_matrix <- matrix(genome_matrix)
genome_numeric <- as.numeric(genome_matrix)

temp <- biclust(genome_matrix_df, method = 'BCPlaid')

cluster <- biclust(genome_matrix, BCPlaid(), back.fit = 2, shuffle = 3, fit.model = ~m
                   + a + b, iter.startup = 5, iter.layer = 30, verbose = TRUE)
