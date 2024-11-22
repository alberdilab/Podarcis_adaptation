#Data integration-pca with functional microbiota data and respirometry data####


#GIFT ELEMENTS

##Include sample ID that appears also in respirometry data
GIFTs_elements_community_pca<-as.data.frame(GIFTs_elements_community)

GIFTs_elements_community_pca <- GIFTs_elements_community_pca %>%
  tibble::rownames_to_column(var = "Tube_code")

#merge gift data with the sample_metadata data frame
gifts_pca<-GIFTs_elements_community_pca %>%
  left_join(sample_metadata, by =join_by(Tube_code==Tube_code)) %>%
  filter(time_point=="1_Acclimation"|time_point=="5_Post-FMT1"|time_point=="6_Post-FMT2")

#change time_point factors to 0,1 and 2
gifts_pca<-gifts_pca %>%
  mutate(
  time_point = fct_recode(time_point, `0` = "1_Acclimation",
                          `1` = "5_Post-FMT1",
                          `2` = "6_Post-FMT2"))


#create interaction variable with individual and time-point
gifts_pca$newID <- paste(gifts_pca$individual, "_", gifts_pca$time_point)
respirometry_resp2_subset$newID<-paste(respirometry_resp2_subset$individual, "_", respirometry_resp2_subset$time_point)

#exclude non numeric columns
gifts_pca_filter<-gifts_pca %>%
  dplyr::select(-c(1,151:161))

# Perform PCA
pca_results <- prcomp(gifts_pca_filter, scale. = TRUE)
pca_scores <- as.data.frame(pca_results$x)
pca_scores$newID <- gifts_pca$newID
gift_resp2 <- inner_join(pca_scores, respirometry_resp2_subset, by = "newID")

ggplot(gift_resp2, aes(x = PC1, y = PC2, color = type, shape=time_point)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "PCA of Microbiota Data",
       x = "PC1",
       y = "PC2",
       color = "type",
       shape="Time_point")
cor.test(gift_resp2$PC1, gift_resp2$QC_normalized, method = "pearson")
cor.test(gift_resp2$PC1, gift_resp2$QC_normalized, method = "spearman")


#GIFT FUNCTIONS

##Include sample ID that appears also in respirometry data
GIFTs_functions_community_pca<-as.data.frame(GIFTs_functions_community)

GIFTs_functions_community_pca <- GIFTs_functions_community_pca %>%
  tibble::rownames_to_column(var = "Tube_code")

#merge gift data with the sample_metadata data frame
gifts_func_pca<-GIFTs_functions_community_pca %>%
  left_join(sample_metadata, by =join_by(Tube_code==Tube_code)) %>%
  filter(time_point=="1_Acclimation"|time_point=="5_Post-FMT1"|time_point=="6_Post-FMT2")

#change time_point factors to 0,1 and 2
gifts_func_pca<-gifts_func_pca %>%
  mutate(
    time_point = fct_recode(time_point, `0` = "1_Acclimation",
                            `1` = "5_Post-FMT1",
                            `2` = "6_Post-FMT2"))


#create interaction variable with individual and time-point
gifts_func_pca$newID <- paste(gifts_func_pca$individual, "_", gifts_func_pca$time_point)
respirometry_resp2_subset$newID<-paste(respirometry_resp2_subset$individual, "_", respirometry_resp2_subset$time_point)

#exclude non numeric columns
gifts_func_pca_filter<-gifts_func_pca %>%
  dplyr::select(-c(1,22:32))

# Perform PCA
pca_func_results <- prcomp(gifts_func_pca_filter, scale. = TRUE)
pca_func_scores <- as.data.frame(pca_func_results$x)
pca_func_scores$newID <- gifts_func_pca$newID
gift_resp2_func <- inner_join(pca_func_scores, respirometry_resp2_subset, by = "newID")

ggplot(gift_resp2_func, aes(x = PC1, y = PC2, color = type, shape=time_point)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "PCA of Microbiota Data",
       x = "PC1",
       y = "PC2",
       color = "type",
       shape="Time_point")
cor.test(gift_resp2_func$PC1, gift_resp2_func$QC_normalized, method = "pearson")
cor.test(gift_resp2_func$PC1, gift_resp2_func$QC_normalized, method = "spearman")

#Data integration-pca with beta div microbiota data and respirometry data####

## Acclimation samples
# Perform classical MDS (PCA) on the distance matrix
pca_beta_div_neutral_accli <- cmdscale(beta_div_neutral_accli$S, k = 2)  # k = 2 for 2 principal components

# Convert to a data frame for easy plotting
pca_beta_div_neutral_accli_df <- as.data.frame(pca_beta_div_neutral_accli)
colnames(pca_beta_div_neutral_accli_df) <- c("PC1", "PC2")

#subset respirometry data to only acclimation measurements
respirometry_resp2_subset_accli<-respirometry_resp2_subset %>%
  filter(time_point=="0")

# Assuming 'respirometry_data' has a 'SampleID' column
pca_beta_div_neutral_accli_df <- pca_beta_div_neutral_accli_df %>%
  tibble::rownames_to_column(var = "Tube_code")
pca_beta_div_neutral_accli_df<-inner_join(pca_beta_div_neutral_accli_df, sample_metadata, by="Tube_code")
resp2_beta_neutral_accli <- inner_join(pca_beta_div_neutral_accli_df, respirometry_resp2_subset_accli, by="individual")

# Plot the PCA with respirometry variable as color
ggplot(resp2_beta_neutral_accli, aes(x = PC1, y = PC2, color = type.y)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "PCA of Beta Diversity",
       x = "PC1",
       y = "PC2",
       color = "Type")

cor.test(resp2_beta_neutral_accli$PC1, resp2_beta_neutral_accli$QC_normalized, method = "pearson")
cor.test(resp2_beta_neutral_accli$PC1, resp2_beta_neutral_accli$QC_normalized, method = "spearman")

pca_result_neutral_accli <- cmdscale(distance_matrix, eig = TRUE, k = 2)

# Calculate the percentage of variance explained
eig_values <- pca_result$eig
var_explained <- eig_values / sum(eig_values) * 100
pc1_var <- round(var_explained[1], 2)
pc2_var <- round(var_explained[2], 2)



## Post1 samples
# Perform classical MDS (PCA) on the distance matrix
pca_beta_div_neutral_post1 <- cmdscale(beta_div_neutral_post1$S, k = 2)  # k = 2 for 2 principal components

# Convert to a data frame for easy plotting
pca_beta_div_neutral_post1_df <- as.data.frame(pca_beta_div_neutral_post1)
colnames(pca_beta_div_neutral_post1_df) <- c("PC1", "PC2")

#subset respirometry data to only acclimation measurements
respirometry_resp2_subset_post1<-respirometry_resp2_subset %>%
  filter(time_point=="1")

# Assuming 'respirometry_data' has a 'SampleID' column
pca_beta_div_neutral_post1_df <- pca_beta_div_neutral_post1_df %>%
  tibble::rownames_to_column(var = "Tube_code")
pca_beta_div_neutral_post1_df<-inner_join(pca_beta_div_neutral_post1_df, sample_metadata, by="Tube_code")
resp2_beta_neutral_post1 <- inner_join(pca_beta_div_neutral_post1_df, respirometry_resp2_subset_post1, by="individual")

# Plot the PCA with respirometry variable as color
ggplot(resp2_beta_neutral_post1, aes(x = PC1, y = PC2, color = type.y)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "PCA of Beta Diversity",
       x = "PC1",
       y = "PC2",
       color = "Type")

cor.test(resp2_beta_neutral_post1$PC1, resp2_beta_neutral_post1$QC_normalized, method = "pearson")
cor.test(resp2_beta_neutral_post1$PC1, resp2_beta_neutral_post1$QC_normalized, method = "spearman")

## Post2 samples
# Perform classical MDS (PCA) on the distance matrix
pca_beta_div_neutral_post2 <- cmdscale(beta_div_neutral_post2$S, k = 2)  # k = 2 for 2 principal components

# Convert to a data frame for easy plotting
pca_beta_div_neutral_post2_df <- as.data.frame(pca_beta_div_neutral_post2)
colnames(pca_beta_div_neutral_post2_df) <- c("PC1", "PC2")

#subset respirometry data to only acclimation measurements
respirometry_resp2_subset_post2<-respirometry_resp2_subset %>%
  filter(time_point=="2")

# Assuming 'respirometry_data' has a 'SampleID' column
pca_beta_div_neutral_post2_df <- pca_beta_div_neutral_post2_df %>%
  tibble::rownames_to_column(var = "Tube_code")
pca_beta_div_neutral_post2_df<-inner_join(pca_beta_div_neutral_post2_df, sample_metadata, by="Tube_code")
resp2_beta_neutral_post2 <- inner_join(pca_beta_div_neutral_post2_df, respirometry_resp2_subset_post2, by="individual")

# Plot the PCA with respirometry variable as color
ggplot(resp2_beta_neutral_post2, aes(x = PC1, y = PC2, color = type.y)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "PCA of Beta Diversity",
       x = "PC1",
       y = "PC2",
       color = "Type")

cor.test(resp2_beta_neutral_post2$PC1, resp2_beta_neutral_post2$QC_normalized, method = "pearson")
cor.test(resp2_beta_neutral_post2$PC1, resp2_beta_neutral_post2$QC_normalized, method = "spearman")

###Functional
# Perform classical MDS (PCA) on the distance matrix
pca_beta_div_func_post2 <- cmdscale(beta_div_func_post2$S, k = 2)  # k = 2 for 2 principal components

# Convert to a data frame for easy plotting
pca_beta_div_func_post2_df <- as.data.frame(pca_beta_div_func_post2)
colnames(pca_beta_div_func_post2_df) <- c("PC1", "PC2")

#subset respirometry data to only acclimation measurements
respirometry_resp2_subset_post2<-respirometry_resp2_subset %>%
  filter(time_point=="2")

# Assuming 'respirometry_data' has a 'SampleID' column
pca_beta_div_func_post2_df <- pca_beta_div_func_post2_df %>%
  tibble::rownames_to_column(var = "Tube_code")
pca_beta_div_func_post2_df<-inner_join(pca_beta_div_func_post2_df, sample_metadata, by="Tube_code")
resp2_beta_func_post2 <- inner_join(pca_beta_div_func_post2_df, respirometry_resp2_subset_post2, by="individual")

# Plot the PCA with respirometry variable as color
ggplot(resp2_beta_func_post2, aes(x = PC1, y = PC2, color = type.y)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "PCA of Beta Diversity",
       x = "PC1",# Perform PCA (MDS) on the Distance Matrix,
       y = "PC2",
       color = "Type")

cor.test(resp2_beta_func_post2$PC1, resp2_beta_func_post2$QC_normalized, method = "pearson")
cor.test(resp2_beta_func_post2$PC1, resp2_beta_func_post2$QC_normalized, method = "spearman")
