############################################################
# Script: Full MOFA Multi-Omics Analysis Workflow
# Author: [Destras et al]
# Date: [2026-03-23]
#
# Modules:
#   1. Prepare input tables for MOFA
#   2. Analyze MOFA model (variance, factors)
#   3. Perform UMAP and K-means clustering (Elbow, Silhouette, Gap)
#   4. Identify specific features associated with clusters
#   4. Run differential expression per cluster
#   5. Conduct GO enrichment and visualize with dot plot
#
# All variables, paths, and data frame names standardized.
############################################################


############################################################
# 0. Load required libraries
############################################################
library(MOFA2)
library(dplyr)
library(tidyr)
library(reshape2)
library(data.table)
library(readr)
library(ggplot2)
library(ggpubr)
library(cluster)
library(factoextra)
library(edgeR)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(magrittr)
library(effsize)

############################################################
# 1. Preparation of tables for MOFA
############################################################

# --- Paths ---
path_root   <- "Verdy"
path_output <- file.path(path_root, "path_output")
dir.create(path_output, showWarnings = FALSE)

# --- Read omics tables counts---
bac_counts     <- read.table("Verdy_16S/normalized_table_16S.csv", sep=";", header=TRUE, row.names=1, dec=",")
RNAbac_counts  <- read.table("Verdy_RNAbac/normalized_table_RNAbac.csv", sep=";", header=TRUE, row.names=1, dec=",")
vir_counts     <- read.table("Verdy_virome/normalized_table_virome.csv", sep=";", header=TRUE, row.names=1, dec=",")

# --- Read clinical metadata (separate temporary clinical metadata files initially created for  mono-omic analyses) ---
meta_bac    <- read.table("Verdy_16S/clinical_temp_vir.txt", sep=" ", header=TRUE, row.names=1)
meta_RNAbac <- read.table("Verdy_metabac/clinical_temp_bacrep.txt", sep=" ", header=TRUE, row.names=1)
meta_vir    <- read.table("Verdy_virome/clinical_temp_otutable.txt", sep="\t", header=TRUE, row.names=1)

# --- Select relevant columns of abundance and diversity indexes (added to separate clinical metadata files herein)---
meta_bac <- meta_bac %>% select(copies_total, index_shannon, index_richness)
meta_RNAbac <- meta_RNAbac %>% select(RPM_total, index_shannon, index_richness)
meta_vir <- meta_vir %>% select(VDR_total, index_shannon, index_richness, virulent, temperate, Undetermined)

# --- Rename columns (avoid duplicates) ---
colnames(meta_bac)    <- c("copies_total_bac", "index_shannon_bac", "index_richness_bac")
colnames(meta_RNAbac) <- c("RPM_total_RNAbac", "index_shannon_RNAbac", "index_richness_RNAbac")
colnames(meta_vir)    <- c("VDR_total_vir", "index_shannon_vir", "index_richness_vir","virulent","temperate","Undetermined")

# --- Merge all count tables ---
bac_counts$sample <- row.names(bac_counts)
RNAbac_counts$sample <- row.names(RNAbac_counts)
vir_counts$sample <- row.names(vir_counts)
meta_bac$sample <- row.names(meta_bac)
meta_RNAbac$sample <- row.names(meta_RNAbac)
meta_vir$sample <- row.names(meta_vir)

merged_counts <- bac_counts %>%
  inner_join(RNAbac_counts, by="sample") %>%
  inner_join(vir_counts, by="sample")

merged_meta <- meta_bac %>%
  inner_join(meta_RNAbac, by="sample") %>%
  inner_join(meta_vir, by="sample")

# --- Filter features present in ≥25% samples ---
merged_counts_filt <- merged_counts[, apply(merged_counts, 2, function(x)
  sum(as.numeric(x) > 0) >= 0.25 * nrow(merged_counts))]

# --- Merge counts + clinical ---
mofa_input <- merge(merged_counts_filt, merged_meta, by="sample")
write.table(merged_meta, file.path(path_output, "merged_meta.tsv"),
            sep="\t", row.names=FALSE)

# --- Standardize features with z scores ---
mofa_scaled <- as.data.frame(lapply(mofa_input[,-1], scale, center=TRUE, scale=TRUE))
row.names(mofa_scaled) <- mofa_input$sample

# --- Transpose and tidy for MOFA ---
mofa_transposed <- data.frame(t(mofa_scaled))
mofa_long <- melt(mofa_transposed, id="sample")
colnames(mofa_long) <- c("feature", "sample", "value")

write.table(mofa_long, file.path(path_output, "MOFA_input_tidy.tsv"),
            sep="\t", row.names=FALSE)


############################################################
# 2. MOFA Model Analysis
############################################################
##A - create MOFA in a linux environment##

library(MOFA2)


dt = read.delim(file.path(path_output,"MOFA_input_tidy.tsv"))


####create mofa object

obj  <- create_mofa(data = dt)

####prepare mofa object for training

data_opts <- get_default_data_options(obj)
head(data_opts)


model_opts <- get_default_model_options(obj)
model_opts$num_factors <- 12 # this should be higher once you add more features in your expression data
head(model_opts)

train_opts <- get_default_training_options(obj)
train_opts$convergence_mode <- "fast"
train_opts$seed <- 42



obj <- prepare_mofa(obj, model_options = model_opts,
                    #mefisto_options = mefisto_opts,
                    training_options = train_opts,
                    data_options = data_opts)

obj <- run_mofa(obj, use_basilisk = F,
                outfile = file.path(path_root, "model_new.hdf5"))


##B - Visualisation of MOFA model (this part does not require a linux environment)##

mofa_model <- load_model(file.path(path_root, "model_new.hdf5"))
samples_metadata(mofa_model) <- merged_meta

# --- Variance explained ---
var_info <- get_variance_explained(mofa_model)
write.table(var_info$r2_per_factor, file.path(path_output, "variance_explained.tsv"), sep="\t", quote=FALSE)

# --- Plot variance explained ---
plot_variance_explained(mofa_model, plot_total=TRUE)[[2]]
plot_variance_explained(mofa_model)

# --- Factors extraction ---
factors_df <- get_factors(mofa_model, as.data.frame=TRUE)
factors_wide <- factors_df %>% spread(key=factor, value=value)
write.table(t(factors_wide), file.path(path_output, "factors_matrix.tsv"), sep="\t", quote=FALSE)


############################################################
# 3. UMAP and Clustering
############################################################

mofa_model <- load_model(file.path(path_root, "model_new.hdf5"))
umap_res <- run_umap(mofa_model, n_neighbors=5, min_dist=0, n_components=2, verbose=FALSE)
umap_df <- as.data.frame(umap_res@dim_red[["UMAP"]])
umap_df$sample <- row.names(umap_df)

# --- Determine optimal cluster count ---
fviz_nbclust(umap_df[, c("UMAP1","UMAP2")], kmeans, method="wss") + ggtitle("Elbow Method")
fviz_nbclust(umap_df[, c("UMAP1","UMAP2")], kmeans, method="silhouette") + ggtitle("Silhouette Method")

# Gap statistic (bootstrap = 500)
set.seed(123)
gap_stat <- clusGap(umap_df[, c("UMAP1","UMAP2")], FUN=kmeans, nstart=25, K.max=10, B=500)
fviz_gap_stat(gap_stat) + ggtitle("Gap Statistic")

# --- Perform K-means clustering ---
set.seed(123)
k <- 4
kmeans_res <- kmeans(umap_df[,c("UMAP1","UMAP2")], centers=k, nstart=10)
cluster_df <- data.frame(sample=umap_df$sample, cluster=factor(kmeans_res$cluster))
write.csv(cluster_df, file.path(path_output, "cluster_assignments.csv"), row.names=FALSE)

# --- Merge with metadata ---
umap_annot <- umap_df %>%
  left_join(cluster_df, by="sample") %>%
  left_join(merged_meta, by="sample")

# --- Plot clusters ---
ggplot(umap_annot, aes(UMAP1, UMAP2, color=cluster, fill=cluster)) +
  stat_ellipse(geom="polygon", alpha=0.3) +
  geom_point(size=3) +
  scale_color_brewer(palette="Dark2") +
  scale_fill_brewer(palette="Dark2") +
  theme_bw(base_size=16) +
  labs(title="UMAP of MOFA Factors", x="UMAP1", y="UMAP2")
ggsave(file.path(path_output, "UMAP_clusters.pdf"), width=7, height=6)


############################################################
# 4. Differential Expression of features per Cluster
############################################################
### Determine the qualitative and continuous variables of interest ###

merged_meta<-read.table(file.path(path_output,"merged_meta,tsv",sep="\t"))

qualitative_vars <- c("BPD","Oxygen_requirement", "Sex", "CAN", "Delivrance_mode")

continuous_vars <- c("GA", "BW", "min", "VDR_total_vir", "lg_bacterial_load", 
                     "RPM_total_RNAbac", "Caudoviricetes",
                     "index_shannon_vir", "index_shannon_bac", "index_shannon_RNAbac",
                     "index_richness_vir", "index_richness_bac", "index_richness_RNAbac",
                     "Pseudomonas", "Staphylococcus","temperate","virulent","Undetermined","Eukaryotic_viruses","Bacteriophages")


compute_all_features <- function(data, cluster_col, continuous_vars, qualitative_vars) {
  
  clusters <- unique(data[[cluster_col]])
  results_list <- list()
  
  ### =========================
  ### 1. CONTINUOUS VARIABLES
  ### =========================
  
  for (var in continuous_vars) {
    
    for (cl in clusters) {
      
      data_ref <- na.omit(data[data[[cluster_col]] == cl, var])
      data_others <- na.omit(data[data[[cluster_col]] != cl, var])
      data_global <- na.omit(data[[var]])
      
      if (length(data_ref) > 1 & length(data_others) > 1) {
        
        test_res <- wilcox.test(data_ref, data_others, exact = FALSE)
        eff <- cliff.delta(data_ref, data_others)$estimate
        
        med_ref <- median(data_ref)
        med_others <- median(data_others)
        
        direction <- ifelse(med_ref > med_others, "higher", "lower")
        median_diff <- med_ref - med_others
        
        z_score <- (mean(data_ref) - mean(data_global)) / sd(data_global)
        
      } else {
        test_res <- list(p.value = NA)
        eff <- NA
        direction <- NA
        median_diff <- NA
        z_score <- NA
      }
      
      results_list[[length(results_list) + 1]] <- data.frame(
        variable = var,
        type = "continuous",
        cluster = cl,
        p_value = test_res$p.value,
        effect_size = eff,
        direction = direction,
        value_cluster = med_ref,
        value_others = med_others,
        z_score = z_score
      )
    }
  }
  
  ### =========================
  ### 2. QUALITATIVE VARIABLES 
  ### =========================
  
  for (var in qualitative_vars) {
    
    for (cl in clusters) {
      
      # Table de contingence
      tab <- table(data[[var]], data[[cluster_col]] == cl)
      
      if (all(dim(tab) >= 2)) {
        
        test_res <- fisher.test(tab)
        
        # Proportions
        prop_cluster <- prop.table(tab, 2)[, "TRUE"]
        prop_others <- prop.table(tab, 2)[, "FALSE"]
        
        # On prend la catégorie dominante
        max_cat <- names(which.max(prop_cluster))
        
        direction <- paste0("over:", max_cat)
        
        effect <- test_res$estimate  # odds ratio
        
      } else {
        test_res <- list(p.value = NA)
        effect <- NA
        direction <- NA
        prop_cluster <- NA
        prop_others <- NA
      }
      
      results_list[[length(results_list) + 1]] <- data.frame(
        variable = var,
        type = "categorical",
        cluster = cl,
        p_value = test_res$p.value,
        effect_size = effect,
        direction = direction,
        value_cluster = NA,
        value_others = NA,
        z_score = NA
      )
    }
  }
  
  ### =========================
  ### 3. COMBINATION
  ### =========================
  
  results <- bind_rows(results_list)
  
  # Correction multiple par variable
  results <- results %>%
    group_by(variable) %>%
    mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
    ungroup()
  
  return(results)
}
results_all <- compute_all_features(
  data = merged_meta,
  cluster_col = "cluster",
  continuous_vars = continuous_vars,
  qualitative_vars = qualitative_vars
)

write.table(results_all,"results_all.tsv")

############################################################
# 5. Differential Expression of genes per Cluster
############################################################
#### Requires the file from RASFLOW ####
expr_df <- read_tsv(file.path(path_output, "Human_transcripts.tsv"))
expr_df <- expr_df %>% tibble::column_to_rownames("Gene")
expr_df <- expr_df[, merged_meta$sample]
expr_filt <- expr_df[rowSums(expr_df > 0) > 0, ]

dge <- DGEList(counts=expr_filt)
dge <- calcNormFactors(dge)


run_DE <- function(cluster_col) {
  design <- model.matrix(reformulate(cluster_col), data=merged_meta)
  y <- estimateDisp(dge, design, robust=TRUE)
  fit <- glmQLFit(y, design, robust=TRUE)
  qlt <- glmQLFTest(fit)
  top <- topTags(qlt, n=nrow(expr_filt))
  return(top)
}

cluster_vars <- paste0("Clus_", 1:k)
DE_results <- lapply(cluster_vars, run_DE)
names(DE_results) <- cluster_vars

lapply(names(DE_results), function(cn)
  write.table(DE_results[[cn]], file.path(path_output, paste0("DE_", cn, ".tsv")),
              sep="\t", quote=FALSE))


############################################################
# 6. GO Enrichment Analysis with Dotplot
############################################################
cat("=== Step 5: GO Enrichment Analysis ===\n")

# --- Create ranked lists ---
ranked_lists <- lapply(DE_results, function(res)
  sort(setNames(res$table$logFC, rownames(res$table)), decreasing=TRUE))

# --- Run GSEA GO enrichment ---
enrich_res <- compareCluster(
  geneClusters = ranked_lists,
  fun = "gseGO",
  OrgDb = org.Hs.eg.db,
  keyType = "ENSEMBL",
  ont = "BP",
  eps = 0,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  minGSSize = 15
)

enrich_res <- enrichplot::pairwise_termsim(enrich_res)
enrich_res <- setReadable(enrich_res, 'org.Hs.eg.db', 'ENSEMBL')

write.table(enrich_res@compareClusterResult,
            file.path(path_output, "GSEA_results.tsv"),
            sep="\t", quote=FALSE)

# --- Dot plot visualization ---
dotplot(enrich_res, showCategory=5, font.size=12) +
ggsave(file.path(path_output, "GSEA_dotplot.pdf"), width=8, height=6)
