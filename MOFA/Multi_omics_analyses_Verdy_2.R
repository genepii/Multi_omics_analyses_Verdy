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
library(vegan)

############################################################
# 1. Computation of diversity indexes
############################################################

setwd("C:\\Users\\destrasgr\\Documents\\Test\\New_Test")
# --- Read omics tables counts---
bac_counts     <- read.table("Bac16S\\16S_data.txt", sep="\t", header=TRUE, row.names=1, dec=",")
RNAbac_counts  <- read.table("RNAbac\\RNAbac_data.txt", sep="\t", header=TRUE, row.names=1, dec=",")
vir_counts <- read.table("Virome\\virome_OTU_data.txt", sep="\t", header=TRUE, row.names=1, dec=",",check.names = FALSE)
metadata       <- read.table("metadata.csv", sep=";", header=TRUE, row.names=1)


# Compute diversity indexes
bac_numeric <- as.data.frame(lapply(bac_counts, as.numeric))
bac_numeric <- bac_numeric[,1:(ncol(bac_numeric)-3)]

meta_bac <- data.frame(
  index_shannon_bac   = diversity(bac_numeric, index = "shannon"),
  copies_total_bac = rowSums(bac_numeric),
  index_richness_bac  = rowSums(bac_numeric > 0),
  stringsAsFactors = FALSE
)
rownames(meta_bac) <- rownames(bac_counts)


RNAbac_numeric <- as.data.frame(lapply(RNAbac_counts, as.numeric))
RNAbac_numeric <- RNAbac_numeric[,1:(ncol(RNAbac_numeric)-3)]

meta_RNAbac <- data.frame(
  index_shannon_RNAbac   = diversity(RNAbac_numeric, index = "shannon"),
  RPM_total_RNAbac = rowSums(RNAbac_numeric),
  index_richness_RNAbac  = rowSums(RNAbac_numeric > 0),
  stringsAsFactors = FALSE
)
rownames(meta_RNAbac) <- rownames(RNAbac_counts)

vir_numeric <- as.data.frame(lapply(vir_counts, as.numeric))
vir_numeric <- vir_numeric[,1:(ncol(vir_numeric)-5)]

meta_vir <- data.frame(
  index_shannon_vir   = diversity(vir_numeric, index = "shannon"),
  VDR_total_vir = rowSums(vir_numeric),
  index_richness_vir  = rowSums(vir_numeric > 0),
  stringsAsFactors = FALSE
)
rownames(meta_vir) <- rownames(vir_counts)

library(stringr)
df_vir<-as.data.frame(t(vir_counts))
df_vir$Family<-unlist(sapply(row.names(df_vir),function(x)str_split(x,";")[[1]][5]))
df_vir<-aggregate(.~Family,df_vir,sum)
row.names(df_vir)<-df_vir$Family
df_vir<-df_vir[,-1]
df_vir<-t(df_vir)

############################################
##Compute virulent et Eukaryotes
##########################################
# 1. Chargement des fichiers (pensez à adapter les chemins d'accès)
# On utilise check.names = FALSE au cas où vos noms de taxons/vOTU contiennent des ";"
phabox    <- read.table("Virome\\phabox.txt", header = TRUE, sep = "\t", check.names = FALSE)
vOTU      <- read.table("Virome\\vOTU.txt", header = TRUE, sep = "\t", check.names = FALSE)
vir_counts <- read.table("Virome\\virome_OTU_data.txt", header = TRUE, sep = "\t", dec=",",check.names = FALSE)

# ==============================================================================
# ÉTAPE 2 : Création de la table de correspondance (Mapping Table)
# ==============================================================================

#colnames(vir_counts)<-c("sample",unlist(sapply(colnames(vir_counts),function(x)str_split(x,";")[[1]][5])))

# On part de la table de mapping 'vOTU' pour lier le nom du taxon (taxonomy) 
# au nom d'accession (OTU_name), puis on va chercher le "TYPE" dans phabox.
mapping_types <- vOTU %>%
  select(Taxonomy, OTU_name) %>%
  left_join(phabox %>% select(Accession, TYPE), by = c("OTU_name" = "Accession")) %>%
  
  # Application stricte de vos règles de priorité :
  mutate(Final_TYPE = case_when(
    # 1. Règle Eucaryotes (Basée sur les familles spécifiques)
    str_detect(Taxonomy, "Anelloviridae|Papillomaviridae|Potyviridae") ~ "Eukaryote",
    
    # 2. Règle Unclassified (Cas exact ou partiel selon votre taxonomie)
    str_detect(Taxonomy, "Viruses_Unclassified_Unclassified_Unclassified_Unclassified") ~ "Unclassified",
    
    # 3. Règles PhaBOX (Uniquement si ce n'est ni Eucaryote ni Unclassified)
    TYPE == "virulent"     ~ "virulent",
    TYPE == "temperate"    ~ "temperate",
    TYPE == "-" ~ "Undetermined",
    
    # 4. Tout le reste devient "Bacteriophage_Autre" (pour les phages non classés par PhaBOX)
    TRUE                   ~ "Bacteriophage_Autre"
  )) %>%
  select(Taxonomy, Final_TYPE)
# ==============================================================================
# ÉTAPE 3 : Pivotement et calcul des abondances
# ==============================================================================

# Si vos patients sont actuellement en LIGNES dans vir_counts (avec une colonne Patient),
# nous allons transformer la table pour faire les calculs par patient.

# Option A : Vos patients sont en LIGNES (Rows), les taxons en COLONNES
# (On suppose que la première colonne s'appelle "Patient")
resultat_final <- vir_counts %>%
  # On suppose que la 1ère colonne s'appelle "Patient" et les autres sont les taxons
  pivot_longer(cols = -sample, names_to = "Taxonomy", values_to = "Abondance") %>%
  left_join(mapping_types, by = "Taxonomy") %>%
  
  # Agrégation par Patient et par la catégorie nettoyée
  group_by(sample, Final_TYPE) %>%
  summarise(Total_Abondance = sum(Abondance, na.rm = FALSE), .groups = "drop") %>%
  
  # Pivot pour mettre chaque catégorie en colonne
  pivot_wider(names_from = Final_TYPE, values_from = Total_Abondance, values_fill = 0)  %>%
  
  # Calculs des 5 colonnes finales demandées
  mutate(
    # 1. Somme de TOUS les bactériophages (les 3 types PhaBOX + les "autres" par défaut)
    bacteriophages = virulent + temperate + Undetermined,
    
    # 2. Somme des Eucaryotes
    eukaryotes     = Eukaryote,
    
    # 3, 4 et 5. Abondances relatives PARMI les bactériophages du patient (exclut Eucaryotes et Unclassified du dénominateur)
    virulent_rel     = if_else(bacteriophages > 0, virulent / bacteriophages, 0),
    temperate_rel    = if_else(bacteriophages > 0, temperate / bacteriophages, 0),
    Undetermined_rel = if_else(bacteriophages > 0, Undetermined / bacteriophages, 0)
  ) 

###################################
##input MOFA
###########
# --- Read clinical metadata (separate temporary clinical metadata files initially created for  mono-omic analyses) ---
bac_counts <- read.table("Bac16S\\16S_data.txt", sep="\t", header=TRUE)
RNAbac_counts <- read.table("RNAbac\\RNAbac_data.txt", sep="\t", header=TRUE)
vir_counts <- read.table("Virome\\Virome_data.txt", sep="\t", header=TRUE)

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

row.names(merged_counts)<-merged_counts$sample
merged_counts<-merged_counts[,-1]

# --- Filter features present in ≥25% samples ---

# 1. Define the list of columns to exclude from the filtering condition
excluded_cols <- c("VDR_total_vir", "index_shannon_vir", "index_richnes_vir","RPM_total_bacrep", "index_shannon_bacrep", "index_richnes_bacrep","copies_total_bac", "index_shannon_bac", "index_richnes_bac","virulent", "temperate", "Undeterminate")

# 2. Select the columns that should actually be tested
candidate_cols <- setdiff(colnames(merged_counts), excluded_cols)

# 3. Apply your condition only to those candidate columns
passed_filter <- apply(merged_counts[, candidate_cols, drop = FALSE], 2, function(x) {
  sum(as.numeric(x) > 0) >= 0.25 * nrow(merged_counts)
})

# 4. Combine the columns that passed the test with the excluded ones
cols_to_keep <- c(names(passed_filter)[passed_filter], excluded_cols)

# 5. Subset the original data frame

merged_counts_filt <- merged_counts[, colnames(merged_counts) %in% cols_to_keep, drop = FALSE]



# --- Standardize features with z scores ---
mofa_scaled <- as.data.frame(lapply(merged_counts_filt, scale, center=TRUE, scale=TRUE))
mofa_scaled$sample<-row.names(merged_counts_filt)


library(reshape2)
mofa_long <- melt(mofa_scaled, id="sample")

colnames(mofa_long) <- c("sample", "feature", "value")
mofa_long$view<-sapply(mofa_long$feature,function(x)ifelse(is.element(x,colnames(vir_counts)),"vir","bac"))
mofa_long<-mofa_long[-grep("Bacteriophages",mofa_long$feature),]

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
