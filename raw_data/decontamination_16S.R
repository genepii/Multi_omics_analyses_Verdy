# ==============================================================================
# 1. PARAMÈTRES, BIBLIOTHÈQUES & CHEMINS
# ==============================================================================

# Chargement des bibliothèques nécessaires
library(stringr)
library(phyloseq)
library(decontam)
library(dplyr)

# Configuration du répertoire de travail
setwd("C:\\Users\\destrasgr\\Documents\\Test2")

# Fichiers d'entrée et de sortie
fileinput           <- "Bac16S_raw.txt" 
file_output          <- "Bac16S" 
fileinput_decontam   <- "decontam_16S.csv"    # Fichier de mapping pour decontam
correspondance       <- "correspondance_16S.csv"# Fichier de mapping pour la soustraction
metadata_file        <- "metadata.csv"    # Uniquement les patients inclus

# Seuils et contrôles du pipeline
threshold            <- 0.01  # Seuil pour considérer un taxon comme positif

# Paramètres de Normalisation
IC                   <- 0     # Nom du contrôle interne (0 si aucun)
ref                  <- "IC"

# Paramètres de Décontamination
deconta              <- "yes" # "yes" ou "no" pour activer la décontamination
paired               <- "no"  # "yes" ou "no" pour une décontamination appariée

# Création du dossier de sortie si nécessaire
dir.create(file_output, showWarnings = FALSE)


# ==============================================================================
# 2. CHARGEMENT DES DONNÉES & NETTOYAGE INITIAL
# ==============================================================================

# Chargement des tables brutes
meta <- read.delim(fileinput, stringsAsFactors = FALSE)
metadata_decontam <- read.csv2(fileinput_decontam, row.names = 1)
metadata <- read.csv2(metadata_file)

# Initialisation de la table d'analyse
cor_analiz <- meta

# Renommer la première colonne et rendre les noms de genres uniques
colnames(cor_analiz)[1] <- "Genera"
cor_analiz$Genera <- sapply(1:nrow(cor_analiz), function(x) paste0(cor_analiz$Genera[x], "_", x))

# Supprimer les lignes dont la somme des abondances est nulle
cor_analiz <- cor_analiz[rowSums(cor_analiz[, -1]) > 0, ]

# Assigner les genres en noms de lignes et supprimer la colonne textuelle
row.names(cor_analiz) <- cor_analiz$Genera
cor_analiz <- cor_analiz[, -grep("Genera", colnames(cor_analiz))]


# ==============================================================================
# 3. DÉCONTAMINATION VIA LE PACKAGE DECONTAM (Méthode Prévalence)
# ==============================================================================

# Création de l'objet Phyloseq
OTU <- otu_table(cor_analiz, taxa_are_rows = TRUE)
sampledata <- sample_data(data.frame(metadata_decontam))
ps <- phyloseq(OTU, sampledata)

# Définition des contrôles négatifs
sample_data(ps)$is.neg <- sample_data(ps)$Sample_or_Control == "Control Sample"

# Détection des contaminants
contamdf.prev <- isContaminant(ps, method = "prevalence", neg = "is.neg", threshold = 0.5, normalize = FALSE)
table(contamdf.prev$contaminant)

# Extraction et filtrage des contaminants
contaminants <- row.names(contamdf.prev[which(contamdf.prev$contaminant), ])
dfconta <- data.frame(conta = contaminants)
print(row.names(contamdf.prev[which(contamdf.prev$contaminant), ]))

contaminants_final <- dfconta$conta

# Retrait des contaminants de la table principale
if (length(contaminants_final) > 0) {
  cor_analiz <- cor_analiz[-grep(paste0(contaminants_final, collapse = "|"), row.names(cor_analiz)), ]
}


# ==============================================================================
# 4. NORMALISATION INITIALE (RPM / CPM)
# ==============================================================================

if (IC != 0) {
  # Normalisation par ratio (Contrôle Interne)
  ic_row_idx <- grep(IC, row.names(cor_analiz))
  valid_cols <- c(1, which(cor_analiz[ic_row_idx, ] > 0))
  cor_analiz <- cor_analiz[, valid_cols]
  
  ic_vector <- as.numeric(cor_analiz[grep(IC, row.names(cor_analiz)), ])
  cor_analiz <- sweep(cor_analiz, 2, ic_vector, FUN = "/")
} else {
  # Normalisation standard en CPM (Counts per Million)
  rpm_scale <- apply(cor_analiz, 2, sum) / 1e6
  cor_analiz <- sweep(cor_analiz, 2, rpm_scale, FUN = "/")
}


# ==============================================================================
# 5. DÉCONTAMINATION PAR SOUSTRACTION DES CONTRÔLES
# ==============================================================================

coord_df <- read.csv2(correspondance, stringsAsFactors = FALSE)
normalized_table <- cor_analiz

if (deconta == "yes") {
  if (paired == "yes") {
    # Décontamination appariée : Soustraction de chaque contrôle (T) de son échantillon (P)
    paired_list <- list()
    for (i in 1:nrow(coord_df)) {
      p_col <- coord_df$P[i]
      t_col <- coord_df$T[i]
      if (p_col %in% colnames(cor_analiz) && t_col %in% colnames(cor_analiz)) {
        paired_list[[p_col]] <- cor_analiz[[p_col]] - cor_analiz[[t_col]]
      }
    }
    normalized_table <- data.frame(Genera = cor_analiz$Genera, as.data.frame(paired_list), check.names = FALSE)
  } else {
    # Décontamination non-appariée : Soustraction de la valeur MAX des contrôles (T)
    t_cols <- coord_df$T[!is.na(coord_df$T)]
    t_cols <- t_cols[t_cols %in% colnames(cor_analiz)]
    
    if (length(t_cols) > 0) {
      max_control <- if (length(t_cols) > 1) apply(cor_analiz[, t_cols], 1, max) else cor_analiz[, t_cols]
      normalized_table <- sweep(cor_analiz, 1, max_control, FUN = "-")
    }
  }
}

# Remplacement des valeurs négatives générées par la soustraction par 0
normalized_table_filtered <- as.data.frame(lapply(normalized_table, function(x) pmax(x, 0)))
row.names(normalized_table_filtered) <- row.names(normalized_table)

# Suppression des lignes devenues vides (uniquement des 0)
final_table <- normalized_table_filtered[rowSums(normalized_table_filtered) > 0, ]

# Filtrage des taxons annotés "NA" sur des niveaux spécifiques
final_table <- final_table[-grep(";  NA", sapply(row.names(final_table), function(x) paste(str_split(x, ";")[[1]][4], ";", str_split(x, ";")[[1]][5]))), ]

# Nettoyage cosmétique des noms de colonnes et retrait des lignes "Equine"
colnames(final_table) <- gsub("X\\.", "", colnames(final_table))
final_table <- final_table[!grepl("Equine", rownames(final_table)), ]


# ==============================================================================
# 6. CALCUL DES ABONDANCES RELATIVES, SEUIL & CHARGE BACTÉRIENNE
# ==============================================================================

# Transformation en abondances relatives
df_rel <- final_table %>%
  mutate(across(everything(), ~ .x / sum(.x)))

# Application du seuil critique
df_rel[df_rel < threshold] <- 0

# Conservation uniquement des échantillons présents dans les métadonnées cliniques
df_rel <- df_rel[, colnames(df_rel) %in% metadata$sample]

# Nettoyage et agrégation par taxonomie cible
df_rel$Genera <- sapply(row.names(df_rel), function(x) paste(str_split(x, "; ")[[1]][5], ";", str_split(x, ";")[[1]][6]))
df_rel2 <- aggregate(. ~ Genera, data = df_rel, FUN = sum)
row.names(df_rel2) <- df_rel2$Genera
df_rel2 <- df_rel2[, -1]

finaltable_filtered <- df_rel2

# Multiplication par la charge bactérienne (bacterial_load)
names <- colnames(finaltable_filtered)
copies_vector <- as.numeric(sapply(names, function(x) metadata$bacterial_load[metadata$sample == x]))
final_table_copies <- sweep(finaltable_filtered, 2, copies_vector, FUN = "*")


# ==============================================================================
# 7. TRANSPOSITION & FILTRAGE FINAL
# ==============================================================================

# Transposition de la table des copies
finaltable_transposed <- as.data.frame(t(final_table_