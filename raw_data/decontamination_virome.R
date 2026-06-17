# ==============================================================================
# 1. PARAMÈTRES, BIBLIOTHÈQUES & CHEMINS
# ==============================================================================

# Chargement des bibliothèques nécessaires
library(stringr)
library(phyloseq)
library(decontam)

# Configuration du répertoire de travail
setwd("C:\\Users\\destrasgr\\Documents\\Test")

# Fichiers d'entrée et de sortie
fileinput           <- "Virome_raw.txt" 
file_output          <- "Virome" 
fileinput_decontam   <- "decontam_virome.csv"      # Fichier de mapping pour decontam
correspondance       <- "correspondance_virome.csv"# Fichier de mapping pour la soustraction
metadata_file        <- "metadata.csv"        # Uniquement les patients inclus

threshold            <- 0.01  # Seuil pour considérer un taxon comme positif

# Paramètres de Normalisation
IC                   <- "mesvirus"  # Nom du contrôle interne (0 si aucun)
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

# Initialisation de la table d'analyse et remplacement des NA par 0
cor_analiz <- meta
cor_analiz[is.na(cor_analiz)] <- 0

# Supprimer les lignes dont la somme des abondances est nulle (sans la 1ère colonne)
cor_analiz <- cor_analiz[rowSums(cor_analiz[, -1]) > 0, ]

# Agrégation des lignes dupliquées pour le contrôle interne (IC)
for (ref_tag in c(IC)) {
  if (ref_tag != 0 && any(duplicated(cor_analiz$SAMPLEID[grep(ref_tag, cor_analiz$SAMPLEID)]))) {
    ref_rows <- cor_analiz[grep(ref_tag, cor_analiz$SAMPLEID), ]
    other_rows <- cor_analiz[-grep(ref_tag, cor_analiz$SAMPLEID), ]
    aggregated_ref <- aggregate(. ~ SAMPLEID, data = ref_rows, FUN = sum)
    cor_analiz <- rbind(other_rows, aggregated_ref)
  }
}

# Renommer la première colonne d'identification
colnames(cor_analiz)[1] <- "Genera"


# ==============================================================================
# 3. NORMALISATION
# ==============================================================================

if (IC != 0) {
  # Normalisation par ratio (Contrôle Interne) : conserve les échantillons avec IC > 0
  ic_row_idx <- grep(IC, cor_analiz$Genera)
  valid_cols <- c(1, which(cor_analiz[ic_row_idx, -1] > 0) + 1)
  cor_analiz <- cor_analiz[, valid_cols]
  
  # Division des comptes de chaque échantillon par le compte de son IC associé
  ic_vector <- as.numeric(cor_analiz[grep(IC, cor_analiz$Genera), -1])
  cor_analiz[, -1] <- sweep(cor_analiz[, -1], 2, ic_vector, FUN = "/")
} else {
  # Normalisation standard en CPM (Counts per Million) si pas d'IC
  rpm_scale <- apply(cor_analiz[, -1], 2, sum) / 1e6
  cor_analiz[, -1] <- sweep(cor_analiz[, -1], 2, rpm_scale, FUN = "/")
}

# Nettoyage du bruit de fond sous le seuil critique
cor_analiz[cor_analiz <= threshold] <- 0


# ==============================================================================
# 4. DÉCONTAMINATION VIA LE PACKAGE DECONTAM (Méthode Prévalence)
# ==============================================================================

# Passage de la colonne Genera en noms de lignes avant création de l'objet Phyloseq
row.names(cor_analiz) <- cor_analiz$Genera
cor_analiz <- cor_analiz[, -1]

# Création de l'objet Phyloseq
OTU <- otu_table(cor_analiz, taxa_are_rows = TRUE)
sampledata <- sample_data(data.frame(metadata_decontam))
ps <- phyloseq(OTU, sampledata)

# Définition des contrôles négatifs
sample_data(ps)$is.neg <- sample_data(ps)$Sample_or_Control == "Control Sample"

# Détection des contaminants
contamdf.prev <- isContaminant(ps, method = "prevalence", neg = "is.neg", threshold = 0.5, normalize = FALSE)
table(contamdf.prev$contaminant)
contamdf.prev[which(contamdf.prev$contaminant), ]

# Extraction des contaminants détectés
contaminants <- row.names(contamdf.prev[which(contamdf.prev$contaminant), ])
dfconta <- data.frame(conta = contaminants)

print(row.names(contamdf.prev[which(contamdf.prev$contaminant), ]))
contaminants_final <- dfconta$conta

# Retrait des contaminants de la table principale
cor_analiz <- cor_analiz[-grep(paste0(contaminants_final, collapse = "|"), row.names(cor_analiz)), ]


# ==============================================================================
# 5. DÉCONTAMINATION PAR SOUSTRACTION DES CONTRÔLES
# ==============================================================================

coord_df <- read.csv2(correspondance, stringsAsFactors = FALSE)
normalized_table <- cor_analiz

if (deconta == "yes") {
  if (paired == "yes") {
    # Décontamination appariée : Soustraction du contrôle (T) de son échantillon (P)
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
      normalized_table[, -1] <- sweep(cor_analiz[, -1], 1, max_control, FUN = "-")
    }
  }
}

# Remplacement des valeurs négatives après soustraction par un plancher à 0
normalized_table[, -1] <- lapply(normalized_table[, -1], function(x) pmax(x, 0))

# Suppression des lignes vides et ré-application du filtre de seuil
final_table <- normalized_table[rowSums(normalized_table[, -1]) > 0, ]
final_table[final_table <= threshold] <- 0


# ==============================================================================
# 6. NETTOYAGE COSMÉTIQUE ET TRANSPOSITION
# ==============================================================================

# Nettoyage des préfixes "X." générés dans les noms de colonnes
colnames(final_table) <- gsub("X\\.", "", colnames(final_table))

# Retrait des lignes taxonomiques indésirables ("Toti" ou "Retro")
final_table <- final_table[!grepl("Toti|Retro", rownames(final_table)), ]

# Transposition de la table (les lignes deviennent colonnes)
finaltable_transposed <- as.data.frame(t(final_table))

# Suppression des colonnes vides (taxons nuls partout après transposition)
finaltable_transposed <- finaltable_transposed[, colSums(finaltable_transposed) > 0]


# ==============================================================================
# 7. FILTRAGE ET CONSERVATION DES PATIENTS INCLUS
# ==============================================================================

# Conservation exclusive des échantillons présents dans la table clinique
finaltable_filtered <- finaltable_transposed[rownames(finaltable_transposed) %in% metadata$sample, ]
finaltable_filtered <- finaltable_filtered[, colSums(finaltable_filtered) > 0]


# ==============================================================================
# 8. ÉCRITURE DES RÉSULTATS FINAUX
# ==============================================================================

write.csv2(finaltable_transposed, paste0(file_output, "\\normalized_decontaminated_table_", file_output, ".csv"), row.names = TRUE)
write.csv2(finaltable_filtered, paste0(file_output, "\\normalized_decontaminated_patients_", file_output, ".csv"), row.names = TRUE)