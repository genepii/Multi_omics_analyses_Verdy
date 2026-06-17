# ==============================================================================
# 1. PARAMÈTRES, BIBLIOTHÈQUES & CHEMINS
# ==============================================================================

# Chargement de l'ensemble des bibliothèques nécessaires
library(stringr)
library(dplyr)
library(tibble)
library(tidyr)
library(phyloseq)
library(decontam)

# Configuration du répertoire de travail
setwd("C:\\Users\\destrasgr\\Documents\\Test")

# Fichiers d'entrée et de sortie
fileinput           <- "RNAbac_raw.txt" 
file_output          <- "RNAbac" 
fileinput_decontam   <- "decotnam_metabac.csv"      # Fichier de mapping pour decontam
correspondance       <- "correspondance_RNAbac.csv" # Fichier de mapping pour la soustraction
metadata_file        <- "metadata.csv"         # Uniquement les patients inclus

threshold            <- 1  # Seuil pour considérer un taxon comme positif

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

# Chargement des tables de comptages et métadonnées
meta <- read.delim(fileinput, stringsAsFactors = FALSE)
metadata_decontam <- read.csv2(fileinput_decontam, row.names = 1)
metadata <- read.csv2(metadata_file)

# Initialisation de la table d'analyse
cor_analiz <- meta

# Suppression des lignes dont toutes les abondances d'échantillons sont à 0
cor_analiz <- cor_analiz[rowSums(cor_analiz[, -1]) > 0, ]

# Repérage et agrégation des lignes associées au Contrôle Interne (IC)
has_tag <- grepl(IC, cor_analiz$SAMPLEID)

if (any(has_tag)) {
  ref_rows <- cor_analiz[has_tag, ]
  other_rows <- cor_analiz[!has_tag, ]
  ref_rows$SAMPLEID <- IC
  
  # Agrégation par somme sur l'identifiant IC unique
  aggregated_ref <- aggregate(. ~ SAMPLEID, data = ref_rows, FUN = sum)
  
  # Reconstruction de la table globale
  cor_analiz <- rbind(other_rows, aggregated_ref)
}

# Renommer la première colonne d'identification
colnames(cor_analiz)[1] <- "Genera"


# ==============================================================================
# 3. PRÉ-FILTRAGE BAS