# Installer le package digest si nécessaire
if (!requireNamespace("digest", quietly = TRUE)) {
  install.packages("digest")
}

library(digest)

# --- Fonction pour calculer le MD5 d'une chaîne ---
md5_string <- function(text) {
  if (!is.character(text) || length(text) != 1) {
    stop("Entrée invalide : veuillez fournir une seule chaîne de caractères.")
  }
  digest(text, algo = "md5", serialize = FALSE)
}

# --- Fonction pour calculer le MD5 d'un fichier ---
md5_file <- function(filepath) {
  if (!file.exists(filepath)) {
    stop("Fichier introuvable : ", filepath)
  }
  digest(file = filepath, algo = "md5")
}

tmp_file <- "..\\Input_Files\\backup_output.rds"

hash_fichier <- md5_file(tmp_file)
cat("MD5 du fichier :", hash_fichier, "\n")