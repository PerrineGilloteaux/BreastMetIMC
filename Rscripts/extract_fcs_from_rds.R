library(flowCore)
library(readxl)
setwd("D:\\GIT\\BreastMetIMC\\Input_files")
# Charger le fichier backup_output.rds
output <- readRDS("backup_output.rds")

# Vérifier si write.flowSet existe dans flowCore
if (exists("write.flowSet", where = as.environment("package:flowCore"))) {
  
  write.flowSet(output$fcs)
} else {
  warning("write.flowSet n'est pas disponible dans flowCore. Utilisez write.FCS ou une autre méthode.")
}