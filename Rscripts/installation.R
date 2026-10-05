# Installation

# Liste des packages à installer
install.packages("BiocManager")

BiocManager::install("flowCore")
BiocManager::install("limma")
BiocManager::install("ComplexHeatmap")

install.packages("remotes")
#remotes::install_github("garretrc/ggvoronoi")
remotes::install_github("garretrc/ggvoronoi", dependencies = TRUE, build_opts = c("--no-resave-data"))

packages <- c(
"readxl",
"stringr",
"matrixcalc",
"Hmisc",
"reshape2",
"dplyr",
"plotrix",
"multcomp",
"flowCore",
"sf",
"clusterSim",
"limma",
"corrplot",
"packcircles",
"cowplot",
"autoimage",
"ggplot2",
"ggpubr",
"ggiraphExtra",
"ggvoronoi",
"ggridges",
"gridExtra",
"igraph",
"qgraph",
"circlize",
"scales",
"RColorBrewer",
"ComplexHeatmap",
"pheatmap",
"pals",
"plot3D",
"akima",
"basetheme")



# Vérifier lesquels ne sont pas encore installés
packages_to_install <- packages[!(packages %in% installed.packages()[,"Package"])]

# Installer uniquement ceux manquants
if(length(packages_to_install) > 0) {
  install.packages(packages_to_install, dependencies = TRUE)
}

# Charger tous les packages
lapply(packages, require, character.only = TRUE)
a