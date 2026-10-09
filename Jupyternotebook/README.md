# Analyse spatiale de données PhenoCycler avec SpatialData

Ce notebook Jupyter analyse des données de **PhenoCycler (CODEX)**. Il part des images, de la segmentation et des intensités par cellule, et les rassemble dans le format standard **[SpatialData](https://spatialdata.scverse.org/en/stable/)**. Il refait ensuite les analyses du pipeline R *BreastMetIMC* (scripts `IMCpipeline_1` à `5`) :

1. **Identification des types cellulaires** : clustering FlowSOM, puis annotation manuelle assistée.
2. **Abondances** : proportions de chaque type cellulaire par échantillon, comparaisons LUM vs TNC et primaire vs métastase.
3. **Voisinage** : marqueurs par type, corrélations, voisins les plus proches des cellules tumorales (CK+).
4. **Distances** : distance de chaque cellule au type cellulaire le plus proche, proximité des macrophages, cartes de proximité.
5. **Réseaux** : représentation en graphe des distances entre types cellulaires.

Aucune connaissance en programmation n'est nécessaire pour l'utiliser, mais il faut **modifier quelques paramètres** (chemins, réglages) et **prendre quelques décisions biologiques**, signalées dans le notebook par le symbole 🖐️.

---

## Sommaire

1. [Ce dont vous avez besoin](#1-ce-dont-vous-avez-besoin)
2. [Installation (une seule fois)](#2-installation-une-seule-fois)
3. [Préparer vos données](#3-préparer-vos-données)
4. [Lancer le notebook](#4-lancer-le-notebook)
5. [Déroulé de l'analyse, étape par étape](#5-déroulé-de-lanalyse-étape-par-étape)
6. [Fichiers produits](#6-fichiers-produits)
7. [Relancer une analyse ou ajouter des régions](#7-relancer-une-analyse-ou-ajouter-des-régions)
8. [Problèmes fréquents](#8-problèmes-fréquents)
9. [Limites et précautions d'interprétation](#9-limites-et-précautions-dinterprétation)

---

## 1. Ce dont vous avez besoin

**Ordinateur**
- Windows 10/11, macOS ou Linux.
- **16 Go de mémoire vive minimum** (32 Go conseillés si vous avez beaucoup de régions ou de cellules).
- De l'espace disque libre : environ la taille de vos TIF, pour la copie au format Zarr.

**Fichiers fournis**
- `PhenoCycler_SpatialData_pipeline_local.ipynb` : le notebook.
- `requirements.txt` : la liste des logiciels à installer.

**Vos données PhenoCycler**, dans l'organisation standard du PhenoCycler (détaillée en section 3)
- Pour chaque ROI : les images TIF (une par canal), le FCS compensé, la carte de segmentation et `channelNames.txt`.
- Pour analyser plusieurs ROI : un tableau Excel qui les liste (section 3.3).

---

## 2. Installation (une seule fois)

L'installation crée un « environnement » Python isolé, qui ne perturbe pas les autres logiciels de l'ordinateur.

### 2.1 Installer Miniconda

1. Téléchargez **Miniconda** pour votre système : <https://docs.conda.io/en/latest/miniconda.html>
2. Installez-le en gardant les options par défaut.
3. Ouvrez un terminal :
   - **Windows** : menu Démarrer → **Anaconda Prompt (miniconda3)** ;
   - **macOS** : application **Terminal** ;
   - **Linux** : votre terminal habituel.

> Si votre institut fournit déjà Anaconda ou un environnement Python, demandez à votre service informatique ou à la plateforme de bio-informatique : les commandes ci-dessous restent valables.

### 2.2 Créer l'environnement et installer les paquets

Dans le terminal, placez-vous dans le dossier qui contient `requirements.txt` et le notebook, par exemple :

```bash
cd D:\Analyses\PhenoCycler          (Windows)
cd ~/Analyses/PhenoCycler           (macOS / Linux)
```

Puis tapez ces trois commandes, une par une :

```bash
conda create -n spatial python=3.12 -y
conda activate spatial
pip install -r requirements.txt
```

L'installation prend 5 à 15 minutes. Elle n'est à faire **qu'une fois**.

> **Pourquoi Python 3.12 ?** Certains paquets scientifiques ne sont pas encore disponibles pour les toutes dernières versions de Python. Les versions 3.10, 3.11 et 3.12 conviennent.

---

## 3. Préparer vos données

### 3.1 Organisation des dossiers

Le notebook lit les dossiers **tels que les produit le PhenoCycler**. Une acquisition (ROI) = un dossier. Plusieurs expériences peuvent être analysées ensemble :

```
D:/ELOISE/
├── regions_metadata.xlsx                    ← liste des ROI à traiter + description (voir 3.3)
├── EXP100/
│   ├── channelNames.txt                     (ici, ou dans chaque dossier ROI)
│   └── ROI_xxx/
└── EXP112/
    └── ROI_9x7/                             ← un dossier ROI
        ├── channelNames.txt                 (facultatif si présent dans EXP112/)
        ├── processed/
        │   ├── stitched/reg001/             reg001_cyc001_ch001_DAPI.tif, … (un TIF par canal)
        │   └── segm/segm-1/fcs/compensated/ …_reg001_compensated.fcs
        └── out/
            └── labels_cy002_ch1.tif         (carte de segmentation)
```

Ces emplacements à l'intérieur d'une ROI sont réglés une fois pour toutes dans la cellule 0.2 :

| Paramètre | Valeur par défaut |
|---|---|
| `IMAGES_SUBDIR` | `processed/stitched/reg001` |
| `FCS_SUBDIR` | `processed/segm/segm-1/fcs/compensated` |
| `LABELS_FILE_REL` | `out/labels_cy002_ch1.tif` |
| `CHANNELS_FILE_REL` | `channelNames.txt`, cherché dans la ROI puis dans les dossiers parents (EXP…) |

Si votre organisation est différente, modifiez ces quatre valeurs. Pour une ROI isolée qui ne suit pas l'organisation commune, utilisez plutôt les colonnes d'exception du tableau Excel (section 3.3).

**Ce qui est obligatoire pour chaque ROI :** le FCS et `channelNames.txt`. Une ROI sans l'un des deux est exclue (signalé à l'étape 0.3).

**Ce qui est facultatif :**
- **Les images :** les analyses utilisent les intensités du FCS ; les images servent à la visualisation et à la mesure sur masque (étape 0.7). Il n'est pas nécessaire d'avoir tous les canaux, sauf pour corriger le débordement d'un marqueur, qui demande son image.
- **La carte de labels :** sans elle, chaque cellule est représentée par un cercle et les mesures de forme sont approximées.

**Un point important :** toutes les acquisitions PhenoCycler portent le même nom interne (`reg001`). Le notebook nomme donc chaque ROI d'après son dossier : `<expérience>_<ROI>`, par exemple `EXP112_ROI_9x7`. C'est cet identifiant qui apparaît dans les graphiques et les fichiers.

### 3.2 La carte de segmentation (labels)

Deux types de fichiers sont acceptés et reconnus automatiquement :

- **Carte d'instances** (recommandé) : image 16 ou 32 bits où chaque cellule a son propre numéro.
- **Carte de contours colorés** : image 8 bits où les contours des cellules sont tracés avec quelques couleurs. C'est le cas de `labels_cy002_ch1.tif` (11 couleurs pour 26 666 cellules). Le notebook **reconstruit** alors les cellules : chaque zone fermée par un contour devient une cellule, puis est associée à la cellule du FCS dont le centre tombe dedans. Sur les données de test, 26 662 cellules sur 26 666 ont été retrouvées.

Si votre logiciel de segmentation peut exporter une vraie carte d'instances, préférez-la : c'est plus fiable.

### 3.3 Le tableau `regions_metadata.xlsx` (liste et description des ROI)

Ce tableau Excel a deux rôles : il dit **quelles ROI traiter**, et il **décrit chaque échantillon** pour les comparaisons (LUM vs TNC, primaire vs métastase, organes, patients appariés). Une ligne par ROI.

| Colonne | Obligatoire | Contenu | Exemple |
|---|---|---|---|
| `roi_path` | **oui** | dossier de la ROI | `D:\ELOISE\EXP112\ROI_9x7` |
| `region` | non | identifiant court (par défaut `<expérience>_<ROI>`) | `EXP112_ROI_9x7` |
| `include` | non | 1 = traiter, 0 = ignorer cette ligne | `1` |
| `sample_id` | non | nom de l'échantillon | `P_12` |
| `Case` | non | identifiant du patient | `SPC04` |
| `Anno` | non | site du prélèvement | `PBC`, `LUNG`, `BRAIN`, `LIVER`, `LN`… |
| `Tumor` | non | type histologique | `IDC`, `ILC` |
| `Phenotype` | non | sous-type tumoral | `LUM`, `TNC`, `HER2`, `CTRL` |
| `Met` | non | 0 = primaire, 1 = métastase | `0` |
| `Paired` | non | 1 si le patient a des prélèvements appariés | `1` |
| `pixel_size_um` | non | taille de pixel propre à cette ROI | `0.325` |
| `images_dir`, `fcs_file`, `labels_file`, `channels_file` | non | chemins d'exception pour une ROI qui ne suit pas l'organisation commune (chemin complet, ou relatif au dossier de la ROI) | `out/autre_labels.tif` |

**Les deux modes de la cellule 0.2 :**
- **Plusieurs ROI** : `REGIONS_METADATA` donne le chemin de ce tableau (par ex. `D:/ELOISE/regions_metadata.xlsx`). Toutes les lignes avec `include = 1` sont traitées.
- **Une seule ROI** : si ce fichier n'existe pas, seul le dossier `RAW_DIR` est traité, sans comparaison de groupes.

Sans les colonnes descriptives, toutes les analyses tournent, mais les comparaisons entre groupes sont ignorées.

> **Astuce :** lancez une première fois le notebook. Il écrit un modèle pré-rempli `regions_metadata_template.xlsx` dans le dossier de résultats (colonnes `region`, `roi_path`, `include` et description). Complétez-le dans Excel, enregistrez-le sous le nom indiqué dans `REGIONS_METADATA`, puis relancez depuis la cellule 0.2.

Les codes `PBC`, `PBC_tx`, `PBC_recur`… sont regroupés sous `PBC` (tumeur primaire) pour les graphiques par site, comme dans le pipeline R.

---

## 4. Lancer le notebook

À chaque nouvelle session de travail :

1. Ouvrez un terminal (Anaconda Prompt sous Windows).
2. Activez l'environnement et lancez JupyterLab :
   ```bash
   conda activate spatial
   jupyter lab
   ```
3. Votre navigateur s'ouvre. Dans le panneau de gauche, double-cliquez sur `PhenoCycler_SpatialData_pipeline_local.ipynb`.
4. **Exécutez les cellules une par une**, de haut en bout, avec **Maj + Entrée**. Le temps de lire les résultats entre les cellules est important, surtout aux étapes marquées 🖐️ et 🔍.

> Évitez « Run All » lors d'une première utilisation : plusieurs étapes demandent votre intervention, et l'annotation des clusters (étape 1.4) doit être faite à la main.

**Symboles utilisés dans le notebook**

| Symbole | Signification |
|---|---|
| 🖐️ **INTERACTIF** | Vous devez modifier un paramètre ou prendre une décision. |
| 🔍 **À INSPECTER** | Regardez le résultat avant de continuer. |
| ⚠️ **DIFFÉRENCE** | Le calcul diffère du pipeline R d'origine (expliqué sur place). |

Pour modifier un paramètre : cliquez dans la cellule, changez la valeur après le signe `=`, puis exécutez la cellule (Maj + Entrée).

---

## 5. Déroulé de l'analyse, étape par étape

### Partie 0 — Préparation des données

| Étape | Que faire |
|---|---|
| **0.1 Installation** | Exécutez la cellule. Si elle installe quelque chose, redémarrez le noyau (menu **Kernel → Restart Kernel**) puis passez à la cellule suivante. |
| **0.2 Paramètres** 🖐️ | Indiquez le tableau des ROI (`REGIONS_METADATA`) ou, pour une seule ROI, son dossier (`RAW_DIR`). Sous Windows, écrivez les chemins avec des `/` (`"D:/ELOISE/regions_metadata.xlsx"`). Vérifiez les sous-chemins des ROI (`IMAGES_SUBDIR`, `FCS_SUBDIR`, `LABELS_FILE_REL`) et la **taille de pixel** `PIXEL_SIZE_UM` (0,325 µm par défaut, à vérifier dans les métadonnées de votre acquisition) : toutes les distances en µm en dépendent. Les Zarr et les résultats sont écrits à côté du tableau Excel (ou dans la ROI unique), sauf si vous renseignez `PROJECT_DIR`. |
| **0.3 Fichiers et canaux** 🔍 | Un tableau indique, pour chaque ROI, le FCS trouvé, le nombre de TIF, la présence des labels et de `channelNames.txt`, et la taille de pixel. Les ROI sans FCS ou sans `channelNames.txt` sont exclues. Si les panels diffèrent entre expériences, seuls les marqueurs communs à toutes les ROI sont analysés ensemble. Le tableau complet des chemins est enregistré dans `regions_table_resolved.xlsx`. |
| **0.4 Conversion Zarr** | Chaque région est convertie **une seule fois** au format Zarr (dossier `spatialdata/`). Comptez quelques secondes à quelques minutes par région selon le nombre de canaux. Le tableau affiché indique le nombre de cellules retrouvées dans la carte de labels. |
| **0.5 Contrôle des labels** 🔍 | Un zoom montre le DAPI à côté des contours des cellules reconstruites. **Vérifiez que les contours entourent bien les noyaux.** Vous pouvez déplacer la fenêtre avec `QC_X0`, `QC_Y0` et `QC_SIZE` (en pixels). |
| **0.7 Mesures sur masque** 🖐️🔍 (facultatif) | Recalcule l'intensité moyenne de chaque cellule à partir des images et du masque, et produit une seconde mesure **corrigée du débordement latéral** (le signal d'une cellule qui « bave » sur ses voisines). Choisissez ensuite la mesure utilisée pour toute l'analyse avec `MEASUREMENT` : `"fcs"` (par défaut), `"mask_mean"` ou `"mask_spillover"`. Voir l'encadré ci-dessous. |

**Mesures sur masque et correction du débordement (étape 0.7)**

- **Pourquoi ?** Pour vérifier les intensités du FCS, et pour réduire les faux « doubles positifs » : par exemple un lymphocyte collé à une cellule tumorale qui paraît faiblement K8+.
- **Prérequis :** la région doit avoir une carte de labels et ses images. **Pour corriger un marqueur, son image TIF doit être présente.** Pour les canaux sans image, la valeur du FCS est conservée.
- **Contrôles à regarder :**
  1. La corrélation FCS ↔ masque doit être proche de 1. Sinon, la carte de labels ne correspond pas aux images.
  2. Le tableau γ (fraction de débordement par marqueur) et α (enrichissement du bord : positif pour un marqueur membranaire, négatif pour un marqueur nucléaire).
  3. Pour des paires de marqueurs de lignées exclusives (CD3/CD20, CD3/K8…), la corrélation et le pourcentage de cellules doubles positives doivent **baisser** avec `mask_spillover`.
- **Limites :** la correction est une approximation inspirée de la méthode REDSEA (Bai *et al.*, 2021). Sur un test simulé, elle réduit d'environ 70 % le signal parasite des cellules les plus contaminées, sans le supprimer totalement. Elle ne corrige ni le bruit de fond ni le recouvrement spectral entre canaux.
- **Après un changement de mesure**, relancez toute la suite et **refaites l'annotation** des clusters.

### Partie 1 — Types cellulaires

| Étape | Que faire |
|---|---|
| **0.6 Liste des types cellulaires** 🖐️ | `CLUSTERLEVELS` contient les noms possibles des types cellulaires. Elle est adaptée au panel PhenoCycler de test (CD31 → `Endo`, pas de marqueur pour les granulocytes ni les cellules dendritiques). **Adaptez-la à votre panel.** Les familles (CK, stroma, lymphocytes, myéloïdes) sont déduites des noms : un nom contenant `CK` est tumoral, `Mac` un macrophage, `Str` ou `Endo` du stroma. |
| **1.1 Métadonnées** | Lecture des colonnes descriptives du tableau des ROI, et écriture du modèle `regions_metadata_template.xlsx`. |
| **1.2 Transformation** 🖐️🔍 | Les intensités sont transformées par `arcsinh(intensité / COFACTOR)`. **Regardez les histogrammes** : pour un bon marqueur de lignée (CD3, CD20, K8, CD68…), on doit voir deux populations, négative et positive, bien séparées. Si tout est écrasé à gauche, diminuez `COFACTOR` ; si tout est étalé sans séparation, augmentez-le. La valeur par défaut est 150 ; testez par exemple 50, 150 et 300. Les canaux DAPI, Blank et Empty sont exclus automatiquement (`EXCLUDE_PATTERN`). |
| **1.3 Clustering** 🖐️ | FlowSOM regroupe les cellules en `N_CLUSTERS` groupes (30 par défaut, 50 dans le pipeline R). Mieux vaut un peu trop de clusters que pas assez : plusieurs clusters pourront recevoir le même nom. |
| **1.4 Annotation** 🖐️ | **L'étape la plus importante**, détaillée juste après ce tableau. |
| **1.5 Heatmap et UMAP** 🔍 | Vérifiez que chaque type annoté a un profil cohérent (par ex. `Tc` : CD3+ CD8+). L'UMAP interactive permet de colorer par type, région ou marqueur (`UMAP_COLOR_BY`) et de zoomer. |
| **1.6 Cartes spatiales** 🖐️ | Les types cellulaires sont affichés sur l'image. Choisissez la région (`MAP_REGION`), le canal de fond (`MAP_CHANNEL`) et la fenêtre (`WIN_X0`, `WIN_Y0`, `WIN_SIZE` en pixels ; `WIN_SIZE = 0` affiche la région entière). C'est un bon contrôle visuel de l'annotation : les cellules CK+ doivent se trouver dans les îlots tumoraux, etc. |

**Comment annoter les clusters (étape 1.4)**

1. Regardez la **heatmap** (marqueurs en colonnes, clusters en lignes, rouge = fort) et le **tableau des marqueurs dominants**.
2. Le notebook propose une **suggestion automatique** fondée sur des règles simples (CD3 + CD8 → `Tc`, K8 + PanCK → `CK8`…). **Ces suggestions sont grossières et doivent être vérifiées une par une.**
3. Dans l'**éditeur à menus déroulants**, choisissez pour chaque cluster le type cellulaire. Utilisez `NA` pour les clusters ininterprétables (débris, double marquages, bruit de fond).
4. Cliquez sur **Enregistrer**, puis exécutez la cellule **« Application de l'annotation »**.
5. Votre annotation est sauvegardée dans `merge_phenocycler.xlsx` et rechargée automatiquement aux exécutions suivantes. Vous pouvez aussi la modifier directement dans Excel.

> Si vous relancez le clustering avec d'autres paramètres, les numéros de clusters changent : il faut refaire l'annotation.

### Partie 2 — Abondances

S'exécute sans intervention. Elle produit :
- les tableaux de comptes et de proportions par échantillon ;
- les barres empilées (tous types, macrophages, stroma, lymphocytes), triées par proportion de macrophages ;
- **si les métadonnées le permettent**, les comparaisons LUM vs TNC (test de Wilcoxon ; seules les différences significatives sont affichées) et les graphiques « lollipop » ;
- la morphologie des cellules (aire en µm², grand axe en µm) par type.

### Partie 3 — Marqueurs et voisinage

| Étape | Que faire |
|---|---|
| **3.1 Heatmaps ciblées et corrélations** 🖐️ | Les listes `MARKERLIST_MAC`, `MARKERLIST_CK` et `MARKERLIST_STR` définissent les marqueurs étudiés pour les macrophages, les cellules tumorales et le stroma. Adaptez-les à votre panel. |
| **3.2 Voisins des cellules CK+** 🖐️ | Pour chaque cellule tumorale, le notebook identifie ses 3 voisines les plus proches et compte combien sont des macrophages (ou un autre type : `IMMUNECELLS`). `MIN_CELLS_CONDITION` (200 par défaut) écarte les conditions avec trop peu de cellules. Si un message « matrice trop petite » apparaît, diminuez ce seuil. |

### Partie 4 — Distances

| Étape | Que faire |
|---|---|
| **4.1–4.2** | Calcul, pour chaque cellule, de la distance (en µm) à la cellule la plus proche de chaque type, dans la même région. Les types très rares (< 0,01 %) sont écartés. |
| **4.3 Carte de proximité** 🖐️ | Chaque cellule est colorée selon sa distance au macrophage le plus proche (plafonnée à 200 µm). Choisissez la région, les types cibles (`PROX_TARGETS`) et la fenêtre. |
| **4.4 Surfaces 3D** 🔍 | Proximité des macrophages en fonction du phénotype des cellules tumorales (E-cadhérine/Vimentine, K8/K14…). La version interactive peut être tournée à la souris. |
| **4.5 Proches vs éloignées** 🖐️ | Expression de Ki67, Vimentine et HLA-DR dans les cellules tumorales proches (< `DIST_THRESHOLD_UM`, 50 µm par défaut) ou éloignées d'un macrophage, par patient (test t apparié). |
| **4.6 Distances moyennes** | Heatmap de la distance moyenne des cellules tumorales à chaque type du micro-environnement. |

### Partie 5 — Réseaux

Chaque type cellulaire est un nœud, dont la taille dépend de son abondance ; plus deux types sont proches, plus le trait qui les relie est épais. Les liens de plus de 250 µm ne sont pas tracés. Si un graphe est illisible, changez `NETWORK_SEED`. La section 5.2 compare les distances depuis un type de départ (`FROM_TYPE`) entre deux sous-ensembles (`GROUP1`, `GROUP2`), par exemple LUM primaire vs LUM métastase.

---

## 6. Fichiers produits

Tous les résultats sont écrits dans **`Results_spatialdata/`**, à côté du tableau des ROI (ou dans la ROI unique, ou dans `PROJECT_DIR` si vous l'avez renseigné).

| Fichier | Contenu |
|---|---|
| `Spillover_QC.pdf` | Contrôle de la correction du débordement (paire de marqueurs exclusifs avant/après) |
| `Diagnostics.pdf` | Nombre de cellules et expression moyenne par région |
| `Clusteringheatmap_all.pdf` | Heatmap des clusters avant/après annotation |
| `merge_phenocycler.xlsx` | **Votre annotation** cluster → type cellulaire |
| `Clusteringheatmap_merged.pdf` | Profil des types cellulaires annotés |
| `Umaps.pdf`, `Umap_interactive.html` | UMAP (le `.html` s'ouvre dans un navigateur) |
| `Spatial_celltypes_*.pdf` | Cartes des types cellulaires |
| `Results_counts.csv`, `Results_props*.csv` | Comptes et proportions par échantillon (ouvrables dans Excel) |
| `Abundance_*.pdf`, `Abundance_*_LvT.csv` | Barres empilées, comparaisons et tests statistiques |
| `Cellgeometricfeatures*.pdf` | Morphologie des cellules |
| `Clusterheatmap_focused_*.pdf`, `Corr_*.pdf` | Heatmaps ciblées et corrélations |
| `NN_*` , `Results_NN_*.csv` | Voisinage des cellules tumorales |
| `Proximity_map_*.pdf` | Carte de distance aux macrophages |
| `Distance_3D_*`, `Expression_byProximity.pdf`, `Distance_MeanCKtoTME_*.pdf` | Analyses de distance |
| `Distance_Relationships.pdf` | Réseaux de distances |
| `Results_distances_*.csv`, `Distance_violinplots.pdf` | Comparaison des distances entre groupes |
| `cells.parquet` | Tableau complet cellule par cellule (pour une réanalyse) |

Le dossier **`spatialdata/`** contient un fichier `<région>.zarr` par ROI (par ex. `EXP112_ROI_9x7.zarr`) : images, segmentation, positions et table annotée. La table contient les trois mesures d'intensité (`raw` = FCS, `mask_mean`, `mask_spillover`). Ces fichiers peuvent être rouverts dans d'autres outils de l'écosystème scverse, par exemple **napari-spatialdata** pour explorer les images en local.

---

## 7. Relancer une analyse ou ajouter des régions

- **Changer un paramètre** (cofacteur, nombre de clusters…) : modifiez-le, puis réexécutez la cellule et toutes les suivantes.
- **Ajouter des ROI** : ajoutez une ligne par nouvelle ROI dans `regions_metadata.xlsx` (ou passez `include` à 1), puis relancez depuis le début. Seules les nouvelles régions sont converties en Zarr, et le clustering est refait sur l'ensemble : **il faudra réannoter les clusters**.
- **Remplacer une carte de labels ou ajouter des images** : la région concernée est reconvertie automatiquement. Pour forcer la reconversion de tout, mettez `REBUILD_ZARR = True`.
- **Reprendre une session interrompue** : réexécutez la partie 0 (rapide, les Zarr existent déjà). Les parties suivantes rechargent les résultats sauvegardés si nécessaire.

---

## 8. Problèmes fréquents

| Symptôme | Cause probable et solution |
|---|---|
| `conda` ou `jupyter` : « commande introuvable » | Vous n'êtes pas dans l'Anaconda Prompt, ou l'environnement n'est pas activé : tapez `conda activate spatial`. |
| `ModuleNotFoundError: No module named …` | Le notebook n'utilise pas le bon environnement. Dans JupyterLab, vérifiez en haut à droite que le noyau est « Python 3 » de l'environnement `spatial`. Sinon, relancez `jupyter lab` depuis le terminal où `spatial` est activé. |
| « Ni tableau Excel ni dossier ROI trouvé » | Chemin mal écrit dans `REGIONS_METADATA` ou `RAW_DIR`. Sous Windows, utilisez des `/` (`"D:/ELOISE/regions_metadata.xlsx"`) ; copiez le chemin depuis l'explorateur de fichiers et remplacez les `\` par des `/`. |
| Une ROI est « exclue » à l'étape 0.3 | Son FCS ou son `channelNames.txt` n'a pas été trouvé : vérifiez `FCS_SUBDIR` et `CHANNELS_FILE_REL`, ou renseignez `fcs_file` / `channels_file` pour cette ROI dans le tableau. |
| « nb TIF = 0 » ou « labels absent » pour une ROI | Vérifiez `IMAGES_SUBDIR` et `LABELS_FILE_REL`, ou renseignez `images_dir` / `labels_file` pour cette ROI dans le tableau. |
| « Identifiants de région en double » | Deux ROI ont le même nom d'expérience et de dossier : ajoutez une colonne `region` avec des noms distincts. |
| Peu de cellules appariées à la carte de labels (étape 0.4) | Vérifiez que la carte correspond bien à la même région et à la même résolution que les images. Si le décalage est d'un pixel, essayez `COORD_ONE_BASED = False`. |
| Les contours ne tombent pas sur les noyaux (étape 0.5) | Carte de labels d'une autre région ou d'une autre taille d'image. |
| « matrice trop petite » en partie 3 | Pas assez de cellules par condition : diminuez `MIN_CELLS_CONDITION`. |
| Les comparaisons sont « ignorées » | Mode « une seule ROI », ou la colonne `Phenotype` du tableau ne contient ni `LUM` ni `TNC`. |
| L'éditeur d'annotation ne s'affiche pas | Installez ipywidgets (`pip install ipywidgets` dans l'environnement), puis redémarrez JupyterLab. Vous pouvez aussi éditer `merge_phenocycler.xlsx` dans Excel. |
| Mémoire saturée / ordinateur très lent | Trop de cellules ou de régions pour la mémoire disponible : réduisez `UMAP_NCELLS`, analysez les régions par lots, ou utilisez une machine avec plus de mémoire. |
| Les graphiques interactifs (UMAP, 3D) restent vides | Rechargez la page du navigateur ; vérifiez que `plotly` est bien installé dans l'environnement. |

---

## 9. Limites et précautions d'interprétation

- **L'annotation des types cellulaires est une étape d'expert.** Tous les résultats en dépendent. Les suggestions automatiques ne sont qu'un point de départ.
- **Le cofacteur de transformation** influence le clustering : choisissez-le en regardant les histogrammes, puis gardez-le identique pour toutes les régions d'une même étude.
- **Les distances sont en µm**, calculées avec `PIXEL_SIZE_UM` : une erreur sur ce paramètre fausse toutes les distances et les seuils (50 µm, 250 µm).
- **Voisins et distances** sont calculés à partir du centre des cellules, et non de leurs contours. Ils sont toujours calculés au sein d'une même région.
- **Les tests statistiques** sont faits par patient ou par échantillon. Avec peu d'échantillons par groupe, un résultat « non significatif » ne signifie pas une absence de différence.
- **Différences avec le pipeline R d'origine** (conçu pour l'IMC) : panel différent, voisins recalculés, morphologie mesurée sur la segmentation, cofacteur adapté. Les résultats ne sont donc pas directement superposables à ceux de l'analyse IMC. Chaque différence est signalée par ⚠️ dans le notebook.
