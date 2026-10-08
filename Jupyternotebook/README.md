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

**Vos données PhenoCycler** (détaillées en section 3)
- Les images TIF, **une par canal**, telles que sorties par le processeur.
- Une carte de segmentation par région (fichier TIF).
- Un fichier FCS par région (intensités moyennes par cellule).
- Le fichier `channelNames.txt`.

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

### 3.1 Organisation du dossier

Créez un dossier de travail (par exemple `D:/PhenoCycler/EXP100`) organisé **exactement** ainsi :

```
EXP100/
├── images/
│   ├── reg001_cyc001_ch001_DAPI.tif
│   ├── reg001_cyc001_ch002_Blank.tif
│   ├── reg001_cyc005_ch002_Vimentin.tif
│   └── …                           (un TIF par canal et par région)
├── labels/
│   ├── reg001_labels.tif           (une carte de segmentation par région)
│   └── reg002_labels.tif
├── fcs/
│   ├── EXP100-…_reg001_compensated.fcs
│   └── EXP100-…_reg002_compensated.fcs
├── channelNames.txt
└── regions_metadata.xlsx           (facultatif, mais nécessaire pour comparer des groupes)
```

### 3.2 Règles de nommage importantes

| Fichier | Règle |
|---|---|
| Images | Le nom doit contenir `regXXX_cycYYY_chZZZ` (c'est le format de sortie du processeur PhenoCycler). La fin du nom (`_DAPI`, `_Vimentin`…) n'est pas utilisée : le marqueur est lu dans `channelNames.txt`. |
| FCS | Le nom doit contenir `regXXX`. **Une région = un FCS.** C'est la liste des FCS qui définit les régions analysées. |
| Labels | Le nom doit contenir `regXXX` pour être associé à la bonne région (par ex. renommez `labels_cy002_ch1.tif` en `reg001_labels_cy002_ch1.tif`). **Exception :** s'il n'y a qu'une seule région et un seul fichier de labels, l'association est automatique. |
| `channelNames.txt` | Un marqueur par ligne, dans l'ordre cycle 1 canal 1, cycle 1 canal 2, … (4 canaux par cycle par défaut). C'est le fichier produit par le PhenoCycler. |

**Il n'est pas obligatoire d'avoir tous les canaux en images.** Les analyses utilisent les intensités du FCS ; les images servent à la visualisation. Vous pouvez donc ne mettre que DAPI et quelques marqueurs d'intérêt pour gagner de la place.

**Les labels sont facultatifs.** Sans carte de segmentation, chaque cellule est représentée par un cercle centré sur sa position, et les mesures de forme (aire, grand axe) sont approximées.

### 3.3 La carte de segmentation (labels)

Deux types de fichiers sont acceptés et reconnus automatiquement :

- **Carte d'instances** (recommandé) : image 16 ou 32 bits où chaque cellule a son propre numéro.
- **Carte de contours colorés** : image 8 bits où les contours des cellules sont tracés avec quelques couleurs. C'est le cas du fichier `labels_cy002_ch1.tif` (11 couleurs pour 26 666 cellules). Le notebook **reconstruit** alors les cellules : chaque zone fermée par un contour devient une cellule. Il associe ensuite chaque zone à la cellule du FCS dont le centre tombe dedans. Sur les données de test, 26 662 cellules sur 26 666 ont été retrouvées.

Si votre logiciel de segmentation peut exporter une vraie carte d'instances, préférez-la : c'est plus fiable.

### 3.4 Le fichier `regions_metadata.xlsx` (description des échantillons)

Ce tableau Excel décrit chaque région. Il est **nécessaire pour les comparaisons** (LUM vs TNC, primaire vs métastase, organes, patients appariés). Sans lui, toutes les analyses tournent, mais sans comparaison entre groupes.

Une ligne par région, avec ces colonnes :

| Colonne | Contenu | Exemple |
|---|---|---|
| `region` | identifiant de la région (comme dans les noms de fichiers) | `reg001` |
| `sample_id` | nom de l'échantillon | `P_12` |
| `Case` | identifiant du patient | `SPC04` |
| `Anno` | site du prélèvement | `PBC`, `LUNG`, `BRAIN`, `LIVER`, `LN`… |
| `Tumor` | type histologique | `IDC`, `ILC` |
| `Phenotype` | sous-type tumoral | `LUM`, `TNC`, `HER2`, `CTRL` |
| `Met` | 0 = primaire, 1 = métastase | `0` |
| `Paired` | 1 si le patient a des prélèvements appariés | `1` |

> **Astuce :** lancez une première fois le notebook sans ce fichier. Il crée un modèle pré-rempli `regions_metadata_template.xlsx` dans le dossier de résultats. Complétez-le dans Excel, enregistrez-le sous le nom `regions_metadata.xlsx` dans votre dossier de données, puis relancez à partir de la partie 1.

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
| **0.2 Paramètres** 🖐️ | Indiquez le chemin de votre dossier dans `RAW_DIR` (sous Windows, écrivez-le avec des `/`, par ex. `"D:/PhenoCycler/EXP100"`). Vérifiez la **taille de pixel** `PIXEL_SIZE_UM` : 0,5068 µm par défaut (PhenoCycler-Fusion, objectif 20×), à corriger selon votre acquisition. Toutes les distances en µm en dépendent. |
| **0.3 Canaux et régions** 🔍 | Vérifiez le tableau des canaux (marqueur ↔ cycle/canal) et le tableau des régions : chaque région doit avoir son FCS, ses images et ses labels. |
| **0.4 Conversion Zarr** | Chaque région est convertie **une seule fois** au format Zarr (dossier `spatialdata/`). Comptez quelques secondes à quelques minutes par région selon le nombre de canaux. Le tableau affiché indique le nombre de cellules retrouvées dans la carte de labels. |
| **0.5 Contrôle des labels** 🔍 | Un zoom montre le DAPI à côté des contours des cellules reconstruites. **Vérifiez que les contours entourent bien les noyaux.** Vous pouvez déplacer la fenêtre avec `QC_X0`, `QC_Y0` et `QC_SIZE` (en pixels). |

### Partie 1 — Types cellulaires

| Étape | Que faire |
|---|---|
| **0.6 Liste des types cellulaires** 🖐️ | `CLUSTERLEVELS` contient les noms possibles des types cellulaires. Elle est adaptée au panel PhenoCycler de test (CD31 → `Endo`, pas de marqueur pour les granulocytes ni les cellules dendritiques). **Adaptez-la à votre panel.** Les familles (CK, stroma, lymphocytes, myéloïdes) sont déduites des noms : un nom contenant `CK` est tumoral, `Mac` un macrophage, `Str` ou `Endo` du stroma. |
| **1.1 Métadonnées** | Chargement de `regions_metadata.xlsx`, ou création du modèle à compléter. |
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

Tous les résultats sont écrits dans **`Results_spatialdata/`**, à l'intérieur de votre dossier de données.

| Fichier | Contenu |
|---|---|
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

Le dossier **`spatialdata/`** contient un fichier `regXXX.zarr` par région : images, segmentation, positions et table annotée. Ces fichiers peuvent être rouverts dans d'autres outils de l'écosystème scverse, par exemple **napari-spatialdata** pour explorer les images en local.

---

## 7. Relancer une analyse ou ajouter des régions

- **Changer un paramètre** (cofacteur, nombre de clusters…) : modifiez-le, puis réexécutez la cellule et toutes les suivantes.
- **Ajouter des régions** : déposez les nouveaux fichiers (FCS, images, labels) dans les dossiers, ajoutez les lignes correspondantes dans `regions_metadata.xlsx`, puis relancez depuis le début. Seules les nouvelles régions sont converties en Zarr, et le clustering est refait sur l'ensemble : **il faudra réannoter les clusters**.
- **Remplacer une carte de labels ou ajouter des images** : la région concernée est reconvertie automatiquement. Pour forcer la reconversion de tout, mettez `REBUILD_ZARR = True`.
- **Reprendre une session interrompue** : réexécutez la partie 0 (rapide, les Zarr existent déjà). Les parties suivantes rechargent les résultats sauvegardés si nécessaire.

---

## 8. Problèmes fréquents

| Symptôme | Cause probable et solution |
|---|---|
| `conda` ou `jupyter` : « commande introuvable » | Vous n'êtes pas dans l'Anaconda Prompt, ou l'environnement n'est pas activé : tapez `conda activate spatial`. |
| `ModuleNotFoundError: No module named …` | Le notebook n'utilise pas le bon environnement. Dans JupyterLab, vérifiez en haut à droite que le noyau est « Python 3 » de l'environnement `spatial`. Sinon, relancez `jupyter lab` depuis le terminal où `spatial` est activé. |
| `RAW_DIR introuvable` | Chemin mal écrit. Sous Windows, utilisez des `/` (`"D:/PhenoCycler/EXP100"`). Copiez le chemin depuis l'explorateur de fichiers et remplacez les `\` par des `/`. |
| Une région n'a pas de labels dans le tableau 0.3 | Le nom du fichier de labels ne contient pas `regXXX` : renommez-le. |
| Peu de cellules appariées à la carte de labels (étape 0.4) | Vérifiez que la carte correspond bien à la même région et à la même résolution que les images. Si le décalage est d'un pixel, essayez `COORD_ONE_BASED = False`. |
| Les contours ne tombent pas sur les noyaux (étape 0.5) | Carte de labels d'une autre région ou d'une autre taille d'image. |
| « matrice trop petite » en partie 3 | Pas assez de cellules par condition : diminuez `MIN_CELLS_CONDITION`. |
| Les comparaisons sont « ignorées » | `regions_metadata.xlsx` absent, ou la colonne `Phenotype` ne contient ni `LUM` ni `TNC`. |
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
