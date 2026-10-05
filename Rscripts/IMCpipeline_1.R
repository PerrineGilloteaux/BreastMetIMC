####PART 1 CLUSTERING

#Won Jin Ho 
#R version 4.2.2 Patched (2022-11-10 r83330)
#Platform: x86_64-pc-linux-gnu (64-bit)
#Running under: Ubuntu 18.04.5 LTS

rm(list = ls())

####READ and CLUSTER FUNCTIONS####
returnfcs <- function(FDR_cutoff=.05,
                      metaDataFile=NULL,
                      panelDataFile=NULL,
                      dataDirectory=NULL,
                      shape_conditions=NULL,
                      color_conditions=NULL){
  ## This function generates an fcs file, subtype_markers, colors and shapes for clustering 
  require(scales);require(readxl);require(dplyr);require(flowCore)
  ## Directory and metadatafile checking
  if(!dir.exists(dataDirectory)) {stop('ERR: cannot find data directory')}
  if(!file.exists(metaDataFile)) {stop('ERR: cannot find metadata.xlsx or .csv file')}
  ## Read-in metadata and clean
  ifelse(grepl(metaDataFile,pattern='.xls'),md <- read_excel(metaDataFile),md <- read.csv(metaDataFile,header = TRUE))#must be in xl format or csv
  md$Phenotype <- factor(md$Phenotype)
  md$Anno <- factor(md$Anno)
  md$Case <- factor(md$Case)
  md$Tumor <- factor(md$Tumor)
  md$Primary <- factor(md$Primary)
  md$Met <- factor(md$Met)
  md$Paired <- factor(md$Paired)
  md$TMA <- factor(md$TMA)
  rownames(md) = md$sample_id
  md$sample_id <- md$sample_id
  
  ## Make sure all files in metadata present in datadirectory
  if(!all(md$file_name %in% list.files(dataDirectory)[grep(list.files(dataDirectory),pattern = '.fcs')])){
    print(paste('ERR: not all filenames in metadata present in data folder - missing',md$file_name[!which(md$file_name %in% list.files(dataDirectory)[grep(list.files(dataDirectory),pattern = '.fcs')])],'Subsetting...'))
    md <- md[-c(!which(md$file_name %in% list.files(dataDirectory)[grep(list.files(dataDirectory),pattern = '.fcs')])),]
  }
 
  ## Read fcs into fcs_raw
  fcs_raw <- read.flowSet(paste0(dataDirectory,"/",md$file_name), transformation = FALSE, truncate_max_range = FALSE)
  panel <- read_excel(paste0(work,'/Config/panel.xlsx'))
  head(data.frame(panel))
  panel$Parameter <- gsub('-', '_', panel$Parameter)
  
  
  ## Export out the parameter/panel data from the flowFrame to edit
  ## use panel$Antigen to fix description in panel_fcs
  ## use metal+isotope as mapping between panel from xlsx and panel from the fcs files
  
  panel_fcs <- pData(parameters(fcs_raw[[1]]))
  panel_fcs$desc <- gsub('-', '_', panel_fcs$desc)
  panel_fcs$desc[is.na(panel_fcs$desc)] <- paste0('NA_',which(is.na(panel_fcs$desc)))  
  
  rownames(panel_fcs) = panel_fcs$name
  
  ## Replace desc with revised Name
  panel_fcs[panel$Parameter,]$desc<-panel$Name
  
  ## Replace parameter data in flowSet with edits
  pData(parameters(fcs_raw[[1]])) <- panel_fcs
  
  ## Assign objects to marker lists
  subtype_markers <- panel$Name[panel$Subtype == 1]
  functional_markers <- panel$Name[panel$Functional == 1]
  otherparameters <- panel$Name[panel$Other ==1]
  cluster_by <- panel$Name[panel$Cluster == 1]
  
  ## Check marker lists
  if(!all(subtype_markers %in% panel_fcs$desc)){stop('ERR: Not all subtype_markers in panel_fcs$desc')}
  if(!all(functional_markers %in% panel_fcs$desc)){stop('ERR: Not all functional_markers in panel_fcs$desc')}
  if(!all(otherparameters %in% panel_fcs$desc)){stop('ERR: Not all otherparameters in panel_fcs$desc')}
  if(!all(cluster_by %in% panel_fcs$desc)){stop('ERR: Not all cluster markers in panel_fcs$desc')}
  
  fcs <- fsApply(fcs_raw, function(x, cofactor = 0.8){
    colnames(x) <- panel_fcs$desc
    expr <- exprs(x)
    expr <- asinh(expr[, union(subtype_markers,functional_markers)] / cofactor)
    exprs(x) <- expr
    x
  })
  
  ## Save out the original rownames from the parameter list from the fcs flowFrame
  panellist <- rownames(pData(parameters(fcs[[1]])))
  
  ## Create dummy list to save all expression data
  exprTr_list<-c()
  
  ## Save arc sin transformed expression + spatial data for each flowFrame
  for(i in 1:length(md$file_name)){
    
    ## Expression data is combined with spatial parameter data
    
    exprRaw<-exprs(fcs_raw[[i]])
    
    colnames(exprRaw)<-panel_fcs$desc
    
    expr<-cbind(exprs(fcs[[i]])[, union(subtype_markers,functional_markers)],exprRaw[,otherparameters])
    
    ## Combine other (spatial) data with the protein data
    colnames(expr)<-c(colnames(exprs(fcs[[i]])),otherparameters)
    
    ## Filter out any event that is 95th percentile for BOTH CD29 and CD45 (antibody aggregates)
    ##expr<-expr[expr[,"CD29"] < quantile(expr[,"CD29"], probs=0.95) & expr[,"CD45"] < quantile(expr[,"CD45"], probs=0.95),]
    
    ## Add to list
    exprTr_list[[i]]<-expr
    
  }
  
  ## Create a flowSet based on the list of expression data
  fcs1<-flowSet(sapply(exprTr_list,flowFrame))
  
  ## Change parameter rownames
  panel_fcs1 <- pData(parameters(fcs1[[1]]))
  rownames(pData(parameters(fcs1[[1]]))) <- rownames(panel_fcs[panel_fcs$desc %in% pData(parameters(fcs1[[1]]))$name,])

  
  
  ###to scale every flowframe
  ## Save out the original rownames from the parameter list from the fcs flowFrame
  panellist <- rownames(pData(parameters(fcs[[1]])))
  
  ## Create dummy list to save all expression data
  exprTr_list<-c()
  
  ## Save arc sin transformed expression + spatial data for each flowFrame
  for(i in 1:length(md$file_name)){
    
    ## Expression data is combined with spatial parameter data
    
    expr<-exprs(fcs[[i]])
    
    expr<-t(scale(t(expr)))
    
    ## Add to list
    exprTr_list[[i]]<-expr
    
  }
  
  ## Create a flowSet based on the list of expression data
  fcs2<-flowSet(sapply(exprTr_list,flowFrame))
  
  
  
    
  ## Get sample ids
  sample_ids <- rep(md$sample_id, fsApply(fcs1, nrow))
  
  ## Return: 
  ## fcs (only marker expressions arcsin transformed), 
  ## fcs1 (arcsin transformed + spatial parameters), 
  ## fcs2 (scaled arcsin expr per flowframe)
  ## fcsraw (all raw data), and all marker/parameter lists
  return(list('fcs'=fcs,'fcs1'=fcs1,'fcs2'=fcs2,'fcsraw'=fcs_raw,'subtype_markers'=subtype_markers,'functional_markers'=functional_markers,'otherparameters'=otherparameters,'cluster_by'=cluster_by,'sample_ids'=sample_ids,'meta_data'=md))
}



clusterfcs <- function(fcs=output$fcs,
                       cluster_by = output$cluster_by,
                       seed=1234,plottitle='consensus_plots',
                       scaleoption=F,
                       numclusters=40){
  ## Cell population identification with FlowSOM and ConsensusClusterPlus
  require(dplyr);require(FlowSOM);require(ConsensusClusterPlus)
  set.seed(seed)
  som <- ReadInput(fcs, transform = FALSE, scale = scaleoption) %>% BuildSOM(colsToUse = cluster_by)
  
  ## Get the cell clustering into 100 SOM codes
  cell_clustering_som <- som$map$mapping[,1]
  
  ## Metaclustering into numclusters with ConsensusClusterPlus
  codes <- som$map$codes
  mc <- ConsensusClusterPlus(t(codes), maxK = numclusters, reps = 100,
                             pItem = 0.9, pFeature = 1, title = plottitle, 
                             plot = "png", clusterAlg = "hc", 
                             innerLinkage = "average", finalLinkage = "average",
                             distance = "euclidean", seed = 1234)
  
  ## Get cluster ids for each cell
  code_clustering <- mc[[numclusters]]$consensusClass#metaclusters consensus
  cell_clustering <- code_clustering[cell_clustering_som]#cell clustering from som
  return(list('code_clustering'=code_clustering,'cell_clustering'=cell_clustering,'metaclusters'=mc))
}


####DIAGNOSTIC HEATMAP FUNCTIONS ####
plot_clustering_heatmap_wrapper <- function(fcs, cell_clustering, nclusters=40,
                                            color_clusters=colorassigned, cluster_merging = NULL, 
                                            cluster_by=output$cluster_by,
                                            clusterMergeFile=NULL,
                                            fileName = 'clusteringheatmap.pdf'){
  require(matrixStats);require(dplyr);require(RColorBrewer);require(pheatmap);require(readxl);require(flowCore);require(scales);require(pals);require(ComplexHeatmap)
  ## Will output the heatmap object and print it 

  #get expression
  expr <- fsApply(fcs, exprs);expr <-expr[,cluster_by]
  ## Scale expression of all markers to values between 0 and 1
  rng <- colQuantiles(expr, probs = c(0.01, 0.99))
  expr01 <- t((t(expr) - rng[, 1]) / (rng[, 2] - rng[, 1]))
  expr01[expr01 < 0] <- 0; expr01[expr01 > 1] <- 1;expr01 <-expr01[,cluster_by]
  
  
  ## Calculate the mean expression##################################################
  
  expr_mean <- data.frame(expr, cell_clustering = cell_clustering, check.names = FALSE) %>%
    group_by(cell_clustering) %>% summarize_all(funs(mean))
  expr01_mean <- data.frame(expr01, cell_clustering = cell_clustering, check.names = FALSE) %>%
    group_by(cell_clustering) %>% summarize_all(funs(mean))
  
  ## Calculate cluster frequencies
  
  clustering_table <- as.numeric(table(cell_clustering))
  clustering_prop <- round(clustering_table / sum(clustering_table) * 100, 2)
  
  ## Sort the cell clusters with hierarchical clustering
  
  d <- dist(expr_mean[, colnames(expr)], method = "euclidean")
  cluster_rows <- hclust(d, method = "average")
  
  expr_heat <- as.matrix(expr01_mean[, colnames(expr01)])
  rownames(expr_heat) <- expr01_mean$cell_clustering
  
  ## Annotation for the original clusters
  
  annotation_row <- data.frame(Cluster = factor(expr01_mean$cell_clustering))
  rownames(annotation_row) <- rownames(expr_heat)
  color_clusters1 <- color_clusters[1:nlevels(annotation_row$Cluster)]
  names(color_clusters1) <- levels(annotation_row$Cluster)
  annotation_colors <- list(Cluster = color_clusters1)
  
  ## Annotation for the merged clusters
  
  if(!is.null(clusterMergeFile)){
    ifelse(grepl(clusterMergeFile,pattern='.xls'),cluster_merging <- read_excel(clusterMergeFile),cluster_merging <- read.csv(clusterMergeFile,header = TRUE))
    cluster_merging$new_cluster <- factor(cluster_merging$new_cluster)
    annotation_row$Merged <- cluster_merging$new_cluster
    color_clusters2 <- color_clusters[1:nlevels(cluster_merging$new_cluster)]
    names(color_clusters2) <- levels(cluster_merging$new_cluster)
    annotation_colors$Merged <- color_clusters2
  }
  
  ## Colors for the heatmap
  
  legend_breaks = seq(from = 0, to = 1, by = 0.2)
  
  clusternames<-sort(unique(cluster_merging$new_cluster))
  
  colorassigned<-kovesi.rainbow_bgyrm_35_85_c69(length(unique(cluster_merging$new_cluster)))

  names(colorassigned)<-clusternames
  
  color_list = list(clusters=colorassigned)
  
  color_list_byoriginal = colorassigned[match(cluster_merging$new_cluster,names(colorassigned))]
  
  cp<-rowAnnotation(clusters=cluster_merging$new_cluster,
                    col=color_list,
                    gp = gpar(col = "white", lwd = .5),
                    prop=anno_barplot(
                      clustering_prop, 
                      gp = gpar(fill=color_list_byoriginal, col=F),
                      border = F,
                      bar_width = 0.75, 
                      width = unit(2,"cm")))
  
  q <- Heatmap(expr_heat, name="scaled",
          col=rev(brewer.rdbu(100)),
          row_order = cluster_merging[order(cluster_merging$new_cluster),]$original_cluster,
          cluster_columns = T,
          cluster_rows = T,
          border = NA,
          rect_gp = gpar(col = "white", lwd = .5),
          right_annotation = cp,
          show_row_names = T,
          row_names_gp = gpar(fontsize=7),
          column_names_gp = gpar(fontsize=10),
          heatmap_legend_param = list(at=seq(from = 0, to = 1, by = 0.2)),
          width = unit(10, "cm"))
  
  print('Colors:')
  print(color_clusters)
  
  pdf(fileName, width=8, height=6) 
  
  return(q)
  
  dev.off() 
  
}


plot_clustering_heatmap_wrapper2 <- function(fcs, cell_clustering, nclusters=40,
                                             color_clusters=colorassigned,
                                             colorbar=rev(brewer.rdbu(100)),
                                             subtype_markers,
                                             fileName = 'clusteringheatmap.pdf'){
  require(matrixStats);require(dplyr);require(RColorBrewer);require(pheatmap);require(readxl);require(flowCore);require(scales)
  ## Will output the heatmap object and print it 
  #get expression
  expr <- fsApply(fcs, exprs);expr <-expr[,subtype_markers]
  ## Scale expression of all markers to values between 0 and 1
  rng <- colQuantiles(expr, probs = c(0.01, 0.99))
  expr01 <- t((t(expr) - rng[, 1]) / (rng[, 2] - rng[, 1]))
  expr01[expr01 < 0] <- 0; expr01[expr01 > 1] <- 1;expr01 <-expr01[,subtype_markers]
  ## Calculate the mean expression##################################################
  pdf(fileName, width=8, height=11) 
  expr_mean <- data.frame(expr, cell_clustering = cell_clustering, check.names = FALSE) %>%
    group_by(cell_clustering) %>% summarize_all(funs(mean))
  expr01_mean <- data.frame(expr01, cell_clustering = cell_clustering, check.names = FALSE) %>%
    group_by(cell_clustering) %>% summarize_all(funs(mean))
  
  ## Calculate cluster frequencies
  
  clustering_table <- as.numeric(table(cell_clustering))
  clustering_prop <- round(clustering_table / sum(clustering_table) * 100, 2)
  
  ## Sort the cell clusters with hierarchical clustering
  
  d <- dist(expr_mean[, colnames(expr)], method = "euclidean")
  cluster_rows <- hclust(d, method = "average")
  
  expr_heat <- as.matrix(expr01_mean[, colnames(expr01)])
  rownames(expr_heat) <- expr01_mean$cell_clustering
  
  ## Colors for the heatmap
  
  #legend_breaks = seq(from = 0, to = 1, by = 0.2)
  #labels_row <- expr01_mean$cell_clustering
  
  labels_row <- paste0(expr01_mean$cell_clustering, " ")
  
  ## Annotation for the original clusters
  
  annotation_row <- data.frame(Cluster = factor(expr01_mean$cell_clustering))
  rownames(annotation_row) <- rownames(expr_heat)
  color_clusters1 <- color_clusters[1:nlevels(annotation_row$Cluster)]
  names(color_clusters1) <- levels(annotation_row$Cluster)
  annotation_colors <- list(Cluster = color_clusters1)
  
  p <- pheatmap(expr_heat, 
                color = colorbar, 
                cluster_cols = F,
                cluster_rows = F, 
                labels_row = labels_row,
                #scale="row",
                display_numbers = F, 
                number_color = "black",
                fontsize = 9, fontsize_number = 6,  
                #legend_breaks = legend_breaks,
                annotation_row = annotation_row, 
                annotation_colors = annotation_colors,
                cellwidth = 8,
                cellheight = 8,
                border_color = "white",
                annotation_legend = F
               )
  dev.off() 
  print('Colors:')
  print(color_clusters)
  print(p);return(p)
}



####UMAP####
#separate UMAP also created in Giotto
do_umap <- function(fcs,subtype_markers,sample_ids,cell_clustering,metadata,
                    clusterMergeFile='~/Desktop/ViralHCC/ViralHCC_merging.xlsx',
                    seed = 1234, ncells=2000,sample_subset=NULL){
  require(umap);require(flowCore);require(readxl)
  expr <- fsApply(fcs, exprs);expr <-expr[,subtype_markers]
  ## Create vector to later find and skip duplicates
  dups <- duplicated(expr[, subtype_markers])
  dups <- which(!(dups))## Find and skip duplicates
  ifelse(grepl(clusterMergeFile,pattern='.xls'),cluster_merging <- read_excel(clusterMergeFile),cluster_merging <- read.csv(clusterMergeFile,header = TRUE))
  ## New clustering1m
  mm <- match(cell_clustering, cluster_merging$original_cluster)
  cell_clustering1m <- cluster_merging$new_cluster[mm]
  ## Create a data frame of sample_ids and cell_clustering1m
  dtf<-data.frame(ids=sample_ids,type=cell_clustering1m)
  #dtf$B<- dtf$type!="B"#add a column that indicates whether the cell is a B cell or not; TRUE is non-B
  ##Why exclude B CELLS?
  ## WE HAVE NO B CELLS bc dtf$B depends on dtf$type depends on cellclustering1m which is just numbers 1:40 so..?
  #should it be the named parts cluster in merge file corresponding to it like if 30 -> grepl(cluster_merging[30,3],pattern='^B')??
  #Im blocking this out till we know why we have to do this
  ## Create subset columns without B cells (ONLY to generate the correct # of columns in inds2 object)
  #sample_ids2 <- dtf$ids[dtf$type!="B"] #sampleids without B cells
  #cell_clustering1m2 <- dtf$type[dtf$type!="B"] #clusters without B cells
  ## Data subsampling: create indices by sample
  inds <- split(1:length(sample_ids), sample_ids) #to get original indexes belonging to each cluster
  #inds2 <- split(1:length(sample_ids2), sample_ids2) #to fill in the original indexes that do not have B cells
  samplenames <- names(inds) #create a name vector of the files
  #FOR THIS TO WORK MUST BE IN FORMAT PRE/POST Tx
  # for (i in 1:(length(samplenames)/2)){#1:15 was here because ids in dtf was 30 factors and could be divided into 0 and 6wk for each so changed it
  #   templength <- length(inds2[[i]])
  #   inds2[[i]] <- inds[[i]][dtf$B[dtf$ids==samplenames[i]]] #subsets the "B cell or not a B cell" column for each sample by the name
  #   inds2[[i]] <- inds2[[i]][1:templength]
  # }
  
  custom.settings = umap.defaults
  custom.settings$seed = seed
  
  #custom.settings$n.neighbors = neighbors
  ####umapindex generation####
  #umap ncells = table of sample ids with how many to downsample to by default col = id, row = ncells
  #sample ids = chr [1:2851129] "4927_0wk" "4927_0wk" "4927_0wk" "4927_0wk" ...
  #^ from ## Generate sample IDs corresponding to each cell in the 'expr' matrix sample_ids <- rep(md$sample_id, fsApply(fcs_raw, nrow))
  #can subset sample_ids and rerun umap 
  #can do timepoint or patient number or pre/post if you know corresp sample ids
  #sample_subset = '02_0wk' or c('02_0wk','2973_0wk') for example picking all 0 wk ones or use regex sample_ids[(grepl(sample_ids,pattern = '0wk'))]
  ifelse(is.null(sample_subset),
         umap_ncells <- pmin(table(sample_ids), ncells),
         umap_ncells <- pmin(table(sample_ids), ncells)[sample_subset]
  )
  if(!is.null(sample_subset)){inds <- inds[sample_subset]}
  umap_inds <- lapply(names(inds), function(i){
    s <- sample(inds[[i]], umap_ncells[i], replace = FALSE)
    intersect(s, dups)
  })
  set.seed(seed)
  umap_inds <- unlist(umap_inds)
  umap_out <- umap(expr[umap_inds, subtype_markers], config = custom.settings, method = 'naive')
  umapRes2D = data.frame(umap1 = umap_out$layout[, 1], umap2 = umap_out$layout[, 2], 
                         expr[umap_inds, subtype_markers],
                         sample_id = sample_ids[umap_inds], cell_clustering = factor(cell_clustering1m[umap_inds]), check.names = FALSE)
  
  #exclude any unassigned cluster post umap if needed--this has to be done by looking at the two columns 
  #metadata$sampleid is just a number in this metadatafile so to make unique ones combine samp_id col with timepoint (0wk)
  #to get in format of umapres2d$sample_id which looks like "02_0wk" do:
  return(umapRes2D)
}


plotUmap <- function(umapRes,seed=1234,neighbors=10,midpoint,color_clusters='auto',code_clustering,subtype_markers=NULL)
{require(umap);require(ggplot2);require(viridis);require(ggrepel)
  if(length(color_clusters) == 1 && color_clusters =='auto'){color_clusters <- hue_pal()(length(unique(code_clustering)))}
  custom.settings = umap.defaults
  custom.settings$seed = seed
  custom.settings$n.neighbors = neighbors
  ggp <- ggplot(umapRes,  aes(x = umap1, y = umap2, color = cell_clustering)) +
    geom_point(size = 1) +
    theme_bw() +
    
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
    ) +
    
    scale_color_manual(values = color_clusters, name="CLUSTERS") +
    guides(color = guide_legend(override.aes = list(size = 3), ncol = 2))
  
  print(ggp)
  
  print(ggp + facet_wrap(~ type, ncol=4)+ggtitle('TYPE'))
  print(ggp + facet_wrap(~ site, ncol=4)+ggtitle('SITE'))
  print(ggp + facet_wrap(~ SPC, ncol=4)+ggtitle('SPC ID'))
  print(ggp + facet_wrap(~ tumor, ncol=4)+ggtitle('TUMOR'))
  print(ggp + facet_wrap(~ TMA, ncol=3)+ggtitle('TMA'))

  ggp2 <- ggplot(umapRes,  aes(x = umap1, y = umap2, color = sample_id)) +
    geom_point(size = 1) +
    theme_bw() +
    
    theme(panel.grid.major = element_blank(),  
          panel.grid.minor = element_blank()
    ) +
    
    guides(color = guide_legend(override.aes = list(size = 3), ncol = 2))
  
  print(ggp2)
  #can specify which markers to display
  if(!is.null(subtype_markers)){
    for(i in subtype_markers)
    {
      ggp <- ggplot(umapRes,  aes(x = umap1, y = umap2, color = umapRes[,i])) +
        geom_point(size = 1) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()) +
        scale_color_gradient2(i, low="dark blue",mid="white",high="dark red", midpoint = mean(unlist(umapRes[,i])))
      print(ggp)
    }
  }
}





####REQUIRED LIBRARIES####
library(readxl)
library(stringr)
library(matrixcalc)
library(matrixStats)
library(Hmisc)
library(reshape2)
library(dplyr)
library(plotrix)
library(multcomp)
library(flowCore)
library(sf)
library(clusterSim)
library(limma)
library(corrplot)
library(packcircles)
library(cowplot)
library(autoimage)
library(ggplot2)
library(ggpubr)
library(ggiraphExtra)
library(ggvoronoi)
library(ggridges)
library(gridExtra)
library(igraph)
library(qgraph)
library(circlize)
library(scales)
library(RColorBrewer)
library(ComplexHeatmap)
library(pheatmap)
library(pals)
library(plot3D)
library(akima)
library(basetheme)



#======================
####RUNNING DATA####
#======================


####DATA LOADING####

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

work<-getwd()

#.libPaths(paste0(work,"/Libraries"))



## If loading from prior run




output<-readRDS('..\\Input_Files\\backup_output.rds')


## Read (skip if previously ran)

output <- returnfcs(metaDataFile = paste0(work,"..\\Input_Files\\metadata.xlsx"),
                    panelDataFile = paste0(work,"..\\Input_Files\\panel.xlsx"),
                    dataDirectory = paste0(work,"/Data"))


## Set up levels

samplevels=c(paste0("P_",1:342))

samplestoexclude <- c(output$meta_data$sample_id[output$meta_data$Tumor=="CTRL"]) # CTRL 

typelevels=c("CTRL","LUM","TNC","HER2")

sitelevels=c("PBC",
             "PBC_1YR",
             "PBC_tx",
             "PBC_recur",
             "GI",
             "PANC",
             "LIVER",
             "BRAIN",
             "OVARY",
             "PLEURA",
             "SPINE",
             "LUNG",
             "LN",
             "CTRL")

tumorlevels=c("IDC","ILC","IMC","CTRL")

TMAlevels=c("TMA788",
            "TMA789",
            "TMA801",
            "TMA974",
            "TMA975")

SPClevels=c("SPC01",
            "SPC02",
            "SPC03",
            "SPC04",
            "SPC06",
            "SPC07",
            "SPC09",
            "SPC10",
            "SPC11",
            "SPC12",
            "SPC13",
            "SPC14",
            "SPC15",
            "SPC16",
            "SPC17",
            "SPC18",
            "SPC19",
            "SPC20",
            "SPC21",
            "SPC22",
            "SPC23",
            "SPC24",
            "SPC25",
            "SPC26",
            "CTRL")


####DIAGNOSTICS####
## Spot check - number of cells per sample
cell_table <- table(output$sample_ids)
ggdf <- data.frame(sample_id = names(cell_table), 
                   cell_counts = as.numeric(cell_table))
ggdf$SPC <- factor(output$meta_data$Case[match(ggdf$sample_id,output$meta_data$sample_id)], levels = SPClevels)
ggdf$tumor <- factor(output$meta_data$Tumor[match(ggdf$sample_id,output$meta_data$sample_id)], levels=tumorlevels)
ggdf$type <- factor(output$meta_data$Phenotype[match(ggdf$sample_id,output$meta_data$sample_id)], levels=typelevels)
ggdf$site <- factor(output$meta_data$Anno[match(ggdf$sample_id,output$meta_data$sample_id)], levels=sitelevels)
ggdf$TMA <- factor(output$meta_data$TMA[match(ggdf$sample_id,output$meta_data$sample_id)], levels=TMAlevels)

ggp<-ggplot(ggdf, aes(x = sample_id, y = cell_counts, fill = type)) + 
  geom_bar(stat = 'identity') + 
  #geom_text(aes(label = cell_counts), angle = 45, hjust = 0.5, vjust = -0.5, size = 2) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=5)) +  
  scale_x_discrete(drop = FALSE)
pdf('Diagnostics_cellcounts.pdf',width=6, height=4);ggp; dev.off()

## Multi-dimensional scaling plot to show similarities between samples
## Get the mean marker expression per sample
expr_mean_sample_tbl <- data.frame(sample_id = output$sample_ids, fsApply(output$fcs,exprs)) %>%
  group_by(sample_id) %>%  summarize_all(funs(mean))
expr_mean_sample <- t(expr_mean_sample_tbl[, -1])
colnames(expr_mean_sample) <- expr_mean_sample_tbl$sample_id
mds <- plotMDS(expr_mean_sample, plot = FALSE)
ggdf <- data.frame(MDS1 = mds$x, MDS2 = mds$y,
                   sample_id = colnames(expr_mean_sample))
ggdf$SPC <- factor(output$meta_data$Case[match(ggdf$sample_id,output$meta_data$sample_id)], levels = SPClevels)
ggdf$tumor <- factor(output$meta_data$Tumor[match(ggdf$sample_id,output$meta_data$sample_id)], levels=tumorlevels)
ggdf$type <- factor(output$meta_data$Phenotype[match(ggdf$sample_id,output$meta_data$sample_id)], levels=typelevels)
ggdf$site <- factor(output$meta_data$Anno[match(ggdf$sample_id,output$meta_data$sample_id)], levels=sitelevels)
ggdf$TMA <- factor(output$meta_data$TMA[match(ggdf$sample_id,output$meta_data$sample_id)], levels=TMAlevels)
ggp<-ggplot(ggdf, aes(x = MDS1, y = MDS2, color = site, shape = type)) +
  geom_point(size = 2.5) +
  #geom_text(aes(label = patient_id)) +
  theme_bw()+
  theme(plot.background = element_rect(fill="black"),
        panel.background = element_rect(fill="black"),
        panel.grid = element_blank(),
        axis.line = element_line(color="white"),
        axis.text = element_text(color="white"),
        axis.title = element_text(color="white"),
        legend.background = element_rect(fill="black"),
        legend.key = element_rect(fill="white"),
        legend.text = element_text(color="white"))
pdf('Diagnostics_MDS_sample.pdf',width=6, height=6);ggp; dev.off()

## Colors for the heatmap
color <- colorRampPalette(brewer.pal(n = 9, name = "YlGnBu"))(100)
pdf('Diagnostics_Heatmap.pdf',width=8, height=20)
pheatmap(expr_mean_sample_tbl[,c(2:18,21:40)], color = color, display_numbers = TRUE,
         number_color = "black", fontsize_number = 3, fontsize_row = 3,
         cluster_rows = F,
         annotation_colors = annotation_colors, clustering_method = "average")
dev.off()


####CLUSTERING####

##Revised loading depending on the diagnostics if needed

output <- returnfcs(metaDataFile = paste0(work,"/Config/metadata.xlsx"),
                    panelDataFile = paste0(work,"/Config/panel.xlsx"),
                    dataDirectory = paste0(work,"/Data"))

##Clustering

output[(length(output)+1):(length(output)+3)] <- clusterfcs(fcs=output$fcs, numclusters=50, scaleoption = F) 
#output$fcs uses just arcsin transformed data; scaleoption scales the entire dataset
#output$fcs2 uses scaled expression data that is scaled for each flowFrame

names(output)[(length(output)-2):(length(output))] <- c('code_clustering','cell_clustering','metaclusters')


####ANNOTATIONS AND VISUALIZATION OF CLUSTERS####

## Load merge file
## Assign colors
clusterMergeFile = paste0(work,"/Config/merge.xlsx") #create dummy merger numbers prior to annotation
cluster_merging <- read_excel(clusterMergeFile)

clusterlevels=c(1:50)

colorassigned<-kovesi.rainbow_bgyrm_35_85_c69(length(unique(clusterlevels)))

#cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

clusternames<-clusterlevels
names(colorassigned)<-clusternames
mm1 <- match(output$cell_clustering, cluster_merging$original_cluster)
cell_clustering1m <- cluster_merging$new_cluster[mm1]
output$cell_clustering1m <- cell_clustering1m

## Metacluster heatmaps
plot_clustering_heatmap_wrapper(fcs=output$fcs,
                                color_clusters = kovesi.rainbow_bgyrm_35_85_c69(50),
                                cell_clustering = output$cell_clustering, 
                                cluster_by=output$cluster_by,
                                clusterMergeFile = clusterMergeFile,
                                fileName = 'Clusteringheatmap_all.pdf'); dev.off()

clusterMergeFile = paste0(work,"/Config/merge.xlsx") #create dummy merger numbers prior to annotation
cluster_merging <- read_excel(clusterMergeFile)

clusterlevels=c("Tc", "ThN", "ThEM", "Treg",
                "NK_T", "NK", "B_I", "B_II",
                "Gran", "DC",
                "Mac_I", "Mac_II", "Mac_III", "Mac_IV",
                "Str_PS", "Str_V", "Str_VC", "Str_VS",
                "CK8", "CK8_E", "CK8_V",
                "CK14", "CK14_E", "CK14_EV",
                "CK8_14_E", "CK8_14_EV",
                "CKlo", "CKlo_EV", "CKother",
                "Neuron", "NA")

clusterlevels_macs=c("CK8", "CK8_E", "CK8_V",
                     "CK14", "CK14_E", "CK14_EV",
                     "CK8_14_E", "CK8_14_EV",
                     "CKlo", "CKlo_EV", "CKother", 
                     "Neuron", "NA",
                     "Str_PS", "Str_V", "Str_VC", "Str_VS",
                     "NK_T", "NK", "B_I", "B_II",                     
                     "Tc", "ThN", "ThEM", "Treg", 
                     "Gran", "DC",
                     "Mac_I", "Mac_II", "Mac_III", "Mac_IV")

colorassigned<-kovesi.rainbow_bgyrm_35_85_c69(length(unique(cluster_merging$new_cluster)))
hex <- hue_pal()(9)
colorassignedbroad <- c(rep(hex[1],3), #ck8
                        rep(hex[2],3), #ck14
                        rep(hex[3],2), #ck8and14
                        rep(hex[4],3), #ckother
                        rep(hex[5],2), #other
                        rep(hex[6],4), #str
                        rep(hex[7],8), #lym
                        rep(hex[8],2), #myeloid1
                        rep(hex[9],4)) #myeloid2


clusternames<-clusterlevels
names(colorassigned)<-clusternames
mm1 <- match(output$cell_clustering, cluster_merging$original_cluster)
cell_clustering1m <- cluster_merging$new_cluster[mm1]
output$cell_clustering1m <- cell_clustering1m

plot_clustering_heatmap_wrapper2(fcs=output$fcs,
                                 colorbar = kovesi.diverging_bwr_40_95_c42(100),
                                 subtype_markers = c("CK8","CK14","PanCK","ECAD","VIM","COL","SMA","PDPN",
                                                     "CD3","CD8","CD4","CD20","CD57","CD45RA","CD45RO","FOXP3","TOX","LAG3","GZMB","PD1",
                                                     "CD33","CD14","CD15","CD16","CD68","CD163","CD206","DCSIGN","DCLAMP","HLADR","CD86","PDL1",
                                                     "ARG1","CLCASP3","KI67"),
                                color_clusters = colorassigned,
                                cell_clustering = factor(output$cell_clustering1m, levels=clusterlevels), 
                                fileName = 'Clusteringheatmap_merged.pdf');dev.off()

plot_clustering_heatmap_wrapper2(fcs=output$fcs,
                                 colorbar = kovesi.diverging_bwr_40_95_c42(100),
                                 subtype_markers = c("CK8","CK14","PanCK","ECAD","VIM","COL","SMA","PDPN",
                                                     "CD3","CD8","CD4","CD20","CD57","CD45RA","CD45RO","FOXP3","TOX","LAG3","GZMB","PD1",
                                                     "CD33","CD14","CD15","CD16","CD68","CD163","CD206","DCSIGN","DCLAMP","HLADR","CD86","PDL1",
                                                     "ARG1","CLCASP3","KI67","PTPN22"),
                                 color_clusters = colorassigned,
                                 cell_clustering = factor(output$cell_clustering1m, levels=clusterlevels), 
                                 fileName = 'Clusteringheatmap_merged2.pdf');dev.off()


## Compute umap
umapRes <- do_umap(fcs=output$fcs,
                   subtype_markers = output$subtype_markers,
                   sample_ids = output$sample_ids,
                   cell_clustering = output$cell_clustering, 
                   metadata=output$metadata,
                   clusterMergeFile=clusterMergeFile,
                   seed = 1234, 
                   ncells=500,
                   sample_subset=NULL)

mm <- match(as.character(umapRes$sample_id), as.character(output[["meta_data"]]$sample_id))
umapRes$sample_id <- factor(output[["meta_data"]]$sample_id[mm], levels=samplevels)
umapRes$SPC <- factor(output[["meta_data"]]$Case[mm], levels = SPClevels)
umapRes$tumor <- factor(output[["meta_data"]]$Tumor[mm], levels=tumorlevels)
umapRes$type <- factor(output[["meta_data"]]$Phenotype[mm], levels=typelevels)
umapRes$site <- factor(output[["meta_data"]]$Anno[mm], levels=sitelevels)
umapRes$TMA <- factor(output[["meta_data"]]$TMA[mm], levels=TMAlevels)

umapRes<-umapRes[umapRes$cell_clustering!="NA",]



pdf('Umaps.pdf',width=8,height=8)
plotUmap(umapRes = umapRes,
         code_clustering=cell_clustering1m,
         color_clusters = colorassigned[names(colorassigned)!="NA"],
         subtype_markers = output$subtype_markers)
dev.off()


## Save output list
saveRDS(output, file="backup_output.rds")
saveRDS(umapRes, file="backup_umap.rds")
umapRes<-readRDS('backup_umap.rds')
