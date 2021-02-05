###################################################################
#
#   Set of functions for formulating SVRs and inferring 
#    a co-splicing module network using WGCNA
#
###################################################################

# Load required R libraries
library(WGCNA)
library(stringr)  
library(ComplexHeatmap)
library(circlize)


# Function to compute splice variant regions from non-redundant LSV sets
make_SVRs <- function(exprMat, lsvDict) {
  
  # For function below, search redundant LSVs to see if it is in PSI matrix
  search_redundant_lsvs <- function (redundant_lsv, exprMat) {
    if (redundant_lsv != 'NaN') {
      # Check if there's multiple pairs
      # If not, see if it's in exprMat and return it
      if (!grepl('|', redundant_lsv, fixed = TRUE)) {
        if (redundant_lsv %in% colnames(exprMat)) {
          return(as.character(redundant_lsv))
        }
      } 
      # If so, we need to parse the string until we find one that is present
      else {
        split_lsvs <- strsplit(as.character(redundant_lsv), '|', fixed = TRUE)[[1]]
        for (i in 1:length(split_lsvs)) {
          lsv <- split_lsvs[i]
          if (lsv %in% colnames(exprMat)) {
            return(as.character(lsv))
          }
        }
      }
    }  
    return(NaN)
  }
  
  # Select LSVs that were initially labeled as non-redundant (python script)
  #   --> Initially selected based on PSI sample variance from MAJIQ
  nonRedundant <- rownames(lsvDict[lsvDict$IS_REDUNDANT == 'keep',])
  
  # If any non-redundant LSVs are absent due to coverage, check if there is a 
  # redundant mate that is present in PSI matrix
  lsv_set <- c()
  for (i in 1:length(nonRedundant)) {
    lsv <- nonRedundant[i]
    # If it's present, use it
    if (lsv %in% colnames(exprMat)) {
      lsv_set <- c(lsv_set, lsv) 
    }
    # If not, see if it exists in redundant pairing
    else {
      redundant_lsv <- lsvDict[lsv,]$REDUNDANT_PAIR
      search_lsv <- search_redundant_lsvs(redundant_lsv, exprMat)
      if (!is.na(search_lsv)) {
        lsv_set <- c(lsv_set, search_lsv)
      }
    }
  }
  
  # Create new PSI matrix of selected LSVs
  #  --> Use unique in case any duplicate LSVs were selected as non-redundant
  useExpr <- exprMat[,unique(lsv_set)]
  
  # Create subset of lsv dictionary with row ordering matching PSI matrix
  lsvDict_set <- lsvDict[colnames(useExpr),]
  
  # Use WGCNA's eigengene function for convenient calculation of PC1
  clusterExpr <- moduleEigengenes(useExpr, as.character(lsvDict_set$SVR))$eigengenes
  names(clusterExpr) <- substring(names(clusterExpr), 3) # Remove "ME"
  
  return(clusterExpr)
}


# Function to merge close modules and get module color labels 
merge_and_get_colors <- function(dColors, exprMat, dist_thresh = 0.25, make_plot = FALSE){
 
  # Plot the current tree if make_plot = TRUE
  if (make_plot == TRUE){
    # Calculate eigengenes
    MEs <- moduleEigengenes(exprMat, colors = dColors)$eigengenes
    # Cluster module eigengenes
    METree = hclust(as.dist(1-cor(MEs)), method = "average")
    # Make plot
    sizeGrWindow(7, 6)
    plot(METree, main = "Co-splicing Module Dendrogram",
         xlab = "", sub = "")
    abline(h=dist_thresh, col = "red")
  }
  
  # Now merge and return color labels
  merge = mergeCloseModules(exprMat, dColors, cutHeight = dist_thresh, verbose = 0)
  return(merge$colors)
}


# Function to create heatmap using correlations between module set
module_network_heatmap <- function(modExpr, modColors = NULL, d_height= unit(2, "cm"), d_width = unit(10, "mm")) {
  
  # Compute module-module correlations and create dendrogram
  moduleCor = cor(modExpr, modExpr, use = "p")
  modTree <- hclust(as.dist(1 - moduleCor), method = "average") 
  # Make correlations a signed network
  moduleCor <- (moduleCor + 1) / 2
  
  #########################################################################
  #
  #   Calcuate top (with module sizes) and left module color annotations
  #
  #########################################################################
  isDark <- function(colr) { 
    if (sum( col2rgb(colr) * c(299, 587,114))/1000 < 123) {
      return("White")
    } else {
      return("Black")
    }
  }
  
  # Top colors with sizes
  modSizes <- data.frame(table(str_to_title(modColors)), row.names = 1)
  modSizes$Colors <- sapply(rownames(modSizes), FUN = isDark)
  modSizes <- modSizes[colnames(modExpr),]
  modSize_annot <- HeatmapAnnotation(ModuleSize = anno_text(modSizes$Freq, rot = 0, location = 0.5, just = "center",
                                                                   gp = gpar(fill = rownames(modSizes), col = modSizes$Colors,
                                                                             fontsize = 9),
                                                                   height = unit(.75, "cm")),
                                     show_annotation_name = TRUE, annotation_name_side = "left")
 
  # Left colors
  modColors_left <- rowAnnotation(ModuleSize = colnames(modExpr), 
                                col = list(ModuleSize = setNames(colnames(modExpr), colnames(modExpr))),
                                show_legend = FALSE, show_annotation_name = FALSE, simple_anno_size = unit(0.5, "cm"))
  
  
  # Create heatmap object
  ht <- Heatmap(t(moduleCor), cluster_columns = modTree, cluster_rows = modTree, 
                  top_annotation = modSize_annot, left_annotation = modColors_left,
                  column_dend_height = d_height, row_dend_width = d_width,
                  show_column_names = FALSE, show_row_names = FALSE,
                  name = "Edge Weight - Signed", col = colorRamp2(c(0, 0.5, 1), c("white", "yellow", "red")),
                  heatmap_legend_param = list(legend_direction = "horizontal",legend_width = unit(8, "cm"),
                                              legend_height = unit(10, "cm"),  title_position = "topcenter"))
  
  return_items = list()
  return_items[['ht']] <- ht
  return_items[['cor']] <- t(moduleCor)
  return(return_items)
}

