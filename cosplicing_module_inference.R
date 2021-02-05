#################################################################################
#
#   Example R script for formulating SVRs from sample LSV PSI data and 
#    performing de novo network inference and module identification
#
#   - Required files are created using create_psi_matrix.py with Python 3
#
#################################################################################

# Load R libraries and functions for SVR formulation and network inference
source ("cosplicing_helper_functions.R")


#############################################################
#
#   LOAD SAMPLE PSI MATRIX AND LSV ANNOTATION DICTIONARY
#
#############################################################

# Load LSV PSI matrix and place in dataframe
lsv_mat = data.frame(read.csv('test_out/psi_matrix.example_data.csv', header = TRUE, row.names = 1, check.names = FALSE), 
                     stringsAsFactors = FALSE, check.names = FALSE)
dim(lsv_mat)

# Load LSV dictionary and place in dataframe
lsv_dict <- data.frame(read.csv('test_out/lsv_data_dictionary.example_data.csv', header = TRUE, row.names = 1), 
                       stringsAsFactors = FALSE)
dim(lsv_dict)


###############################################
#
#   SELECT LSVS AND FORMULATE SVRS
#
###############################################

# Select LSVs that were sucessfully quantified by MAJIQ in all samples
lsv_mat.select <- lsv_mat[,colSums(is.na(lsv_mat)) == 0]
dim(lsv_mat.select)

# Now formulate SVRs using LSV to SVR assignments and place in new dataframe
svr_mat <- make_SVRs(lsv_mat.select, lsv_dict)
dim(svr_mat)


##############################################
#
#  INFER CO-SPLICING NETWORK OF SVRS
#
##############################################

# Pick soft-thresholding value based on scale-freeness (WGCNA)
sft = pickSoftThreshold(svr_mat, powerVector = 1:20, verbose = 5, networkType = 'signed', 
                        corFnc = 'bicor', RsquaredCut = 0.9)

# Compute adjacency matrix for SVR network (WGCNA)
svr.adj <- adjacency(svr_mat, power = sft$powerEstimate, corFnc = 'bicor', type = 'signed')

# Convert adjacency matrix to topological overlap matrix (WGCNA)
svr.tom <- TOMsimilarity(svr.adj)
rownames(svr.tom) <- rownames(svr.adj)
colnames(svr.tom) <- colnames(svr.adj)


###################################################################
#
#   IDENTIFY CO-SPLICING NETWORK MODULES
#
#   - Paper used hierarchical clustering followed by spectral
#      clustering to identify co-splicing modules. As an example
#      we just use hierarchical clustering for simplicity.
#
####################################################################

# Convert TOM to distance
dissTOM.svr = 1-svr.tom

# Create network tree from TOM distance
svrTree = hclust(as.dist(dissTOM.svr), method = "average")

# Get modules from tree using dynamic tree cutting
hcColors = labels2colors(cutreeDynamic(dendro = svrTree, distM = dissTOM.svr, deepSplit = 2, 
                                       pamRespectsDendro = FALSE, minClusterSize = 100))

# Merge close modules and use for final set
svrMods <- merge_and_get_colors(hcColors, svr_mat, 0.25, FALSE)


#############################################################
#
#   COMPUTE MODULE SPLICING VALUES & SAVE DATA TO FILE
#
#############################################################

# Calculate module splicing values
modSplicing = moduleEigengenes(svr_mat, colors = svrMods, excludeGrey = TRUE)$eigengenes
names(modSplicing) <- str_to_title(substring(names(modSplicing), 3)) # Remove "ME" from WGCNA colors

# Plot co-splicing module network
spliceNet <- module_network_heatmap(modSplicing, svrMods)
draw(spliceNet$ht, heatmap_legend_side = "bottom", column_title = 'Co-splicing Module Network')



