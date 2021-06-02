#########
### this script is intented to analyse Nanostring data though NACHO
#########
##################
### loading required packages
##################
if (!require("BiocManager")) {
  install.packages("BiocManager", ask =FALSE)
  library("BiocManager")
}
if (!require("NACHO")) {
  install.packages("NACHO", ask =FALSE)
  library("NACHO")
}
if (!require("GEOquery")) {
  BiocManager::install("GEOquery", dependencies = TRUE)
  library("GEOquery")
}
if (!require("dplyr")) {
  BiocManager::install("dplyr", dependencies = TRUE)
  library("dplyr")
}
if (!require("tidyr")) {
  BiocManager::install("tidyr", dependencies = TRUE)
  library("tidyr")
}
if (!require("tibble")) {
  BiocManager::install("tibble", dependencies = TRUE)
  library("tibble")
}
if (!require("limma")) {
  BiocManager::install("limma", dependencies = TRUE)
  library("limma")
}
if (!require("reshape")) {
  BiocManager::install("reshape", dependencies = TRUE)
  library("reshape")
}
if (!require("BiocManager")) {
  install.packages("BiocManager", ask =FALSE)
  library("BiocManager")
}
if (!require("pheatmap")) {
  BiocManager::install("pheatmap", ask =FALSE)
  library("pheatmap")
}
if (!require("ggplot2")) {
  BiocManager::install("ggplot2", ask =FALSE)
  library("ggplot2")
}
if (!require("argparse")) {
  install.packages("argparse", ask =FALSE)
  library("argparse")
}
if (!require("ggfortify")) {
  install.packages("ggfortify", ask =FALSE)
  library("ggfortify")
}
if (!require("limma")) {
  BiocManager::install("limma", ask =FALSE)
  library("limma")
}
if (!require("smooth")) {
  install.packages("smooth", ask =FALSE)
  library("smooth")
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", ask =FALSE)
  library("RColorBrewer")
}
if (!require("plotrix")) {
  install.packages("plotrix", ask =FALSE)
  library("plotrix")
}
if (!require("datasets")) {
  install.packages("datasets", ask =FALSE)
  library("datasets")
}
if (!require("reshape2")) {
  install.packages("reshape2", ask =FALSE)
  library("reshape2")
}
if (!require("ggridges")) {
  install.packages("ggridges", ask =FALSE)
  library("ggridges")
}
if (!require("preprocessCore")) {
  BiocManager::install("preprocessCore", ask =FALSE)
  library("preprocessCore")
}
if (!require("affy")) {
  BiocManager::install("affy", ask =FALSE)
  library("affy")
}
if (!require("Rtsne")) {
  BiocManager::install("Rtsne", ask =FALSE)
  library("Rtsne")
}


################
### Data given by the user
#####################
myhousekeeping_Genes <- ""
my_NACHO_normalisation_method <- "GLM"
or "GEOM"
myIDcolname <- "Unique_ID"
code_path <- "/media/rmejia/mountme88/code/Ncounter_RCC_processing/"
code_path <- normalizePath(code_path)
annot_path <- "/media/rmejia/mountme88/Projects/Maja-covid/Data/ssheet_annot_log2.csv"
outputfolder <- "/media/rmejia/mountme88/Projects/Maja-covid/Results/Preprocessing_through_Log2/NachoNorm"
your_main_groups  <- "Tissue"
label <- "Log2_NACHO_4HK_GLM"

#############
## Required functions
##########
source( paste0( code_path,"/libraries/","Nacho2Matrix.R" ) )

#########
## Body of the program
#########
dir.create(outputfolder, recursive = TRUE)
outputfolder <- normalizePath(outputfolder)

# Manual # vignette( "NACHO-analysis" )
####
input_RCCs <- load_rcc(data_directory = "/media/rmejia/mountme88/Projects/Maja-covid/Data/Original_RCC_log2/",
                       ssheet_csv ="/media/rmejia/mountme88/Projects/Maja-covid/Data/ssheetlog2_csv.csv",
                       id_colname = "Unique_ID" )

annot <- read.table( file=annot_path , sep="\t", header=TRUE)

# input_RCCs_normalized <- normalise(nacho_object= input_RCCs,
#                                    housekeeping_norm = TRUE,
#                                    normalisation_method = "GLM",
#                                    housekeeping_genes = c("NRDE2","MRPS7","GUSB","PGK1"),
#                                    remove_outliers = TRUE)
# 
# "MRPS7"
? load_rcc
? normalise
input_RCCs_normalized <- normalise(nacho_object= input_RCCs,
                                   housekeeping_norm = TRUE,
                                   normalisation_method = my_NACHO_normalisation_method ,
                                   housekeeping_genes = c("OAZ1","PGK1","SDHA","MRPS7"),
                                   remove_outliers = TRUE)

# visualise( input_RCCs_normalized )
input_RCCs_DefaultNorm_Mat <- NachoNorm2matrix( input_RCCs )
input_RCCs_4HK_GLM_Norm_Mat <- NachoNorm2matrix( input_RCCs_normalized  )
input_RCCs_Original_Counts_Mat <- Nacho_Orig_count_2matrix( input_RCCs)

Matrix_norm <- input_RCCs_4HK_GLM_Norm_Mat

# Exploring log2 matrix 
source( paste0( code_path,"/libraries/","matrix_N_annotdf_2_melteddf.R" ) )
annot4NachoOrginalcounts <- annot[ which( annot[,myIDcolname] %in%  colnames( input_RCCs_Original_Counts_Mat )    ) , ]
melted_Matrix_nacho_Orig <- matrix_N_annotdf_2_melteddf( input_RCCs_Original_Counts_Mat , annot4NachoOrginalcounts  ) 

source( paste0( code_path,"/libraries/","PCA_box_density_tsnes_plots.R" ) )
annot4NachoOrginalcounts[ , your_main_groups ] <- as.factor( annot4NachoOrginalcounts[ , your_main_groups ] )
PCA_box_density_tsnes_plots( paste0(  outputfolder,"/PCA_2D" )  ,
                             input_RCCs_Original_Counts_Mat , annot4NachoOrginalcounts , melted_Matrix_nacho_Orig , "NACHO_Orig_counts_log2" , 3, your_main_groups  )

# Exploring the normalized matrix (Over log 2)
source( paste0( code_path,"/libraries/","matrix_N_annotdf_2_melteddf.R" ) )
annot4Nachonorm <- annot[ which( annot[,myIDcolname] %in%  colnames( Matrix_norm )    ) , ]
melted_Matrix_nacho_norm <- matrix_N_annotdf_2_melteddf( Matrix_norm , annot4Nachonorm ) 

source( paste0( code_path,"/libraries/","PCA_box_density_tsnes_plots.R" ) )
annot_4_plotting_pca <- annot4Nachonorm
annot_4_plotting_pca[ , your_main_groups ] <- as.factor( annot_4_plotting_pca[ , your_main_groups ] )
PCA_box_density_tsnes_plots( paste0(  outputfolder,"/PCA_2D" )  ,
                             Matrix_norm ,  annot_4_plotting_pca , melted_Matrix_nacho_norm , "NACHO_log2_then_HK4_GEO" , 3, your_main_groups  )

rownames(annot4Nachonorm) <- annot4Nachonorm[,myIDcolname]
annot4Nachonorm <- annot4Nachonorm[colnames(Matrix_norm), ]
colors_from_rainbow <- colors_4_plotDensities( annot4Nachonorm , your_main_groups )

heatmap(Matrix_norm, ColSideColors = colors_from_rainbow )

##### quantiles over raw
ExpMat_input_RCCs_Original_Counts_plus_quantiles <- normalizeQuantiles(input_RCCs_Original_Counts_Mat )
heatmap(ExpMat_input_RCCs_Original_Counts_plus_quantiles)

annot4NachoOrginalcounts_plusquantiles <- annot[ which( annot[,myIDcolname] %in%  colnames( ExpMat_input_RCCs_Original_Counts_plus_quantiles )    ) , ]
melted_Matrix_nacho_Orig_plusquantiles <- matrix_N_annotdf_2_melteddf( ExpMat_input_RCCs_Original_Counts_plus_quantiles , annot4NachoOrginalcounts_plusquantiles  ) 

annotOrginalcounts_plusquantiles_4_plotting_pca <- annot4NachoOrginalcounts_plusquantiles
annotOrginalcounts_plusquantiles_4_plotting_pca[ , your_main_groups ] <- as.factor( annotOrginalcounts_plusquantiles_4_plotting_pca[ , your_main_groups ] )
PCA_box_density_tsnes_plots( paste0(  outputfolder,"/PCA_2D" )  ,
                             ExpMat_input_RCCs_Original_Counts_plus_quantiles ,  annotOrginalcounts_plusquantiles_4_plotting_pca , melted_Matrix_nacho_Orig_plusquantiles  , "NACHO_log2_plus_quantiles" , 3, your_main_groups  )

rownames( annotOrginalcounts_plusquantiles_4_plotting_pca) <- annotOrginalcounts_plusquantiles_4_plotting_pca[, myIDcolname]
annotOrginalcounts_plusquantiles_4_plotting_pca <- annotOrginalcounts_plusquantiles_4_plotting_pca[ colnames(ExpMat_input_RCCs_Original_Counts_plus_quantiles ),]

colors_log2_plus_quantiles <- colors_4_plotDensities(annotOrginalcounts_plusquantiles_4_plotting_pca  , your_main_groups )
heatmap(ExpMat_input_RCCs_Original_Counts_plus_quantiles , ColSideColors = colors_log2_plus_quantiles )

# quantiles without outlier
colnames_no_Outlier <- setdiff( colnames(input_RCCs_Original_Counts_Mat) , "20210413_covidnma_01_01_log2.RCC")
annotOrginalcounts_plusquantiles_4_plotting_pca 

log2_No_Outlier <- input_RCCs_Original_Counts_Mat[,colnames_no_Outlier]
rownames( annotOrginalcounts_plusquantiles_4_plotting_pca) <- annotOrginalcounts_plusquantiles_4_plotting_pca[, myIDcolname]
annotOrginalcounts_plusquantiles_4_plotting_pca <- annotOrginalcounts_plusquantiles_4_plotting_pca[ colnames(ExpMat_input_RCCs_Original_Counts_plus_quantiles ),]




dev.off()






