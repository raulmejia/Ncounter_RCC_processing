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

#############
## Required functions
##########
relocate_row_names <- function(somedf,somecolname){
  rownames(somedf) <- somedf[, somecolname]
  resultdf <- somedf[ , - which( colnames(somedf) %in% somecolname) ]
  return(resultdf)
}

relocate_count_norm_samplesID <- function(anotherdf){
  resultdf <- as.data.frame(anotherdf[,"count_norm_vec"] )
  rownames( resultdf) <- rownames( anotherdf)
  colnames( resultdf ) <- as.character(anotherdf[1, "count_norm_samplesID"])
  return(resultdf)
}

NachoNorm2matrix <- function(NachoObj){
  GeneNames <- NachoObj$nacho$Name
  count_norm_vec <- NachoObj$nacho$Count_Norm
  count_norm_samplesID <- NachoObj$nacho[ , myIDcolname]
  Long_df_norm <- data.frame( GeneNames ,
                              count_norm_samplesID , 
                              count_norm_vec )
  Long_df_norm$count_norm_samplesID <- as.factor( Long_df_norm$count_norm_samplesID )
  List_dfs_norm <- split(Long_df_norm, Long_df_norm$count_norm_samplesID)
  List_dfs_norm_rownames <- lapply( List_dfs_norm, relocate_row_names, "GeneNames" )
  
  List_dfs_norm_rownanmes_colnames <- lapply(List_dfs_norm_rownames , relocate_count_norm_samplesID )
  
  result_mat <- do.call( cbind, List_dfs_norm_rownanmes_colnames)
  result_mat <- as.matrix( result_mat )
  return( result_mat)
}
################
### Data given by the user
#####################
myhousekeeping_Genes <- ""
my_normalisation_method <- ""
myIDcolname <- "Unique_ID"
code_path <- "/media/rmejia/mountme88/code/Ncounter_RCC_processing/"
code_path <- normalizePath(code_path)
annot_path <- "/media/rmejia/mountme88/Projects/Maja-covid/Data/ssheet_annot.csv"
outputfolder <- "/media/rmejia/mountme88/Projects/Maja-covid/Results/NachoNorm"
your_main_groups  <- "Tissue"
label <- "NACHO_4HK_GLM"


#########
## Body of the program
#########
dir.create(outputfolder, recursive = TRUE)
outputfolder <- normalizePath(outputfolder)

# Manual
vignette( "NACHO-analysis" )
####
input_RCCs <- load_rcc(data_directory = "/media/rmejia/mountme88/Projects/Maja-covid/Data/Original_RCC_files/",
                       ssheet_csv ="/media/rmejia/mountme88/Projects/Maja-covid/Data/ssheet_csv.csv",
                       id_colname = "Unique_ID" )

annot <- read.table( file=annot_path , sep="\t", header=TRUE)

head(input_RCCs$nacho$Count_Norm)
head(input_RCCs_normalized$nacho$Count_Norm)

# Delete one HK STK11IP
# Delete one kidney
# visualise( Maja_Covid)
#MyHKG <- setdiff( input_RCCs$housekeeping_genes, "STK11IP")
#MyHKG2 <- setdiff( input_RCCs$housekeeping_genes, c("STK11IP","ALAS1","NRDE2") )
?normalise
# Maja_Covid_normalized <- normalise(nacho_object=Maja_Covid,
#                                    housekeeping_norm = TRUE,
#                                    normalisation_method = "GEO",
#                                    housekeeping_genes = MyHKG2,
#                                    remove_outliers = TRUE)

# Maja_Covid_normalized <- normalise(nacho_object=Maja_Covid,
#                                    housekeeping_norm = TRUE,
#                                    normalisation_method = "GLM",
#                                    housekeeping_genes = MyHKG2,
#                                    remove_outliers = TRUE)

input_RCCs_normalized <- normalise(nacho_object= input_RCCs,
                                   housekeeping_norm = TRUE,
                                   normalisation_method = "GLM",
                                   housekeeping_genes = c("NRDE2","MRPS7","GUSB","PGK1"),
                                   remove_outliers = TRUE)

#visualise( input_RCCs_normalized )



input_RCCs_DefaultNorm_Mat <- NachoNorm2matrix( input_RCCs )


Matrix_norm <- as.matrix(df_norm)

# Exploring the normalized matrix

source( paste0( code_path,"/libraries/","matrix_N_annotdf_2_melteddf.R" ) )
annot4Nachonorm <- annot[ which( annot[,myIDcolname] %in%  colnames( Matrix_norm )    ) , ]
melted_Matrix_nacho_norm <- matrix_N_annotdf_2_melteddf( Matrix_norm , annot4Nachonorm ) 

source( paste0( code_path,"/libraries/","PCA_box_density_tsnes_plots.R" ) )
annot_4_plotting_pca <- annot4Nachonorm
annot_4_plotting_pca[ , your_main_groups ] <- as.factor( annot_4_plotting_pca[ , your_main_groups ] )
PCA_box_density_tsnes_plots( paste0(  outputfolder,"/PCA_2D" )  ,
                             Matrix_norm ,  annot_4_plotting_pca , melted_Matrix_nacho_norm , label , 3, your_main_groups  )

rownames(annot4Nachonorm) <- annot4Nachonorm[,myIDcolname]
annot4Nachonorm <- annot4Nachonorm[colnames(Matrix_norm), ]
colors_from_rainbow <- colors_4_plotDensities( annot4Nachonorm , your_main_groups )

heatmap(Matrix_norm, ColSideColors = colors_from_rainbow )


?heatmap

ColSideColors 



?heatmap


# Plotting
PCA_box_density_tsnes_plots( paste0(  outputfolder,"/PCA_2D" )  ,
                             mymatrix , annot_4_plotting_pca , meltedrawdata , paste0( label ) )







######

#######

data_directory <- file.path("/media/rmejia/mountme88/Projects/Maja-covid/Data/NACHO_test", "GSE70970", "Data")
dir.create(data_directory,recursive = TRUE)
gse <- getGEO("GSE70970")
targets <- pData(phenoData(gse[[1]]))
getGEOSuppFiles(GEO = "GSE70970", baseDir = "/media/rmejia/mountme88/Projects/Maja-covid/Data/NACHO_test")
untar(
  tarfile = file.path( "/media/rmejia/mountme88/Projects/Maja-covid/Data/NACHO_test", "GSE70970", "GSE70970_RAW.tar"), 
  exdir = data_directory
)
targets$IDFILE <- list.files(data_directory)                
GSE70970 <- load_rcc(data_directory, targets, id_colname = "IDFILE")


### GEt the phenotypes
selected_pheno <- GSE70970[["nacho"]] %>% 
  select(IDFILE, `age:ch1`, `gender:ch1`, `chemo:ch1`, `disease.event:ch1`) %>% 
  distinct() %>% 
  mutate_all(~ na_if(.x, "NA")) %>% 
  drop_na()

## get the normalised counts
expr_counts <- GSE70970[["nacho"]] %>% 
  filter(grepl("Endogenous", CodeClass)) %>% 
  select(IDFILE, Name, Count_Norm) %>% 
  pivot_wider(names_from = "Name", values_from = "Count_Norm") %>% 
  column_to_rownames("IDFILE") %>% 
  t()
## Select the phenotypes and counts
samples_kept <- intersect(selected_pheno[["IDFILE"]], colnames(expr_counts))
expr_counts <- expr_counts[, samples_kept]
selected_pheno <- filter(selected_pheno, IDFILE %in% !!samples_kept)

## Design matrix
design <- model.matrix(~ `disease.event:ch1`, selected_pheno)

# limma
results <- eBayes(lmFit(expr_counts, design))

topTable(results)
