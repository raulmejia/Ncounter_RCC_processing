#########
## Nacho_default_Normalization.R
#########
### This script receive a folder with RCC files and retreives:
###     A) Expression matrix of the InputCounts (as expressed in the NACHO object)
###     B) Expression matrix with the default normalization given by NACHO
###     C) Nacho Object
###
###     Document that describes Nacho normalization
###
### Your annotation file should contain a col with the name equal to the content of this parameter myIDcolname
### for example "Unique_ID" if the content is myIDcolname <- "Unique_ID"
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

############################# 
## Data given by the user
##############################
# create parser object
parser <- ArgumentParser()
# specify our desired options 
# by default ArgumentParser will add an help option 
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
                    help="Print extra output [default]")
parser$add_argument("-q", "--quietly", action="store_false", 
                    dest="verbose", help="Print little output")
parser$add_argument("-d", "--datadirectory", type="character", 
                    help="path to the directory with the RCC files")
parser$add_argument("-s", "--ssheet", type="character", 
                    help="path to your ssheet")
parser$add_argument("-i", "--idcolname", type="character", 
                    help="colanme used as id")
parser$add_argument("-a", "--annotation", type="character", 
                    help="Annotation File")
parser$add_argument("-n", "--normalizationmethod", type="character", 
                    help="normalization method only for housekeeping?")
parser$add_argument("-c", "--codepath", type="character", 
                    help="Path where you have your code")
parser$add_argument("-g", "--maingroups", type="character", 
                    help="The most reelevant column in your annotation file")
parser$add_argument("-l", "--label", type="character", 
                    help="label for your results")
parser$add_argument("-o", "--outputfolder", type="character", 
                    help="Folder to store the outputs")

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args <- parser$parse_args( )
# print some progress messages to stderr if "quietly" wasn't requested

###########################
## Reading the parameters
###########################
data_directory_path <- args$datadirectory
# data_directory_path <- "/media/rmejia/mountme88/Projects/Maja-covid/Data/Original_RCC_log2"
data_directory_path <- normalizePath( data_directory_path)

ssheet_path <- args$ssheet
#  ssheet_path <-"/media/rmejia/mountme88/Projects/Maja-covid/Data/ssheetlog2_csv.csv"

myIDcolname <- args$idcolname
# myIDcolname <- "Unique_ID"

annotation_path <- args$annotation
# annotation_path <- "/media/rmejia/mountme88/Projects/Maja-covid/Data/ssheet_annot_log2.tsv"

mynormalizationmethod <- args$normalizationmethod
# mynormalizationmethod <- "GEO" # it could be GLM as well

code_path <- args$codepath
# code_path <- "/media/rmejia/mountme88/code/Ncounter_RCC_processing/"
code_path <- normalizePath( code_path )

your_main_groups  <- args$maingroups
# your_main_groups  <- "Tissue"

label <- args$label
# label <- "Log2_NACHO_4HK_GLM"

outputfolder <- args$outputfolder
# outputfolder <- "/media/rmejia/mountme88/Projects/Maja-covid/Results/Preprocessing_through_Log2/NachoNorm"
dir.create(outputfolder, recursive = TRUE) ; outputfolder <- normalizePath( outputfolder )

#####################################
## Loading some required libraries  #
#####################################
source( paste0( code_path,"/libraries/","Nacho2Matrix.R" ) )

#########
## Body of the program
#########
# Manual # vignette( "NACHO-analysis" )
####
input_RCCs <- load_rcc(data_directory = data_directory_path ,
                       ssheet_csv = ssheet_path ,
                       id_colname = myIDcolname  )
annot <- read.table( file = annotation_path , sep="\t", header=TRUE)

#####
## Extracting the Input data (Counts) stored in the NACHO object
#####
input_RCCs_Original_Counts_Mat <- Nacho_Orig_count_2matrix( input_RCCs)

path2save_orig <- paste0( outputfolder , "/","ExpMat_as_input_from_the_RCCs_in_the_folder--",basename(data_directory_path),"--.tsv")
write.table( input_RCCs_Original_Counts_Mat, file=path2save_orig,  sep="\t", row.names = TRUE, col.names = TRUE)

# Saving the Annot for related to The Original matrix
path2save_annot_from_origExpMat <- paste0( outputfolder , "/","Annot_from_ExpMat_as_input_from_the_RCCs_in_the_folder--",basename(data_directory_path),"--.tsv")

rownames(annot) <- annot[, myIDcolname]
order4annot <- colnames(input_RCCs_Original_Counts_Mat)  
write.table(  annot[order4annot,] , file= path2save_annot_from_origExpMat ,  sep="\t", row.names = TRUE, col.names = TRUE, quote=FALSE)

#####
## Extracting the normalized Counts from NACHO object using the default Normalization
#####
input_RCCs_DefaultNorm_Mat <- NachoNorm2matrix( input_RCCs )

path2save_NormDefaultNacho <- paste0( outputfolder , "/","ExpMatNorm_NACHO_defaults_from_the_RCCs_in_the_folder--",basename(data_directory_path),"--.tsv")
write.table( input_RCCs_DefaultNorm_Mat , file=path2save_NormDefaultNacho ,  sep="\t", row.names = TRUE, col.names = TRUE)

# Saving the Annot for related to The Norm matrix
path2save_annot_from_NachoDefNormExpMat <- paste0( outputfolder , "/","Annot_from_ExpMatNorm_NACHO_defaults_from_the_RCCs_in_the_folder",basename(data_directory_path),"--.tsv")

rownames(annot) <- annot[, myIDcolname]
order4annot_DfNorm <- colnames(input_RCCs_DefaultNorm_Mat)  
write.table(  annot[order4annot_DfNorm, ] , file= path2save_annot_from_NachoDefNormExpMat ,  sep="\t", row.names = TRUE, col.names = TRUE, quote=FALSE)

##################### 
## Saving the Nacho obj in an RDS
###################
path2save_NachoObj <- paste0( outputfolder , "/","NACHO_Obj_from_the_RCCs_in_the_folder--",basename(data_directory_path),"--.RDS")
saveRDS( input_RCCs , file=path2save_NachoObj) 


