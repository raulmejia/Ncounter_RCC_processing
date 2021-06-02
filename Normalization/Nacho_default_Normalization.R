#########
## Nacho_default_Normalization.R
#########
### This script receive a folder with RCC files and retreives:
###     A) Expression matrix of the InputCounts (as expressed in the NACHO object)
###     B) Expression matrix with the default normalization given by NACHO
###     C) Nacho Object
###
###     Document that describes Nacho normalization
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
parser$add_argument("-n", "--normalizationmethod", type="character", 
                    help="normalization method only for housekeeping?")
parser$add_argument("-c", "--codepath", type="character", 
                    help="Path where you have your code")
parser$add_argument("-o", "--outputfolder", type="character", 
                    help="Folder to store the outputs")
parser$add_argument("-g", "--maingroups", type="character", 
                    help="The most reelevant column in your annotation file")
parser$add_argument("-l", "--label", type="character", 
                    help="label for your results")


# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args <- parser$parse_args( )
# print some progress messages to stderr if "quietly" wasn't requested

###########################
## Reading the parameters
## Substitute them by hand bc NACHO::visualise should be run in an interactive session
###########################
data_directory_path <- args$datadirectory
#data_directory_path <- "/media/rmejia/mountme88/Projects/Maja-covid/Data/Original_RCC_log2"
data_directory_path <- normalizePath( data_directory_path)

ssheet_path <- args$ssheet
#  ssheet_path <-"/media/rmejia/mountme88/Projects/Maja-covid/Data/ssheetlog2_csv.csv"

myIDcolname <- args$idcolname
# myIDcolname <- "Unique_ID"

mynormalizationmethod <- args$normalizationmethod
# mynormalizationmethod <- "GEO" # it could be GLM as well

code_path <- args$codepath
# code_path <- "/media/rmejia/mountme88/code/Ncounter_RCC_processing/"
code_path <- normalizePath( code_path )

outputfolder <- args$outputfolder
# outputfolder <- "/media/rmejia/mountme88/Projects/Maja-covid/Results/Preprocessing_through_Log2/NachoNorm"

your_main_groups  <- "Tissue"


label <- "Log2_NACHO_4HK_GLM"

#########
## Body of the program
#########
# Manual # vignette( "NACHO-analysis" )
####

input_RCCs <- load_rcc(data_directory = data_directory_path ,
                       ssheet_csv = ssheet_path ,
                       id_colname = myIDcolname  )