# This script gives a boxplot visualization of selected genes according conditions defined in a annotation file
# The structure of your matrix should be 
#          sample1 sample2 ...
# rowname1  0.42    45
# rowname2  1       0
# rowname3  45      0
# NEG_Prob1 0       12
# POS_E     2       1
## Notes
##  The annotation file should have a column called "Unique_ID"

############################## 
## Required libraries
##############################
if (!require("BiocManager")) {
  install.packages("BiocManager", ask =FALSE)
  library("BiocManager")
}
if (!require("ggplot2")) {
  BiocManager::install("ggplot2", ask =FALSE)
  library("ggplot2")
}
if (!require("argparse")) {
  install.packages("argparse", ask =FALSE)
  library("argparse")
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", ask =FALSE)
  library("RColorBrewer")
}
if (!require("reshape2")) {
  install.packages("reshape2", ask =FALSE)
  library("reshape2")
}
if (!require("affy")) {
  BiocManager::install("affy", ask =FALSE)
  library("affy")
}

############################## 
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
parser$add_argument("-m", "--matrix", type="character", 
                    help="path to your expression matrix")
parser$add_argument("-a", "--annotation", type="character", 
                    help="path to your annotation file")
parser$add_argument("-c", "--code", type="character", 
                    help="path to your code")
parser$add_argument("-l", "--label", type="character", 
                    help="label to your results")
parser$add_argument("-g", "--maingroups", type="character", 
                    help="the name of your column to correct / make intrabatch normalization")
parser$add_argument("-s", "--selectedgenesfromrownamesofadataframe", type="character", 
                    help="A data frame that contains the genes that you want to plot in the row names, only the rownames will be used for that tash")
parser$add_argument("-o", "--outputfolder", type="character", 
                    help="output folder where you want to store your results")

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args <- parser$parse_args( )
# print some progress messages to stderr if "quietly" wasn't requested

#############################
## Reading or preparing the inputs
#############################
mymatrix <-read.table( file=args$matrix, stringsAsFactors = FALSE , check.names = FALSE)
#  mymatrix <-read.table(file="/media/rmejia/mountme88/Projects/Maja-covid/Results/Normalizations/NK_Geo/Majalog2_OAZ1_HPRT1_ABCF1.tsv", stringsAsFactors = FALSE, check.names = FALSE)

annotdf <-read.table( file=args$annotation, stringsAsFactors = FALSE, header = TRUE )
# annotdf <-read.table(file="/media/rmejia/mountme88/Projects/Maja-covid/Results/Preprocessing_through_Log2/NachoNorm/Annot_from_ExpMat_as_input_from_the_RCCs_in_the_folder--Original_RCC_log2--.tsv", stringsAsFactors = FALSE, , header = TRUE)

code_path <- args$code
# code_path <- "/media/rmejia/mountme88/code/Ncounter_RCC_processing/"

label <- args$label # label <- "First_run"
your_main_groups <- args$maingroups # your_main_groups <- "Tissue"

DataFrame_with_the_selected_genes_in_the_rownames <-read.table( file=args$selectedgenesfromrownamesofadataframe, stringsAsFactors = FALSE, header = FALSE )
#


outputfolder <- args$outputfolder
#  outputfolder <- "/media/rmejia/mountme88/Projects/Maja-covid/Results/Normalizations/NK_Geo/Majalog2_OAZ1_HPRT1_ABCF1/BoxPlot"

dir.create(outputfolder, recursive = TRUE)
outputfolder <- normalizePath(outputfolder)

code_path <- normalizePath(code_path)

##############################
## The program starts
#############################
if(all(colnames( mymatrix) ==  annotdf$Unique_ID) == TRUE ){ 
  print("your annotation and colnames match")
}
if(all(colnames( mymatrix) ==  annotdf$Unique_ID) != TRUE ){ 
  print("ERROR: your annotation and colnames DON´T match")
  break()
}

# Extracting a submatrix that contains only the genes of interest
Submatrix_with_only_the_genes_of_interest <- mymatrix[ DataFrame_with_the_selected_genes_in_the_rownames[,1] , ]

## Melting the data
source( paste0( code_path ,"/libraries/" , "matrix_N_annotdf_2_melteddf.R") )
melteddf <- matrix_N_annotdf_2_melteddf( Submatrix_with_only_the_genes_of_interest , annotdf )

############
## graphs Boxplots and saving
############
pdf( file=paste0(outputfolder,"/" , label ,".pdf"),
     width = 10, height = 7)
ggplot(melteddf, aes_string(y="value", fill = your_main_groups, x = "variable")) + 
  geom_boxplot( )+ ggtitle(paste( label ))  # Density plots
ggplot(melteddf, aes_string(y="value", fill = "variable", x = your_main_groups)) + 
  geom_boxplot( )+ ggtitle(paste( label ))  # Density plots
dev.off()

