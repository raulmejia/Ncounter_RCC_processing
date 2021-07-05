#########
## This script transforms a matrix adding it up one and then applying log2  Tranformed_matrix = log2 ( matrix + 1)
#########
### This script receives:
###     A) Expression matrix
###
##################
### loading required packages
##################
if (!require("BiocManager")) {
  install.packages("BiocManager", ask =FALSE)
  library("BiocManager")
}
if (!require("dplyr")) {
  BiocManager::install("dplyr", dependencies = TRUE)
  library("dplyr")
}
if (!require("tidyr")) {
  BiocManager::install("tidyr", dependencies = TRUE)
  library("tidyr")
}
if (!require( "psych" )) {
  BiocManager::install( "psych" , dependencies = TRUE)
  library( "psych" )
}
if (!require("argparse")) {
  install.packages("argparse", ask =FALSE)
  library("argparse")
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
parser$add_argument("-m", "--matrix", type="character", 
                    help="path to your normalization matrix")
parser$add_argument("-o", "--outputfile", type="character", 
                    help="file to save the result matrix")

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args <- parser$parse_args( )
# print some progress messages to stderr if "quietly" wasn't requested

###########################
## Reading the parameters
###########################
my_expression_matrix_path <- args$matrix
# my_expression_matrix_path <- "/data/.../subset_kidney.tsv"

outfile_path <- args$outputfile
# outfile_path <- "/data/.../subset_kidney_OAZ1_PGK1_SDHA_MRPS7.tsv"

###########################
## Reading the parameters
###########################
mymatrix <- read.table( file= my_expression_matrix_path , stringsAsFactors = FALSE , check.names = FALSE, header = TRUE)

###########################
## The calculation start
###########################
log2matrix <- log2( mymatrix + 1)

############
## save the results
############

write.table( log2matrix , file = outfile_path,  sep="\t", row.names = TRUE, col.names = TRUE, quote=FALSE)
