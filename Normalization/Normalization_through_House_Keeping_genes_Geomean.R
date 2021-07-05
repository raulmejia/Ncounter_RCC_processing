#########
## House keeping Norm_Geomean
#########
### This script receives:
###     A) Expression matrix
###     B) A file with the House Keeping (HK) genes to be used.
###
###
#########
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
parser$add_argument("-k", "--housekeepingfile", type="character", 
                    help="File with the House Keeping Genes")
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

my_house_keeping_genes_file_path <- args$housekeepingfile
# my_house_keeping_genes_file_path <-"/data/.../OAZ1_PGK1_SDHA_MRPS7.tsv"

outfile_path <- args$outputfile
# outfile_path <- "/data/.../subset_kidney_OAZ1_PGK1_SDHA_MRPS7.tsv"

###########################
## Reading the parameters
###########################
mymatrix <- read.table( file= my_expression_matrix_path , stringsAsFactors = FALSE , check.names = FALSE)

house_keeping_genes <- read.table( file = my_house_keeping_genes_file_path , stringsAsFactors = FALSE , header=FALSE) ; hk_vec<- house_keeping_genes[,1]  # Getting the HK genes as a character vector

###########################
## The calculation start
###########################
# Getting the submatrix of hk genes
hkgenes_from_matrix <- intersect( hk_vec , rownames( mymatrix )  ) 
Submatix_of_hk_genes_4_normalization <- mymatrix[ hk_vec , ]

# getting the normalization factors
HK_GeoMean_per_sample <- apply( Submatix_of_hk_genes_4_normalization , 2 , geometric.mean, na.rm = TRUE   ) 
#Aritmetic_mean_of_GeoMeans <- mean( HK_GeoMean_per_sample )
HK_Normalization_factors <-  mean( HK_GeoMean_per_sample ) /  HK_GeoMean_per_sample

hk_geom_normalized_matrix <- t( t( mymatrix )  * HK_Normalization_factors)

############
## save the results
############

write.table( hk_geom_normalized_matrix , file = outfile_path,  sep="\t", row.names = TRUE, col.names = TRUE, quote=FALSE)
