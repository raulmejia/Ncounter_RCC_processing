##################
## ComBat.R
###################
### This script receive a expression matrix and a annotation file.
###     A) Expression matrix 
###     B) An annotation file:
###            (the entrances in the Unique_ID column) should exactly match with the colnames of your expression matrix
###     C) The coulmn name from the annot file to be used as batch deffinition (please indicated the batches as numbers) 
###     D) If you expect that your biological groups show markedly different variances or not. 
###           DifferentBiological_Variances_expected = TRUE or FALSE
###     E) Batch (given in a numerical way) that will be considered as refence
###         Batch_used_as_reference = 1, 2, etc... (if none = NULL) indeed that is the default 
###     
###
###
####  Example of use
##
## Further extensions:  mod = model4combat you could change it, the intention is:
##     Model matrix for outcome of interest and other covariates besides batch
##          Thus, at leat you could ask for the name of the column of your biological variable of interest
##################
### loading required packages
##################

if (!require("BiocManager")) {
  install.packages("BiocManager", ask =FALSE)
  library("BiocManager")
}
if (!require("sva")) {
  BiocManager::install("sva", dependencies = TRUE)
  library("sva")
}
if (!require("bladderbatch")) {
  BiocManager::install("bladderbatch", dependencies = TRUE)
  library("bladderbatch")
}
if (!require("pamr")) {
  BiocManager::install("pamr", dependencies = TRUE)
  library("pamr")
}
if (!require("limma")) {
  BiocManager::install("limma", dependencies = TRUE)
  library("limma")
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
parser$add_argument("-m", "--inputmatrix", type="character", 
                    help="The Matix given as input")
parser$add_argument("-a", "--annotation", type="character", 
                    help="Annotation File")
parser$add_argument("-b", "--batchcolumn", type="character", 
                    help="column from your annotation file that contains your batches")
parser$add_argument("-d", "--differentvariances", type="character", 
                    help="Are Diffent variances in your biological groups expected? (TRUE / FALSE)")
parser$add_argument("-r", "--reference", type="character", 
                    help="Which bacht (in numeric format) will be used as reference, (or NULL)")
parser$add_argument("-o", "--output", type="character", 
                    help="Path for your output matrix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args <- parser$parse_args( )
# print some progress messages to stderr if "quietly" wasn't requested

###########################
## Reading the parameters
###########################
path_to_input_matix <- args$inputmatrix
#  path_to_input_matix <- "/media/rmejia/mountme88/Projects/Maja-covid/Results/Normalizations/NK_Geo/MajaPreprocessed_GlomTubPreprocessed_OAZ1_HPRT1_ABCF1.tsv"

path_to_annotation <- args$annotation
# path_to_annotation <- "/media/rmejia/mountme88/Projects/Maja-covid/Results/Annot_MajaWhole_GSE113342/Annot_MajaWholw_GSE113342_NumericBatch.tsv"

Colname_with_batches <- args$batchcolumn
# Colname_with_batches <- "NumBatch"

DifferentBiological_Variances_expected <- args$differentvariances
# DifferentBiological_Variances_expected <- "TRUE"
DifferentBiological_Variances_expected <- as.logical(DifferentBiological_Variances_expected)

Batch_used_as_reference <- args$reference
# Batch_used_as_reference <- "NULL"
if(Batch_used_as_reference == "NULL") {
    Batch_used_as_reference <- NULL
} 

path_output  <- args$output
# path_output  <-"/media/rmejia/mountme88/Projects/Maja-covid/Results/WholeMaja_GSE113342/HK_then_Combat/MajaPreprocessed_GlomTubPreprocessed_OAZ1_HPRT1_ABCF1_ComBat.tsv"

#########
## Body program
#########
inputmatrix <-read.table(file=path_to_input_matix
                       , stringsAsFactors = FALSE, check.names = FALSE)

annot <- read.table(file = path_to_annotation, stringsAsFactors = FALSE, check.names = FALSE, header = TRUE,  sep = "\t")

########
### Combat!
#######
batches = as.factor(annot[ , Colname_with_batches ])
model4combat = model.matrix ( ~1 , data = annot)
# model4combat_biol = model.matrix ( ~ Tissue , data = annot)

matrix_combated <- ComBat( dat = inputmatrix , batch=batches , mod = model4combat , par.prior = TRUE, prior.plots = FALSE,
ref.batch = Batch_used_as_reference , mean.only = DifferentBiological_Variances_expected)

dir.create( dirname(path_output) , recursive = TRUE)
write.table( matrix_combated , file= path_output ,  sep="\t", row.names = TRUE, col.names = TRUE)
