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
###         Number_of_Batch_In_yout_batch_description_that_you_want_to_use_as_reference = 1, 2, etc... (if none = NULL) indeed that is the default 
###     
###
####  Example of use

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
                    help="Are Diffent variances in your biological groups expected?")
parser$add_argument("-r", "--reference", type="character", 
                    help="Which bacht (in numeric format) will be used as reference")
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

<- args$annotation
  <- args$batchcolumn
  <- args$differentvariances
  <- args$reference
  <- args$output
  <- args$
  <- args$
  <- args$


