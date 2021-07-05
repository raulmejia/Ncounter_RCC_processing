#########
## Graphical Nacho from interactive session
#########
#########
### Shamely NACHO::visualise should be run in an interactive session
### Thus I couldn't run it only from the terminal
#############
##################
### loading required packages
##################
if (!require( "BiocManager" )) {
  install.packages("BiocManager", ask =FALSE)
  library("BiocManager")
}
if (!require("NACHO")) {
  install.packages("NACHO", ask =FALSE)
  library("NACHO")
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
parser$add_argument("-d", "--datadirectory", type="character", 
                    help="path to the directory with the RCC files")
parser$add_argument("-s", "--ssheet", type="character", 
                    help="path to your ssheet")
parser$add_argument("-i", "--idcolname", type="character", 
                    help="colanme used as id")
parser$add_argument("-n", "--normalizationmethod", type="character", 
                    help="normalization method only for housekeeping?")

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args <- parser$parse_args( )
# print some progress messages to stderr if "quietly" wasn't requested

###########################
## Reading the parameters
## Substitute them by hand bc NACHO::visualise should be run in an interactive session
###########################
data_directory_path <- args$datadirectory
#data_directory_path <- "/data/.../Original_RCC_log2"
data_directory_path <- normalizePath( data_directory_path)

ssheet_path <- args$ssheet
#  ssheet_path <-"/data/.../ssheetlog2_csv.csv"

myIDcolname <- args$idcolname
# myIDcolname <- "Unique_ID"

mynormalizationmethod <- args$normalizationmethod
# mynormalizationmethod <- "GEO" # it could be GLM as well

#########
## Body of the program
#########
# Manual # vignette( "NACHO-analysis" )
####

input_RCCs <- load_rcc(data_directory = data_directory_path ,
                       ssheet_csv = ssheet_path ,
                       id_colname = myIDcolname  )

visualise( input_RCCs)
