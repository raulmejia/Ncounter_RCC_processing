#####################
# this program 
# takes TSV derived from RCC and convert its numerical values (fourth column) to log2 
####################

######
# Required libraries
######
if (!require("argparse")) {
  install.packages("argparse", ask =FALSE)
  library("argparse")
}


############################## 
## Data given by the user
##############################
# create parser object
parser <- ArgumentParser()
# specify our desired options 
# by default ArgumentParser will add an help option 
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE ,
                    help="Print extra output [default]")
parser$add_argument("-q", "--quietly", action="store_false", 
                    dest="verbose", help="Print little output")
parser$add_argument("-t", "--tsv", type="character", 
                    help="TSV to take and transform")

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args <- parser$parse_args( )
# print some progress messages to stderr if "quietly" wasn't requested

#############################
## Reading or preparing the inputs
#############################
tsv_path <- args$tsv 
#tsv_path <- "/media/rmejia/mountme88/Projects/Maja-covid/Data/Original_RCC_files/Original_RCC_in_TSV_format/20210413_covidnma_01_12.tsv"
mytsv <-read.table( file = tsv_path , stringsAsFactors = FALSE , check.names = FALSE, header = TRUE)
 
##############################
## The program starts
#############################
result <- mytsv
result$Count <- log2(mytsv$Count +1)

write.table(result , paste0(args$tsv,".log2") , sep=",", row.names = FALSE,col.names = TRUE , quote = FALSE)
