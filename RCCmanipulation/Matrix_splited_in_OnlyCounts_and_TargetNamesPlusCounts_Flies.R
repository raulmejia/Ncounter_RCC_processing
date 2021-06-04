#########
## Matrix to Column counts & Name plus Column counts
#########
### This script receive a matrix and a output folder. It will retreive:
###     A) A set of files only with the numeric column from the matrix (the colname of that column will be "Count") keeping the name of the matrix column in the name of the file "Columnname.count"
###     B) Similar files but with a second column, column name = Name and the content the Matrix rownames 
###
############
## 
############
if (!require("ggplot2")) {
  BiocManager::install("ggplot2", ask =FALSE)
  library("ggplot2")
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
parser$add_argument("-m", "--inputmatrix", type="character", 
                    help="path to the input matrix")
parser$add_argument("-o", "--outputfolder", type="character", 
                    help="Folder to store the outputs")

########################
## Reading the data
########################
inputmatrix_path <- args$inputmatrix
# inputmatrix_path <- "/media/rmejia/mountme88/Projects/Maja-covid/Results/Preprocessing_through_Log2/NachoNorm/ExpMat_as_input_from_the_RCCs_in_the_folder--Original_RCC_log2--.tsv"

outputfolder <- args$outputfolder
# outputfolder <- "/media/rmejia/mountme88/Projects/Maja-covid/Results/Preprocessing_through_Log2/Just_ComBat_Equal_variances_assumed_No_Quantiles/"
dir.create(outputfolder, recursive = TRUE) ; outputfolder <- normalizePath( outputfolder )

inputmatrix <- read.table( file = inputmatrix_path , sep="\t", header=TRUE, check.names = FALSE)


dir.create( paste0( outputfolder,"/OnlyCountColumn/") , recursive = TRUE)
dir.create( paste0( outputfolder,"/Names_n_CountColumn/") , recursive = TRUE)

            
for( k in 1:dim(inputmatrix)[2]   ){ 
  Result_only_count_column <- data.frame( Count = inputmatrix[ , colnames(inputmatrix)[k] ]   )
  write.table(  Result_only_count_column, 
              file = paste0( outputfolder,"/OnlyCountColumn/",colnames(inputmatrix)[k],".OnlyCountColumn"),
              row.names = FALSE, col.names = TRUE, quote = FALSE) 
  
  Result_NAmes_n_Count_column <- data.frame( Names= rownames( inputmatrix), Count = inputmatrix[ , colnames(inputmatrix)[k] ]  )
  write.table(  Result_NAmes_n_Count_column , 
                file = paste0( outputfolder,"/Names_n_CountColumn/",colnames(inputmatrix)[k],".Names_n_CountColumn"),
                row.names = FALSE, col.names = TRUE, quote = FALSE) 
}





