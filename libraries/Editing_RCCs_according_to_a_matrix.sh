#########
## This script reads the RCC files from a folder and edit their values according to a Matrix.
##  Inputs
##   input_path = folder hat contain the RCCs
##   input_matrix_path = path to the matrix that contain the values that you want to substitute in your RCCs
##   code_path = path that contains the libraries/Matrix_splited_in_OnlyCounts_and_TargetNamesPlusCounts_Flies.R  script . That one that generates intermediate files (that contain just the counts) from the matrix
##   label = some label for your resulting RCC files
##   output_folder = some folder to store your results
##
#########


input_path=$1
input_matrix_path=$2
code_path=$3
label=$4
output_folder=$5

echo $input_path
echo $input_matrix_path
echo $code_path
echo $label
echo $output_folder

### Generating the intermediate files 
Rscript $code_path/libraries/Matrix_splited_in_OnlyCounts_and_TargetNamesPlusCounts_Flies.R \
	-m $input_matrix_path \
	-o $output_folder

#echo $output_folder/OnlyCountColumn/*


mkdir -p $output_folder/body
mkdir -p $output_folder/header
# Extracting the "body" from the RCC files
for sample in $input_path/*.RCC;do
    name=${sample%.RCC}
    #echo $(basename $name)
    grep -E '>' $sample  | awk '{print $1 $2 $3}' > $output_folder/body/$( basename $name ).body123.csv
    paste $output_folder/body/$( basename $name ).body123.csv  $output_folder/OnlyCountColumn/$( basename $name ).RCC.OnlyCountColumn >  $output_folder/body/$( basename $name ).body123_transformed4

    sed '1,/<Code_Summary>/!d' $sample >  $output_folder/header/$( basename $name ).header
done















