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

########
## Program starts
#######

############
### Building the ingredients
############

###############
### Generating the intermediate files that contain the counts 
###############
Rscript $code_path/libraries/Matrix_splited_in_OnlyCounts_and_TargetNamesPlusCounts_Flies.R \
	-m $input_matrix_path \
	-o $output_folder

##########
## Creating the generic tail for all the files
##########
mkdir -p $output_folder/tail

rm $output_folder/tail/GenericTail.txt
touch $output_folder/tail/GenericTail.txt
echo \</Code_Summary\> >> $output_folder/tail/GenericTail.txt
echo \\n >> $output_folder/tail/GenericTail.txt
echo \<Messages\> >> $output_folder/tail/GenericTail.txt
echo \</Messages\> >> $output_folder/tail/GenericTail.txt
echo \\n >> $output_folder/tail/GenericTail.txt

###########
### Generating the body and header for the new RCC files
##########
mkdir -p $output_folder/body
mkdir -p $output_folder/header
mkdir -p $output_folder/RCCs_with_Transformed_values

# Extracting the "body" from the RCC files
for sample in $input_path/*.RCC;do
    name=${sample%.RCC}
    #echo $(basename $name)
#    grep -E ',' $sample  | awk '{print $1 $2 $3}' > $output_folder/body/$( basename $name ).body123.csv
  grep -E 'CodeClass|Endogenous|Positive|Negative|Housekeeping' $sample  | awk '{print $1 $2 $3}' > $output_folder/body/$( basename $name ).body123.csv
  
    
    paste -d ',' $output_folder/body/$( basename $name ).body123.csv  $output_folder/OnlyCountColumn/$( basename $name ).RCC.OnlyCountColumn >  $output_folder/body/$( basename $name ).body123_transformed4  # here is the body

    sed '1,/<Code_Summary>/!d' $sample >  $output_folder/header/$( basename $name ).header # Generating the header 
    
    cat $output_folder/header/$( basename $name ).header \
	$output_folder/body/$( basename $name ).body123_transformed4 \
	$output_folder/tail/GenericTail.txt \
	> $output_folder/RCCs_with_Transformed_values/$( basename $name ).$label.RCC  # appending all the elements (header, body and tail)
	
  
done


#CodeClass | Endogenous | Positive | Negative | Housekeeping 













