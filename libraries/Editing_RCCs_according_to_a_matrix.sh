#########
## This script reads the RCC files from a folder and edit their values according to a given Matrix.
##  Inputs:
##   input_path = folder hat contain the RCCs
##   input_matrix_path = path to the matrix that contain the values that you want to substitute in your RCCs (The colnames of the matrix should match exactly with the names of the RCC files)
##   code_path = path that contains the libraries/Matrix_splited_in_OnlyCounts_and_TargetNamesPlusCounts_Flies.R  script . That one that generates intermediate files which contain just the extracted counts from the matrix
##   label = some label for your resulting RCC files
##   output_folder = some folder to store your results
##
#########
##########################################
## Data given by the user
##########################################
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

#########################################
## Program starts
#########################################
###############
### Building the ingredients
##############

###############
### A) Generating the intermediate files that contain the counts 
###############
Rscript $code_path/libraries/Matrix_splited_in_OnlyCounts_and_TargetNamesPlusCounts_Flies.R \
	-m $input_matrix_path \
	-o $output_folder

##########
## B) Creating the generic tail for all the files
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
### C) Generating the body and header for the new RCC files
##########
mkdir -p $output_folder/body
mkdir -p $output_folder/header
mkdir -p $output_folder/RCCs_with_Transformed_values

# Extracting the "body" from the RCC files
for sample in $input_path/*.RCC;do
    name=${sample%.RCC}
  grep -E 'CodeClass|Endogenous|Positive|Negative|Housekeeping' $sample  | awk -F ',' '{print $1,",",$2,",",$3}' | sed 's/ //g' > $output_folder/body/$( basename $name ).body123.csv # building the body
    paste -d ',' $output_folder/body/$( basename $name ).body123.csv  $output_folder/OnlyCountColumn/$( basename $name ).RCC.OnlyCountColumn >  $output_folder/body/$( basename $name ).body123_transformed4  # here is the body

    sed '1,/<Code_Summary>/!d' $sample >  $output_folder/header/$( basename $name ).header # Generating the header 
    
    cat $output_folder/header/$( basename $name ).header \
	$output_folder/body/$( basename $name ).body123_transformed4 \
	$output_folder/tail/GenericTail.txt \
	> $output_folder/RCCs_with_Transformed_values/$( basename $name ).$label.RCC  # appending all the elements (header, body and tail)

    ######
    ## deleting the intermediate files
    #######
    rm $output_folder/body/$( basename $name ).body123.csv
    rm $output_folder/body/$( basename $name ).body123_transformed4
    rm $output_folder/header/$( basename $name ).header
    
done

##########################################
#### Deleting intermediate files
##########################################
 rmdir $output_folder/body
 rmdir $output_folder/header
 rm $output_folder/tail/GenericTail.txt
 rmdir $output_folder/tail













