# 
# This scripts creates an expression matrix in a tsv format
# Inputs a Folder that contains RCC files (there shouldn't be any space in the RCC file names)
# 
# Example of use :
# sh /path/RCC_2_tsv-expression_matrix.sh /path/with/your/.RCC/files /path/.RCC_file/to/take_the_rownames.RCC /path/outcome/result/path/GSE115989_matrix.tsv

input_path=$1
file_extract_rownames=$2
output_path=$3
#file_extract_rownames=/media/rmejia/mountme88/Projects/Maja-covid/Data/Controls/Ncounter_Platform/Kidney/GSE115989/RCCs/GSM3204351_20160729_072916_01_122_ABMR.RCC

mkdir $input_path/your_matrix

# Getting only the numeric part of the matrix
for sample in $input_path/*.RCC;do
    name=${sample%.RCC}
    grep -E 'CodeClass|Endogenous|Positive|Negative|Housekeeping' $sample | sed 's/,/\t/g' | awk -F '\t' '{print $4}' | grep -E '[0-9]'  > $input_path/your_matrix/$( basename $name ).just_numbers.tsv
done

paste -d '\t' $input_path/your_matrix/*.just_numbers.tsv | sed -e "s/\r//g" > $input_path/your_matrix/matrix.tsv

### Rownames
grep -E 'CodeClass|Endogenous|Positive|Negative|Housekeeping' $file_extract_rownames | sed 's/,/\t/g' | awk -F '\t' '{print $2}' | sed "1 d" >  $input_path/your_matrix/rownames.rownames

# paste numeric content and Rownames
paste -d '\t' $input_path/your_matrix/rownames.rownames $input_path/your_matrix/matrix.tsv | sed -e "s/\r//g" > $input_path/your_matrix/matrix_with_rownames.tsv 

# The colnames
for sample in $input_path/*.RCC;do
    basename $sample >> $input_path/your_matrix/columnames.names
done

sed 's/\.RCC//g' $input_path/your_matrix/columnames.names > $input_path/your_matrix/columnames_no_extension.names
tr "\n" "\t" < $input_path/your_matrix/columnames_no_extension.names > $input_path/your_matrix/line_colnames

#Paste colnames and matrix
awk 'NF' $input_path/your_matrix/line_colnames $input_path/your_matrix/matrix_with_rownames.tsv  > $output_path

# Creaning the work space
rm $input_path/your_matrix/*.just_numbers.tsv
rm $input_path/your_matrix/columnames.names
rm $input_path/your_matrix/rownames.rownames
rm $input_path/your_matrix/matrix.tsv
rm $input_path/your_matrix/line_colnames
rm $input_path/your_matrix/matrix_with_rownames.tsv
rm $input_path/your_matrix/columnames_no_extension.names
rmdir $input_path/your_matrix/
