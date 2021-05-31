# This program converts from RCC to TSVs (with gene symbols as rownames)

# Warning!! check that your are not missing rows this script only retrieves
# rows that contain these words CodeClass|Endogenous|Positive|Negative|Housekeeping

input_path=$1
#input_path=/media/rmejia/mountme88/Projects/Maja-covid/Data/Controls/Ncounter_Platform/Kidney/GSE115989/RCCs

mkdir $input_path/your_TSVs
for sample in $input_path/*.RCC;do
    name=${sample%.RCC}
    grep -E 'CodeClass|Endogenous|Positive|Negative|Housekeeping' $sample | sed 's/,/\t/g' | awk -F '\t' '{print $2,$4}'  > $input_path/your_TSVs/$( basename $name ).tsv
done

