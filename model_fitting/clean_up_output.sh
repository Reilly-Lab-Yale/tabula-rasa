temp_dir=$1
stat_suffix=$2
num_files=$3
out_file_stats=$4

cd $temp_dir

echo -e "model_choice\tformula\tmaxiter\tregfit\tconverged\taic\tbic\tloglike\tllr_chi2\tllr_chi2_pval\tpsuedor2\tresids" > $out_file_stats
for i in $(seq 1 $num_files) ; do
    head -n1 $i*$stat_suffix >> $out_file_stats
done