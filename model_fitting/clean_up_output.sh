temp_dir=$1
stat_suffix=$2
num_files=$3
out_file_stats=$4

cd $temp_dir

echo -e "model_choice\tformula\tmaxiter\tzi_param\tconverged\taic\tbic\tloglike\tllr_chi2\tllr_chi2_pval\tpsuedor2\tresids" > $out_file_stats
echo -e "model_choice\tformula\tmaxiter\tzi_param\tconverged\taic\tbic\tloglike\tllr_chi2\tllr_chi2_pval\tpsuedor2\tresids" > ${out_file_stats}_fix
for i in $(seq 1 $num_files) ; do
    head -n1 $i*$stat_suffix >> $out_file_stats
    sed -n 1  $i*$stat_suffix >> ${out_file_stats}_fix
done