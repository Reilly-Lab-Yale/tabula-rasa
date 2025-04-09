temp_dir=$1
stat_suffix=stats.txt
num_files=$3
out_file_stats=$4

cd $temp_dir

echo -e "idx\tmodel_choice\tformula\tmaxiter\tzi_param\tconverged\taic\tbic\tloglike\tllr_chi2\tllr_chi2_pval\tpsuedor2" > $out_file_stats
echo -e "idx\tmodel_choice\tformula\tmaxiter\tzi_param\tconverged\taic\tbic\tloglike\tllr_chi2\tllr_chi2_pval\tpsuedor2" > ${out_file_stats}_fix
for i in $(seq 1 $num_files) ; do
    #head -n1 $i_*$stat_suffix >> $out_file_stats
    sed -n "1p"  ${i}_*$stat_suffix >> ${out_file_stats}_fix
done