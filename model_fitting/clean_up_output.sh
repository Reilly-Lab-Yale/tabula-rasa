temp_dir=$1
stat_suffix=$2
num_files=$3
out_file_stats=$4
out_file_resid=$5

for i in {1..$num_files} ; do
    head -n1 $i*$stat_suffix >> $out_file_stats
    
