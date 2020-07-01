in_dir=$1
out_dir=$2

linename=$(echo $(basename $in_dir) | sed 's/_chunks$//')
mkdir -p $out_dir

find $in_dir | grep .fa$ | while read f
do
	o=$(echo $(basename $f) | sed 's/.fa$//')
	sbatch -N 1 -c 1 -J ${o} --mem=5gb -p jic-long -o $out_dir/$o.out -e $out_dir/$o.err --wrap "python seq_parser.py -a $f -l $linename -k 51 -c jellies.cfg -o $out_dir/$o.db"
done