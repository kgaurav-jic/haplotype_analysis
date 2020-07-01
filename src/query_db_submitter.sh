in_dir=$1
jelly_file=$2
out_dir=$3

linename=$(echo $(basename $jelly_file) | sed 's/.jf$//')
dir=$out_dir/$linename
mkdir -p $dir

find $in_dir | grep .db$ | while read d
do
	chunk=$(echo $(basename $d) | sed -e 's/.db$//' -e 's/database_//')
	sbatch -N 1 -c 1 -J ${linename}_$chunk --mem=5gb -p jic-long -o $dir/${linename}_$chunk.out -e $dir/${linename}_$chunk.err --wrap "python query_db.py -d $d -j $jelly_file -k 51 -o $dir/${linename}_$chunk.txt"
done