# symlink a bunch of files from a list of paths
# eg
for i in `cat $1`; do 
	x=$i/*Aligned.Quality.Sorted.bam.junc
	base=$(basename $x)
	ln -s $x $base
done
