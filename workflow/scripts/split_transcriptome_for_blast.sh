INPUT=$1
OUTPUTDIR=$2
length=$(wc -l < $INPUT)
count=1
N=100000
lines=$((N+N))
i=$lines
echo "Splitting transcriptome $(basename $INPUT) into batches of $N transcripts for blast processing."  
while [ $i -lt $length ]
do
	echo "Saving transcript batch $count ($N transcripts) to $OUTPUTDIR/blast_input_$count.fasta..." 
	head -n $i $INPUT | tail -n $lines > "$OUTPUTDIR/blast_input_$count.fasta"
	count=$(($count+1))
	i=$(($i+$lines))
done
spillover=$(($length-$i+$lines))
echo "Saving transcript batch $count ($(($spillover / 2)) transcripts) to $OUTPUTDIR/blast_input_$count.fasta..." 
tail -n $spillover $INPUT > "$OUTPUTDIR/blast_input_$count.fasta"
echo "Finished."
