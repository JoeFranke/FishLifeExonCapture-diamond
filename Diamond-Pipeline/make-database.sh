cd references/

# if reference files exist already remove them
if [ -f "references2.fasta" ]; then
    rm references2.fasta
fi

if [ -f "references1.fasta" ]; then
    rm references1.fasta
fi

# This creates two seperate reference files based on which genetic code they use
cut -f 1 references.tsv | while read i; do
	j=$(grep $i references.tsv | cut -f 2)
	echo $j
	if [ $j -eq 2 ]; then
	    cat references.fasta.split/references.id_$i.fasta >> references2.fasta
	elif [ $j -eq 1 ]; then
	    cat references.fasta.split/references.id_$i.fasta >> references1.fasta
	fi
done
# makes database out of each reference file
diamond makedb --in references2.fasta -d reference2
diamond makedb --in references1.fasta -d reference1

cd ../
