for directory in *;
do
if [ -d $directory ];
then
cd $directory
    for species in *R1.trimmed.fastq; do
       	j=$( echo $(basename $species) | cut -d "." -f 1 | cut -d "_" -f 1-3)
	if [ -f "$j".reference.iter1.fasta ]; then
		rm "$j".reference.iter1.fasta
	    fi
        for trinities in "$j"_trinity_paired*.fasta; do
	    # Finds the reference sequence name
	    i=$( echo $(basename $trinities) | cut -d "." -f 2-3)
	    # If reference file already exists, remove it
	    # finds the genetic code needed for exonerate and translation
	    k=$( grep $i ../references/references.tsv | cut -f 2)
	    # Runs exonerate
            exonerate --model protein2genome ../references/references.fasta.split/references.id_$i.fasta "$j"_trinity_paired.$i.Trinity.fasta --geneticcode $k --ryo ">%qi %qab-%qae iteration1 \n%tas" --showcigar F --showvulgar F --showalignment F --showsugar F --showquerygff F --showtargetgff F | sed '$d' | sed "1,2d" > "$j".$i.exonerate.fasta
	    # Runs translation
	    seqkit translate -T $k "$j".$i.exonerate.fasta > "$j".$i.exonerate_translated.fasta
	    #Creates new reference
	    cat "$j".$i.exonerate_translated.fasta >> "$j".reference.iter1.fasta
        done
    done

cd ../
fi
done
