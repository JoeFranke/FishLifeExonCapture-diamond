#!/bin/bash
# Change the 5 to any number to change no. of iterations
for counter in `seq 2 5`; do
    for directory in *;
    do
    if [ -d $directory ];
    then
    Reference="references"
    if [[ $directory == $Reference ]]; then
	continue
    fi
    cd $directory/
        echo $counter
# This mirrors make-database.sh
	seqkit split --by-id $directory.reference.iter$(($counter - 1)).fasta
        cut -f 1 ../references/references.tsv | while read i; do
	    j=$(grep $i ../references/references.tsv | cut -f 2)
	    if [ $j -eq 2 ]; then
		cat $directory.reference.iter$(($counter - 1)).fasta.split/$directory.reference.iter$(($counter - 1)).id_$i.fasta >> $directory.references2.iter$counter.fasta
	    elif [ $j -eq 1 ]; then
		cat $directory.reference.iter$(($counter - 1)).fasta.split/$directory.reference.iter$(($counter - 1)).id_$i.fasta >> $directory.references1.iter$counter.fasta
	fi
	done

	diamond makedb --in $directory.references2.iter$counter.fasta -d reference2_iter$counter
	diamond makedb --in $directory.references1.iter$counter.fasta -d reference1_iter$counter

# This is the run-diamond.sh portion
	for file in *.trimmed.fastq; do
	    j=$( echo $(basename $file) | cut -d "." -f 1 | cut -d "_" -f 1-4)
	    diamond blastx -q $file --query-gencode 2 -d reference2_iter$counter -o "$j"_iter"$counter"_out.tsv
	    diamond blastx -q $file --query-gencode 1 -d reference1_iter$counter >> "$j"_iter"$counter"_out.tsv
	done

# Split-file portion
	for species in *R1.trimmed.fastq; do
	    j=$( echo $(basename $species) | cut -d "." -f 1 | cut -d "_" -f 1-3)
	    cut -f 2 "$j"_R1_iter"$counter"_out.tsv | sort -u | while read i; do
		cut -f 1-2 "$j"_R1_iter"$counter"_out.tsv | grep $i - | sort | uniq > "$j"_temp_1.$i.tsv
	    done
	    cut -f 2 "$j"_R2_iter"$counter"_out.tsv | sort -u | while read i; do
		cut -f 1-2 "$j"_R2_iter"$counter"_out.tsv | grep $i - | sort | uniq > "$j"_temp_2.$i.tsv
	    done
	    cat "$j"_R1_iter"$counter"_out.tsv "$j"_R2_iter"$counter"_out.tsv | cut -f 2 | sort -u | while read i; do
		if [ -f "$j"_temp_1.$i.tsv ] && [ -f "$j"_temp_2.$i.tsv ]; then
		    cat "$j"_temp_1.$i.tsv "$j"_temp_2.$i.tsv | sort > "$j"_temp_3.$i.tsv
		    uniq -d -c "$j"_temp_3.$i.tsv | sort > "$j"_temp_4.$i.tsv
		    sed 's/      2 //' "$j"_temp_4.$i.tsv > "$j"_paired_reads.$i.tsv
		    comm -13 "$j"_paired_reads.$i.tsv "$j"_temp_1.$i.tsv >> "$j"_R1_single_reads.$i.tsv
		    comm -13 "$j"_paired_reads.$i.tsv "$j"_temp_2.$i.tsv >> "$j"_R2_single_reads.$i.tsv
		    cut -f 1 "$j"_paired_reads.$i.tsv | seqkit grep -f - "$j"_R1.trimmed.fastq | sed '1~4 s/$/\/1/g' > "$j"_R1.iter$counter.$i.fq
		    cut -f 1 "$j"_paired_reads.$i.tsv | seqkit grep -f - "$j"_R2.trimmed.fastq | sed '1~4 s/$/\/2/g' > "$j"_R2.iter$counter.$i.fq
		    cut -f 1 "$j"_R1_single_reads.$i.tsv | seqkit grep -f - "$j"_R1.trimmed.fastq | sed '1~4 s/$/\/1/g' > "$j"_R1.iter$counter.$i.fq
		    cut -f 1 "$j"_R2_single_reads.$i.tsv | seqkit grep -f - "$j"_R2.trimmed.fastq | sed '1~4 s/$/\/1/g' > "$j"_R2.iter$counter.$i.fq
		elif [ -f "$j"_temp_1.$i.tsv ]; then
		    echo 'R1: "$i"'
		    cp "$j"_temp_1.$i.tsv "$j"_R1_single_reads.$i.tsv
		    cut -f 1 "$j"_R1_single_reads.$i.tsv | seqkit grep -f - "$j"_R1.trimmed.fastq | sed '1~4 s/$/\/1/g' > "$j".single.iter$counter.$i.fq
		elif [ -f "$j"_temp_2.$i.tsv ]; then
		    echo 'R2: "$j"'
		    cp "$j"_temp_2.$i.tsv "$j"_R2_single_reads.$i.tsv
		    cut -f 1 "$j"_R2_single_reads.$i.tsv | seqkit grep -f - "$j"_R2.trimmed.fastq | sed '1~4 s/$/\/1/g' > "$j".single.iter$counter.$i.fq
		else
		    echo "Error: No file found"
		fi
		done
		rm *temp*
	done

# Trinity portion
	for species in *R1.trimmed.fastq; do
	    j=$( echo $(basename $species) | cut -d "." -f 1 | cut -d "_" -f 1-3)
	    for paired in "$j"_R1.iter$counter.*.fq; do
		i=$( echo $(basename $paired) | cut -d "." -f 3-4)
		Trinity --seqType fq --left "$j"_R1.iter$counter.$i.fq --right "$j"_R2.iter$counter.$i.fq --max_memory 8G --output "$j"_trinity_paired.iter$counter."$i" --full_cleanup
	    done
	    for single in "$j".single.*.fq; do
		i=$( echo $(basename $single) | cut -d "." -f 4-5)
		Trinity --seqType fq --single "$j".single.iter$counter.$i.fq --max_memory 8G --output "$j"_trinity_single.iter$counter."$i" --full_cleanup
	    done
	done

# Exonerate, translation, and creating reference portion
	for species in *R1.trimmed.fastq; do
	    j=$( echo $(basename $species) | cut -d "." -f 1 | cut -d "_" -f 1-3)
	    if [ -f "$j".reference.iter$counter.fasta ]; then
		    rm "$j".reference.iter$counter.fasta
		fi
	    for trinities in "$j"_trinity_paired.iter$counter*.fasta; do
		i=$( echo $(basename $trinities) | cut -d "." -f 3-4)
		echo $i
		k=$( grep $i ../references/references.tsv | cut -f 2)
		exonerate --model protein2genome ../references/references.fasta.split/references.id_$i.fasta "$j"_trinity_paired.iter$counter.$i.Trinity.fasta --geneticcode $k --ryo ">%qi %qab-%qae iteration${counter} \n%tas" --showcigar F --showvulgar F --showalignment F --showsugar F --showquerygff F --showtargetgff F | sed "1,2d" | sed '$d' > "$j".iter$counter.$i.exonerate.fasta
		seqkit translate -T $k "$j".iter$counter.$i.exonerate.fasta > "$j".iter$counter.$i.exonerate_translated.fasta
		cat "$j".iter$counter.$i.exonerate_translated.fasta >> "$j".reference.iter$counter.fasta
	    done

	    for trinities in "$j"_trinity_single*.fasta; do
		i=$( echo $(basename $trinities) | cut -d "." -f 3-4)
		exonerate --model protein2genome --geneticcode $k ../references/references.fasta.split/references.id_$i.fasta "$j"_trinity_single.iter$counter.$i.Trinity.fasta --geneticcode $k --ryo ">%qi %qab-%qae iteration${counter} \n%tas" --showcigar F --showvulgar F --showalignment F --showsugar F --showquerygff F --showtargetgff F > "$j".iter$counter.$i.exonerate.fasta | sed '$d' | sed "1,2d" | seqkit translate -T $k > "$j".iter$counter.$i.exonerate.fasta
	    done
	done
cd ../
fi
done
rm references/*.iter*
done
