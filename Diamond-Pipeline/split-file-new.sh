#!/bin/bash
for directory in *;
do
if [ -d $directory ];
then
cd $directory/
for species in *R1.trimmed.fastq; do
# Cuts out the Family_Genus_species_S1234 part of the name
    j=$( echo $(basename $species) | cut -d "." -f 1 | cut -d "_" -f 1-3)
# Loops through every reference sequence present in each tsv
    cut -f 2 "$j"_R1_out.tsv | sort -u | while read i; do
	# cuts the first two fields out, greps to find all lines that share a reference sequence and puts them in temp_1
        cut -f 1-2 "$j"_R1_out.tsv | grep $i - | sort | uniq > "$j"_temp_1.$i.tsv
    done
 # Repeats previous step with R2
    cut -f 2 "$j"_R2_out.tsv | sort -u | while read i; do
        cut -f 1-2 "$j"_R2_out.tsv | grep $i - | sort | uniq > "$j"_temp_2.$i.tsv
    done
# Adds both tsvs together and loop through the reference sequences
    cat "$j"_R1_out.tsv "$j"_R2_out.tsv | cut -f 2 | sort -u | while read i; do
	# If both R1 and R2 have corresponding files for a reference sequence, both hit that reference sequence
	if [ -f "$j"_temp_1.$i.tsv ] && [ -f "$j"_temp_2.$i.tsv ]; then
	    # Make temp_3
            cat "$j"_temp_1.$i.tsv "$j"_temp_2.$i.tsv | sort > "$j"_temp_3.$i.tsv
	    # Find all values that appeared twice
            uniq -d -c "$j"_temp_3.$i.tsv | sort > "$j"_temp_4.$i.tsv
	    # Get rid of the number that shows how many times they appeared and put them into paired_reads file
            sed 's/      2 //' "$j"_temp_4.$i.tsv > "$j"_paired_reads.$i.tsv
	    # Compare paired_reads file to temp 1 to see all values unique to temp1 meaning only appeared in R1
            comm -13 "$j"_paired_reads.$i.tsv "$j"_temp_1.$i.tsv >> "$j"_R1_single_reads.$i.tsv
	    # Find all the reads that only appeared in R2
            comm -13 "$j"_paired_reads.$i.tsv "$j"_temp_2.$i.tsv >> "$j"_R2_single_reads.$i.tsv
	    # Put all the paired reads into R1 and R2 respectively
            cut -f 1 "$j"_paired_reads.$i.tsv | seqkit grep -f - "$j"_R1.trimmed.fastq | sed '1~4 s/$/\/1/g' > "$j"_R1.$i.fq
            cut -f 1 "$j"_paired_reads.$i.tsv | seqkit grep -f - "$j"_R2.trimmed.fastq | sed '1~4 s/$/\/2/g' > "$j"_R2.$i.fq
	    # Put the single end reads into R1 and R2 respectively
            cut -f 1 "$j"_R1_single_reads.$i.tsv | seqkit grep -f - "$j"_R1.trimmed.fastq | sed '1~4 s/$/\/1/g' >> "$j"_R1.$i.fq
            cut -f 1 "$j"_R2_single_reads.$i.tsv | seqkit grep -f - "$j"_R2.trimmed.fastq | sed '1~4 s/$/\/1/g' >> "$j"_R2.$i.fq
	# This handles if one end aligned onto a reference sequence that the other end did not
	elif [ -f "$j"_temp_1.$i.tsv ]; then
	    echo 'R1: "$i"'
	    cp "$j"_temp_1.$i.tsv "$j"_R1_single_reads.$i.tsv
	    cut -f 1 "$j"_R1_single_reads.$i.tsv | seqkit grep -f - "$j"_R1.trimmed.fastq | sed '1~4 s/$/\/1/g' > "$j".single.$i.fq
	elif [ -f "$j"_temp_2.$i.tsv ]; then
	    echo 'R2: "$j"'
	    cp "$j"_temp_2.$i.tsv "$j"_R2_single_reads.$i.tsv
	    cut -f 1 "$j"_R2_single_reads.$i.tsv | seqkit grep -f - "$j"_R2.trimmed.fastq | sed '1~4 s/$/\/1/g' > "$j".single.$i.fq
	else
	    echo "Error: No file found"
	fi
    done
done
cd ../
fi
done
