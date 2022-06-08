for directory in *;
do
if [ -d $directory ];
then
cd $directory
    for species in *R1.trimmed.fastq; do

        j=$( echo $(basename $species) | cut -d "." -f 1 | cut -d "_" -f 1-3)
        for paired in "$j"_R1.*.fq; do
            i=$( echo $(basename $paired) | cut -d "." -f 2-3) 
            Trinity --seqType fq --left "$j"_R1.$i.fq --right "$j"_R2.$i.fq --max_memory 8G --output "$j"_trinity_paired."$i" --full_cleanup
        done
        for single in "$j".single.*.fq; do
	    i=$( echo $(basename $single) | cut -d "." -f 3-4)
	    Trinity --seqType fq --single "$j".single.$i.fq --max_memory 8G --output "$j"_trinity_single."$i" --full_cleanup
        done
    done
cd ../
fi
done
