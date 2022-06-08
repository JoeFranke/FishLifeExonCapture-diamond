for directory in *;
do
if [ -d $directory ];
then
cd $directory/
    for file in *.trimmed.fastq; do
        j=$( echo $(basename $file) | cut -d "." -f 1 | cut -d "_" -f 1-4)
        diamond blastx -q $file --query-gencode 2 -d ../references/reference2 -o "$j"_out.tsv
        diamond blastx -q $file --query-gencode 1 -d ../references/reference1 >> "$j"_out.tsv
    done
fi
cd ../
done
