#!/bin/bash

samples=($(ls raw_reads))

for i in "${samples[@]}"
do
	TEST1=$(ls raw_reads/"$i"/*R1* | tr '\n' ' ' | \
	awk '{for(i=1;i<=NF;i++){printf $i;if(i<NF)printf ","}}')

	TEST2=$(ls raw_reads/"$i"/*R2* | tr '\n' ' ' | \
	awk '{for(i=1;i<=NF;i++){printf $i;if(i<NF)printf ","}}')

	run_rnacocktail.py quantify --quantifier_idx human_transcriptome_19_index/ --1 $TEST1 --2 $TEST2 \
	--libtype ISR --salmon_k 19 --outdir salmon_out --workdir salmon_work --sample $i --unzip

	echo '\n'
	echo '\n'
	echo "Done with $i"
	echo '\n'
	echo '\n'
done
