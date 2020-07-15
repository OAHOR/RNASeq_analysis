#!/bin/bash
#Prefetch and especially fastq-dump can be fairly slow due to their single-threaded nature.
#Here, I use a simple bash multi-threading method with an n=4 threads. This will spawn circa
#4-6 workers. Setting this number higher while using a typical SATA HDD will not see any
#speed advantage.
N=4
echo "This bash script should be run in a (conda) environment which contains"
echo "the SRA-toolkit. You can install this via conda install sra-tools"

#First go to the SRA run selector (https://www.ncbi.nlm.nih.gov/Traces/study/) and download
#the metadata file (StaRunTable.txt) from a set of samples.
#Extract the SRR identifiers:
VAR=$(tail -n +2 SraRunTable.txt | cut -d ',' -f 1)

# This is a loop for downloading and processing the data
for i in ${VAR}
do
	(
	#First print prefetch version. This cannot be combined with other options
	prefetch --version
	#Use prefetch to download. Resume interrupted download and verify download
	prefetch --resume yes --verify yes ${i}

	#Now use the obtained SRA file to write fastq file. Check if file exists
	if [ -f ${i}.sra ]
	    then
		echo "${i} already downloaded"
	else
		echo "(o) Converting SRA entry: ${i}"
    #For personal preference, I set the quality identifier line simply to '+' as it is
    #of no use to me while taking up space. 
		fastq-dump --gzip --defline-qual '+' ${i}
		echo "(o) Done converting ${i}"
	fi
	) &

	# allow only to execute $N jobs in parallel
	if [[ $(jobs -r -p | wc -l) -gt $N ]]; then
	# wait only for first job
	wait -n
    fi
done

# wait for pending jobs
wait
echo "all done" 
