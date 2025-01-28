#!/bin/bash
#SBATCH -A Research_Project-MRC164847 # research project to submit under.
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq
#SBATCH --time=15:00:00 # Maximum wall time for the job
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors.
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=m.kouhsar@exeter.ac.uk # email address

#########################################################################################
#########################################################################################

out_dir=./circRNA
fastq_dir=./circRNA
genome_fasta=./Ref_Genome/GRCh38.p14.genome.fa
genome_gtf=./Ref_Genome/gencode.v47.chr_patch_hapl_scaff.annotation.gtf
thread=16
#######################################################################################
#######################################################################################

out_dir_call=${out_dir}/CIRI.long.Call
out_dir_collapse=${out_dir}/CIRI.long.Collapse
fastq_files=(${fastq_dir}/*.fastq)

echo Output directory: $out_dir
echo Fastq files directory: $fastq_dir
echo Genome fasta file: $genome_fasta
echo Genome fasta file: $genome_gtf
echo Number of samples: ${#fastq_files[@]}
echo Number of threads: $thread

echo "##########################################################################"
echo -e '\n'

mkdir -p $out_dir_call
mkdir -p $out_dir_collapse

for f in ${fastq_files[@]}
do
	f_name=$(basename $f)
	f_name=${f_name%".fastq"}
	echo "*********************************************************************************"
	echo "*********************************************************************************"
	echo "                Running CIRI-long Call on ${f_name}:                             "
	echo "*********************************************************************************"
	echo "*********************************************************************************"
	CIRI-long call -i $f \
		-o ${out_dir_call}/${f_name} \
		-r $genome_fasta \
		-p $f_name \
		-a $genome_gtf \
		-t $thread
	echo $f_name ${out_dir_call}/${f_name}/${f_name}.cand_circ.fa > ${out_dir_call}/${f_name}/${f_name}.lst
	echo "*********************************************************************************"
	echo "*********************************************************************************"
	echo "                Running CIRI-long Collapse on ${f_name}:                         "
	echo "*********************************************************************************"
	echo "*********************************************************************************"
	CIRI-long collapse -i ${out_dir_call}/${f_name}/${f_name}.lst \
		-o ${out_dir_collapse}/${f_name} \
		-p $f_name \
		-r $genome_fasta \
		-a $genome_gtf \
		-t $thread
done

echo "All the process is done!"


