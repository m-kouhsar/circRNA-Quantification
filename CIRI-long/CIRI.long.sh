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

out_dir./Results
fastq_dir=./Raw
genome_fasta=./GRCh38.p14.genome.fa
genome_gtf=./gencode.v47.chr_patch_hapl_scaff.annotation.gtf
circ_annot=
thread=16
collapse_only=yes
#######################################################################################
#######################################################################################

out_dir_call=${out_dir}/CIRI.long.Call
out_dir_collapse=${out_dir}/CIRI.long.Collapse
fastq_files=(${fastq_dir}/*.fastq)
collapse_only=$(echo $collapse_only | xargs)
collapse_only=$(echo $collapse_only | tr '[:upper:]' '[:lower:]')
circ_annot=$(echo $circ_annot | xargs)

echo Output directory: $out_dir
echo Fastq files directory: $fastq_dir
echo Genome fasta file: $genome_fasta
echo Genome fasta file: $genome_gtf
echo circRNA anootation file (optional): $circ_annot
echo Running collapse mode only? $collapse_only
echo Number of CPU cores: $thread

echo "##########################################################################"
echo -e '\n'

mkdir -p $out_dir_call
mkdir -p $out_dir_collapse

for f in ${fastq_files[@]}
do
	f_name=$(basename $f)
	f_name=${f_name%".fastq"}

	if [ "$collapse_only" = "yes" ]
	then
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
	fi

	echo $f_name ${out_dir_call}/${f_name}/${f_name}.cand_circ.fa > ${out_dir_call}/${f_name}/${f_name}.lst
	echo "*********************************************************************************"
	echo "*********************************************************************************"
	echo "                Running CIRI-long Collapse on ${f_name}:                         "
	echo "*********************************************************************************"
	echo "*********************************************************************************"
	if [ "$circ_annot" != "" ]
	then
		CIRI-long collapse -i ${out_dir_call}/${f_name}/${f_name}.lst \
			-o ${out_dir_collapse}/${f_name} \
			-p $f_name \
			-r $genome_fasta \
			-a $genome_gtf \
			-c $circ_annot \
			-t $thread
	else
		CIRI-long collapse -i ${out_dir_call}/${f_name}/${f_name}.lst \
			-o ${out_dir_collapse}/${f_name} \
			-p $f_name \
			-r $genome_fasta \
			-a $genome_gtf \
			-t $thread
   fi
done

echo "All the process is done!"
