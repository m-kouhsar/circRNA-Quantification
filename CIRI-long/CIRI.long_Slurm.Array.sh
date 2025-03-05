#!/bin/bash
#SBATCH -A Research_Project-MRC164847 # research project to submit under.
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq
#SBATCH --time=120:00:00 # Maximum wall time for the job
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors.
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=m.kouhsar@exeter.ac.uk # email address
#SBATCH --array=0-5

#########################################################################################
#########################################################################################

out_dir=./Results
fastq_dir=./Raw
genome_fasta=./GRCh38.p14.genome.fa
genome_gtf=./gencode.v47.chr_patch_hapl_scaff.annotation.gtf
circ_annot=./circAtlas_human_bed_v3.0.bed
file_identifier=fastq   # fastq fq.gz fastq.gz
call_only=no
thread=15

#######################################################################################
#######################################################################################

out_dir_call=${out_dir}/CIRI.long.Call
out_dir_collapse=${out_dir}/CIRI.long.Collapse
fastq_files=(${fastq_dir}/*.${file_identifier})
circ_annot=$(echo $circ_annot | xargs)
call_only=$(echo $call_only | xargs)
call_only=$(echo "$call_only" | tr '[:upper:]' '[:lower:]')

Num_samp=${#fastq_files[@]}

# Ref: https://stackoverflow.com/questions/2395284/round-a-divided-number-in-bash
denom_2=$(( SLURM_ARRAY_TASK_COUNT / 2 ))
window_size=$(( ( Num_samp + denom_2 ) / SLURM_ARRAY_TASK_COUNT ))

lower=$(( SLURM_ARRAY_TASK_ID * window_size ))

fastq_files1=(${fastq_files[@]:${lower}:${window_size}})

if [ "$SLURM_ARRAY_TASK_ID" -eq "$SLURM_ARRAY_TASK_MAX" ]
then
    fastq_files1=(${fastq_files[@]:$lower})
fi

echo Output directory: $out_dir
echo Fastq files directory: $fastq_dir
echo Genome fasta file: $genome_fasta
echo Genome fasta file: $genome_gtf
echo "Genome circRNA anootation file (optional): $circ_annot"
echo "Run Call mode only? $call_only"
echo Number of CPU cores: $thread
echo Start array index: $SLURM_ARRAY_TASK_MIN
echo End array index : $SLURM_ARRAY_TASK_MAX
echo numer of arrays: $SLURM_ARRAY_TASK_COUNT
echo current array index: $SLURM_ARRAY_TASK_ID
echo Total number of samples: $Num_samp
echo number of samples in the current array: ${#fastq_files1[@]}

echo "##########################################################################"
echo -e '\n'

mkdir -p $out_dir_call
mkdir -p $out_dir_collapse

j=0

for f in ${fastq_files1[@]}
do
	f_name=$(basename $f)
	f_name=${f_name%".fastq"}

  j=$(( j + 1 ))

	echo "*********************************************************************************"
	echo "*********************************************************************************"
	echo "                Running CIRI-long Call on sample $j: ${f_name}                   "
	echo "*********************************************************************************"
	echo "*********************************************************************************"
	if [ "$circ_annot" != "" ]
	then
			CIRI-long call -i $f \
			-o ${out_dir_call}/${f_name} \
			-r $genome_fasta \
			-p $f_name \
			-a $genome_gtf \
			-c $circ_annot \
			-t $thread
	else
		CIRI-long call -i $f \
			-o ${out_dir_call}/${f_name} \
			-r $genome_fasta \
			-p $f_name \
			-a $genome_gtf \
			-t $thread
	fi

	if [ -d "${out_dir_call}/${f_name}/tmp" ]
	then
		rm -r ${out_dir_call}/${f_name}/tmp
	fi

	if [ "$call_only" != "yes" ]
	then
		echo $f_name ${out_dir_call}/${f_name}/${f_name}.cand_circ.fa > ${out_dir_call}/${f_name}/${f_name}.lst
		echo "*********************************************************************************"
		echo "*********************************************************************************"
		echo "                Running CIRI-long Collapse on sample $j: ${f_name}                   "
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

		if [ -d "${out_dir_collapse}/${f_name}/tmp" ]
		then
			rm -r ${out_dir_collapse}/${f_name}/tmp
		fi
   fi
done

echo "All the process is done!"
