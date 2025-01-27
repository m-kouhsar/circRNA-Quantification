#!/bin/bash
#SBATCH -A Research_Project-MRC164847 # research project to submit under.
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq
#SBATCH --time=10:00:00 # Maximum wall time for the job
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors.
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=m.kouhsar@exeter.ac.uk # email address
#SBATCH --array=1-24

#########################################################################################
#########################################################################################

out_dir=/lustre/projects/Research_Project-191391/Morteza/circRNA/Results
fastq_dir=/lustre/projects/Research_Project-191391/Morteza/circRNA/Raw/Rosie_longRead
genome_fasta=/lustre/projects/Research_Project-191391/Morteza/circRNA/Ref_Genome/GRCh38.p14.genome.fa
genome_gtf=/lustre/projects/Research_Project-191391/Morteza/circRNA/Ref_Genome/gencode.v47.chr_patch_hapl_scaff.annotation.gtf
thread=16
#######################################################################################
#######################################################################################

out_dir_call=${out_dir}/CIRI.long.Call.${SLURM_ARRAY_TASK_ID}
out_dir_collapse=${out_dir}/CIRI.long.Collapse.${SLURM_ARRAY_TASK_ID}
fastq_files=(${fastq_dir}/*.fastq)

Num_samp=${#fastq_files[@]}
window_size=$(( Num_samp / SLURM_ARRAY_TASK_COUNT + 1 ))

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
echo Number of CPU cores: $Num_samp
echo Start array index: $SLURM_ARRAY_TASK_COUNT
echo End array index : $SLURM_ARRAY_TASK_COUNT
echo numer of arrays: $SLURM_ARRAY_TASK_COUNT
echo current array index: $SLURM_ARRAY_TASK_ID

echo "##########################################################################"
echo -e '\n'

mkdir -p $out_dir_call
mkdir -p $out_dir_collapse

for f in ${fastq_files1[@]}
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


