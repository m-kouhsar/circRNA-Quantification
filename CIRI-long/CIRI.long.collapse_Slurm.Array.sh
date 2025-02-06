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
#SBATCH --array=0-9

#########################################################################################
#########################################################################################

out_dir=/lustre/projects/Research_Project-191391/Morteza/circRNA/Results/circAtlas_Annot
call_dir=/lustre/projects/Research_Project-191391/Morteza/circRNA/Results/CIRI.long.Call
genome_fasta=/lustre/projects/Research_Project-191391/Morteza/circRNA/Ref_Genome/GRCh38.p14.genome.fa
genome_gtf=/lustre/projects/Research_Project-191391/Morteza/circRNA/Ref_Genome/gencode.v47.chr_patch_hapl_scaff.annotation.gtf
circ_annot=/lustre/projects/Research_Project-191391/Morteza/circRNA/Ref_Genome/circAtlas_human_bed_v3.0.bed
thread=15

#######################################################################################
#######################################################################################

out_dir_collapse=${out_dir}/CIRI.long.Collapse
cd $call_dir
call_samples=($(ls -d */ | sed 's#/##'))
cd -
circ_annot=$(echo $circ_annot | xargs)

Num_samp=${#call_samples[@]}
window_size=$(( Num_samp / SLURM_ARRAY_TASK_COUNT + 1 ))

lower=$(( SLURM_ARRAY_TASK_ID * window_size ))

call_samples1=(${call_samples[@]:${lower}:${window_size}})

if [ "$SLURM_ARRAY_TASK_ID" -eq "$SLURM_ARRAY_TASK_MAX" ]
then
    call_samples1=(${call_samples[@]:$lower})
fi

echo Output directory: $out_dir
echo CIRI-long call results directory: $call_dir
echo Genome fasta file: $genome_fasta
echo Genome fasta file: $genome_gtf
echo "Genome circRNA anootation file (optional): $circ_annot"
echo Number of CPU cores: $thread
echo Start array index: $SLURM_ARRAY_TASK_MIN
echo End array index : $SLURM_ARRAY_TASK_MAX
echo numer of arrays: $SLURM_ARRAY_TASK_COUNT
echo current array index: $SLURM_ARRAY_TASK_ID
echo Total number of samples in call directory: $Num_samp
echo Window size: $window_size
echo number of samples in the current array: ${#call_samples1[@]}

echo "##########################################################################"
echo -e '\n'

mkdir -p $out_dir_collapse

for f in ${call_samples1[@]}
do
	echo "*********************************************************************************"
	echo "*********************************************************************************"
	echo "                Running CIRI-long Collapse on $f:                         "
	echo "*********************************************************************************"
	echo "*********************************************************************************"
	echo $f ${call_dir}/${f}/${f}.cand_circ.fa > ${call_dir}/${f}/${f}.lst
	if [ "$circ_annot" != "" ]
	then
		CIRI-long collapse -i ${call_dir}/${f}/${f}.lst \
			-o ${out_dir_collapse}/${f} \
			-p $f \
			-r $genome_fasta \
			-a $genome_gtf \
			-c $circ_annot \
			-t $thread
	else
		CIRI-long collapse -i ${call_dir}/${f}/${f}.lst \
			-o ${out_dir_collapse}/${f} \
			-p $f \
			-r $genome_fasta \
			-a $genome_gtf \
			-t $thread
   fi
done

echo "All the process is done!"


