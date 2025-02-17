#!/bin/bash
#SBATCH -A Research_Project1 # research project to submit under.
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel test queue
#SBATCH --time=20:00:00 # Maximum wall time for the job
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors.
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=m.kouhsar@exeter.ac.uk # email address
#SBATCH --array=0-9
#SBATCH --job-name=DCC

DataDir=./RNASeq.fastp
OutDir=./Results.DCC
StarIndex=./GenomeIndex/STAR
RefDir=./GenomeRef

module load STAR
module load Python
module load SAMtools 

samples=($(ls ${DataDir}/*R1*.gz))
Num_samp=${#samples[@]}
denom_2=$(( SLURM_ARRAY_TASK_COUNT / 2 ))
window_size=$(( ( Num_samp + denom_2 ) / SLURM_ARRAY_TASK_COUNT ))

lower=$(( SLURM_ARRAY_TASK_ID * window_size ))
next=$(( SLURM_ARRAY_TASK_ID + 1 ))
upper=$(( next * window_size ))

if [ "$SLURM_ARRAY_TASK_ID" -eq "$SLURM_ARRAY_TASK_MAX" ]; then
    upper=$Num_samp
fi

char="/"
index=$(echo "${DataDir}" | awk -F"${char}" '{print NF-1}')
index=$(( index + 2 ))

for((i=$lower;i<$upper;i++))
do

    name=${samples[$i]}
    R1=$name
    R2=${name/R1/R2}
    
    sample_name=$(echo $name| cut -d'/' -f $index)
    sample_name=$(echo ${sample_name%_R1*.fastq.gz})
	
	mkdir -p ${OutDir_DCC}/star/${sample_name}/both
    out_dir_both=${OutDir_DCC}/star/${sample_name}/both/${sample_name}
    
    mkdir  -p  ${OutDir_DCC}/star/${sample_name}/m1
    out_dir_mate1=${OutDir_DCC}/star/${sample_name}/m1/${sample_name}
    
    mkdir -p ${OutDir_DCC}/star/${sample_name}/m2
    out_dir_mate2=${OutDir_DCC}/star/${sample_name}/m2/${sample_name}
    
    temp_dir=${OutDir_DCC}/star/${sample_name}/temp
	
    echo "start sample  $sample_name"
	echo -e '\n'
    echo "STAR both reads alignment"
    
    STAR --runThreadN 16 \
		 --genomeDir $StarIndex \
		 --outSAMtype BAM SortedByCoordinate \
		 --outFileNamePrefix $out_dir_both \
		 --readFilesIn $R1 $R2 \
		 --readFilesCommand zcat \
		 --outReadsUnmapped Fastx \
		 --outSJfilterOverhangMin 15 15 15 15 \
		 --alignSJoverhangMin 15 \
		 --alignSJDBoverhangMin 15 \
		 --seedSearchStartLmax 30 \
		 --outFilterMultimapNmax 20 \
		 --outFilterScoreMin 1 \
		 --outFilterMatchNmin 1 \
		 --outFilterMismatchNmax 2 \
		 --chimSegmentMin 15 \
		 --chimScoreMin 15 \
		 --chimScoreSeparation 10 \
		 --chimJunctionOverhangMin 15 \
		 --genomeLoad LoadAndKeep \
		 --limitBAMsortRAM 50000000000 \
		 --outTmpDir $temp_dir \

    echo "Samtools indexing"
    
    samtools index  ${out_dir_both}Aligned.sortedByCoord.out.bam
#mate 1
    echo "STAR mate 1 alignment"
    STAR --runThreadN 16 \
		 --genomeDir $StarIndex \
		 --outSAMtype None \
		 --outFileNamePrefix $out_dir_mate1 \
		 --readFilesIn $R1 \
		 --readFilesCommand zcat \
		 --outReadsUnmapped Fastx \
		 --outSJfilterOverhangMin 15 15 15 15 \
		 --alignSJoverhangMin 15 \
		 --alignSJDBoverhangMin 15 \
		 --seedSearchStartLmax 30 \
		 --outFilterMultimapNmax 20 \
		 --outFilterScoreMin 1 \
		 --outFilterMatchNmin 1 \
		 --outFilterMismatchNmax 2 \
		 --chimSegmentMin 15 \
		 --chimScoreMin 15 \
		 --chimScoreSeparation 10 \
		 --chimJunctionOverhangMin 15 \
		 --genomeLoad LoadAndKeep \
		 --limitBAMsortRAM 50000000000 \
		 --outTmpDir $temp_dir \

#mate 2
     echo "STAR mate2 alignment"
     
    STAR --runThreadN 16 \
		 --genomeDir $StarIndex \
		 --outSAMtype None \
		 --outFileNamePrefix $out_dir_mate2 \
		 --readFilesIn $R2 \
		 --readFilesCommand zcat \
		 --outReadsUnmapped Fastx \
		 --outSJfilterOverhangMin 15 15 15 15 \
		 --alignSJoverhangMin 15 \
		 --alignSJDBoverhangMin 15 \
		 --seedSearchStartLmax 30 \
		 --outFilterMultimapNmax 20 \
		 --outFilterScoreMin 1 \
		 --outFilterMatchNmin 1 \
		 --outFilterMismatchNmax 2 \
		 --chimSegmentMin 15 \
		 --chimScoreMin 15 \
		 --chimScoreSeparation 10 \
		 --chimJunctionOverhangMin 15 \
		 --genomeLoad LoadAndKeep \
		 --limitBAMsortRAM 50000000000 \
		 --outTmpDir $temp_dir \


#DCC
	mkdir -p ${OutDir_DCC}/${sample_name}
	cd ${OutDir_DCC}/${sample_name}

	if [ -f samplesheet ]; then
	rm samplesheet
	fi

	if [ -f mate1 ]; then
	rm mate1
	fi

	if [ -f mate2 ]; then
	rm mate2
	fi

	if [ -f bam_files ]; then
	rm bam_files
	fi

	filename1='samplesheet'
	filename2='mate1'
	filename3='mate2'
	filename4='bam_files'

	echo ${out_dir_both}Chimeric.out.junction  >> $filename1
	echo ${out_dir_mate1}Chimeric.out.junction  >> $filename2
	echo ${out_dir_mate1}Chimeric.out.junction  >> $filename3
	echo ${out_dir_both}Aligned.sortedByCoord.out.bam  >> $filename4
	echo "DCC circRNA calling"
	DCC @samplesheet \
		 -mt1 @mate1 \
		 -mt2 @mate2 \
		 -T 16 \
		 -D \
		 -R ${RefDir}/combine_repeat1.gtf \
		 -an ${RefDir}/gencode.v38.primary_assembly.annotation.gtf \
		 -Pi \
		 -F \
		 -M \
		 -Nr 1 1 \
		 -fg \
		 -G \
		 -A ${RefDir}/GRCh38.primary_assembly.genome.fa \
		 -B @bam_files
      
   rm -r ${OutDir_DCC}/${sample_name}/_tmp_DCC/
   rm -r ${OutDir_DCC}/star/${sample_name}

done

echo "Done!"



