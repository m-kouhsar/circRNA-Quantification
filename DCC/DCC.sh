#!/bin/bash
#SBATCH -A  # research project to submit under.
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel test queue
#SBATCH --time=20:00:00 # Maximum wall time for the job
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors.
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=m.kouhsar@exeter.ac.uk # email address

genome_index_dir=./genome_index1
data_dir=./RNA_Seq
repeat_gtf=./DCC/combine_repeat1.gtf
genome_annotation=./references/genome_anno/gencode.v38.primary_assembly.annotation.gff3
genome_fasta=./references/genome_fasta/GRCh38.primary_assembly.genome.fa

module load Salmon
module load STAR
module load Python
module load SAMtools

for i in ${data_dir}*r1*.gz
do

    R1=$i
    R2=${i/r1/r2}
    R2=${R2/val_1/val_2}
    sample_name=$(echo $i| cut -d'/' -f 7)
    sample_name=$(echo $sample_name | awk -F'.' '{print $1}')
    mkdir -p ${data_dir}star/both/
    out_dir_both=${data_dir}star/both/${sample_name}
    mkdir  -p  ${data_dir}star/m1/
    out_dir_mate1=${data_dir}star/m1/${sample_name}
    mkdir -p ${data_dir}star/m2/
    out_dir_mate2=${data_dir}star/m2/${sample_name}
    temp_dir=${data_dir}star/temp
    out_salmon=$salmon_dir${sample_name}

    echo "start sample  $sample_name"

    #alighmnet both reads

    STAR --runThreadN 16 \
		 --genomeDir $genome_index_dir \
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

    echo "STAR both reads alignment was done"
    rm -r $temp_dir
    samtools index  ${out_dir_both}Aligned.sortedByCoord.out.bam
#mate 1

    STAR --runThreadN 16 \
		 --genomeDir $genome_index_dir \
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

    echo "STAR mate1 alignment was done"
    rm -r $temp_dir
#mate 2
 
    STAR --runThreadN 16 \
		 --genomeDir $genome_index_dir \
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

    echo "STAR mate2 alignment was done"
    rm -r $temp_dir

#DCC
	mkdir -p ${data_dir}DCC/${sample_name}
	cd ${data_dir}DCC/${sample_name}

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
	
	DCC @samplesheet \
		 -mt1 @mate1 \
		 -mt2 @mate2 \
		 -T 16 \
		 -D \
		 -R $repeat_gtf \
		 -an $genome_annotation \
		 -Pi \
		 -F \
		 -M \
		 -Nr 1 1 \
		 -fg \
		 -G \
		 -A $genome_fasta \
		 -B @bam_files
      
   rm -r ${data_dir}DCC/${sample_name}/_tmp_DCC/
   rm -r ${data_dir}star/
   #rm $R1
   #rm $R2
   echo "${sample_name} is done"
done

echo "All process has been done!"



