#!/bin/bash
#SBATCH -A Research_Project-MRC164847 # research project to submit under.
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel test queue
#SBATCH --time=200:00:00 # Maximum wall time for the job
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors.
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=m.kouhsar@exeter.ac.uk # email address


cd /lustre/projects/Research_Project-191391/
module load Salmon
module load STAR
module load Python
module load SAMtools

# genome indexing

#STAR --runThreadN 16 \
      #--runMode genomeGenerate \
     #--genomeDir /lustre/projects/Research_Project-191391/STAR/references/genome_index/ \
     #--genomeFastaFiles /lustre/projects/Research_Project-191391/STAR/references/genome_fasta/GRCh38.primary_assembly.genome.fa \
     #--sjdbGTFfile /lustre/projects/Research_Project-191391/STAR/references/genome_anno/gencode.v38.primary_assembly.annotation.gff3 \
     #--sjdbOverhang 150\

genome_index_dir=/lustre/projects/Research_Project-191391/STAR/references/genome_index1/
data_dir=/lustre/projects/Research_Project-191391/Download_Synapse/ROSMAP_RNAseq1/

mkdir -p ${data_dir}salmon/
salmon_dir=${data_dir}salmon/

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

    #echo "R1 file is $R1"
    #echo "R2 file is $R2"
    echo "start sample  $sample_name"
    #echo "out dir both is $out_dir_both"
    #echo "out dir m1 is $out_dir_mate1"
    #echo "out dir m2 is $out_dir_mate2"
    #echo "temp dir is $temp_dir"

   #salmon 
    #salmon quant -i /lustre/projects/Research_Project-191391/salmon/ref/gencode.v38.transcripts_index \
			#-l A -1 $R1 -2 $R2 \
			#-p 16 --validateMappings \
			#--geneMap /lustre/projects/Research_Project-191391/salmon/ref/salmon_gene_transcript_v38.txt \
			#-o $out_salmon \
   
    #echo "salmon alignbment was done"

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
		 -R /lustre/projects/Research_Project-191391/DCC/combine_repeat1.gtf \
		 -an /lustre/projects/Research_Project-191391/STAR/references/genome_anno/gencode.v38.primary_assembly.annotation.gff3 \
		 -Pi \
		 -F \
		 -M \
		 -Nr 1 1 \
		 -fg \
		 -G \
		 -A /lustre/projects/Research_Project-191391/STAR/references/genome_fasta/GRCh38.primary_assembly.genome.fa \
		 -B @bam_files
      
   rm -r ${data_dir}DCC/${sample_name}/_tmp_DCC/
   rm -r ${data_dir}star/
   #rm $R1
   #rm $R2
   echo "${sample_name} is done"
done

echo "All process has been done!"



