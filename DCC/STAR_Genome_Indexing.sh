OutDir=./genome_index
grnome_fasta=./genome_ref/GRCh38.primary_assembly.genome.fa
grnome_annotation=./genome_ref/gencode.v38.primary_assembly.annotation.gff3

STAR --runThreadN 16 \
     --runMode genomeGenerate \
     --genomeDir $OutDir \
     --genomeFastaFiles $grnome_fasta \
     --sjdbGTFfile $grnome_annotation \
     --sjdbOverhang 150\
