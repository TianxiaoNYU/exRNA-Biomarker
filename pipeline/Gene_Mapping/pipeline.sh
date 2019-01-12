mkdir /home/zhaotianxiao/CRC-lncRNA/GSE57325-RKO/output

path0=/home/zhaotianxiao/CRC-lncRNA/GSE57325-RKO/output
path1=/BioII/lulab_b/shared/projects/collaboration/CRC_lncRNA/cell-line-GSE57325-RKO

mkdir $path0/02.mapping
mkdir $path0/03.counts
mkdir $path0/04.matrix

for dir in $(ls $path1/02.mapping)
do j=${dir};
echo $j 
mkdir $path0/02.mapping/${j}

#bowtie2 -p 4 -x /BioII/lulab_b/shared/genomes/human_hg38/index/transcriptome_rsem_bowtie2/lncRNA  -1 $path1/sra/${j}_1.fastq -2 $path1/sra/${j}_2.fastq -S $path0/02.mapping/${j}/${j}_lncRNA.sam

#bowtie2 -p 4 -x /BioII/lulab_b/shared/genomes/human_hg38/index/transcriptome_rsem_bowtie2/rRNA --un $path0/02.mapping/${j}/${j}_norRNA.fastq -1 $path1/sra/${j}_1.fastq -2 $path1/sra/${j}_2.fastq -S $path0/02.mapping/${j}/${j}_rRNA.sam

mkdir $path0/02.mapping/${j}/STAR
cd $path0/02.mapping/${j}/STAR

STAR --genomeDir /BioII/lulab_b/shared/genomes/human_hg38/index/STAR_hg38_index --readFilesIn $path1/sra/${j}_1.fastq $path1/sra/${j}_2.fastq --sjdbGTFfile /BioII/lulab_b/shared/genomes/human_hg38/anno/gtf/lncRNA.gencode27.gtf --runThreadN 4 --outSAMtype BAM SortedByCoordinate

#STAR --genomeDir /BioII/lulab_b/shared/genomes/human_hg38/index/STAR_hg38_index --readFilesIn $path1/02.mapping/${j}/bowtie2_rm_rRNA/${j}.norRNA.1.fq $path1/02.mapping/${j}/bowtie2_rm_rRNA/${j}.norRNA.2.fq --sjdbGTFfile /BioII/lulab_b/shared/genomes/human_hg38/anno/gtf/lncRNA.gencode27.gtf --runThreadN 4 --outSAMtype BAM SortedByCoordinate

mkdir $path0/03.counts/${j}
featureCounts  -t exon -g transcript_id -s 1 -a /BioII/lulab_b/shared/genomes/human_hg38/anno/gtf/lncRNA.gencode27.gtf -o $path0/03.counts/${j}/${j}.lncRNA.featureCounts.counts $path0/02.mapping/${j}/STAR/Aligned.sortedByCoord.out.bam

cut -f 7 $path0/03.counts/${j}/${j}.lncRNA.featureCounts.counts > $path0/04.matrix/${j}_temp

echo "${j} is finished!!"

done;

