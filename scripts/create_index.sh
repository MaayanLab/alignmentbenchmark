#####################################
#   AUTHOR: ALEXANDER LACHMANN      #
#   DATE: 5/16/2019                 #
#   Mount Sinai School of Medicine  #
#####################################

num_threads=10

# this script will download all files needed to build indeces for alignment and transcript quantification
# downloading and building all the files will take a while so I recommend getting lunch in the mean time

# a function to upload the index to the cloud, this will be helpful if we don't want to rebuilt the index all the time. Downloading is faster
# if computation is in AWS amazon instances no additional cost is incurred and download speed is high, especially for high memory instances (~ 150MB/sec)
# this will require a AWS account and the instance needs AWS credentials
saveIndexS3(){
    tar cf - ${1} | pigz > ${1}.tar.gz
    aws s3 cp ${1}.tar.gz s3://${3}/${2}.tar.gz
    rm ${1}.tar.gz
}

# 1) Download all data required for alignment or transcript quantification
#    We are using the Ensembl genomes and annotations in this project

conda activate alignmnet-benchmark

mkdir reference
mkdir -p index/salmon
mkdir -p index/kallisto
mkdir -p index/star/human_96
mkdir -p index/star/mouse_96
mkdir -p index/hisat2/human_96
mkdir -p index/hisat2/mouse_96
mkdir -p index/bwa/human_96
mkdir -p index/bwa/mouse_96

# set number of open files, needed for high thread count
ulimit -n 512

# fetch 96 release cDNA sequences.
# This is the reference used by the transcript quantification methods
# this will take some time depending on the network speed
# the downloads are forked, sometimes the FTP server of ensembl is slow
curl ftp://ftp.ensembl.org/pub/release-96/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz -o reference/human_cdna_96.fa.gz &
curl ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz -o reference/mouse_cdna_96.fa.gz &

# fetch 96 whole genome sequences
# This is the reference used by true aligners such as STAR and HISAT2
curl ftp://ftp.ensembl.org/pub/release-96/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz -o reference/human_dna_96.fa.gz &
curl ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz -o reference/mouse_dna_96.fa.gz &

# fetch GTF files annotationg the raw DNA sequences
# they are needed by the alignment algorithms STAR and HISAT2 to count transcript reads
curl ftp://ftp.ensembl.org/pub/release-96/gtf/homo_sapiens/Homo_sapiens.GRCh38.96.gtf.gz -o reference/human_96.gtf.gz &
curl ftp://ftp.ensembl.org/pub/release-96/gtf/mus_musculus/Mus_musculus.GRCm38.96.gtf.gz -o reference/mouse_96.gtf.gz &

# the script is now waiting for the downloads to finish before continuing
# the jobs -p command lists forked jobs (commands followed by &)
for job in `jobs -p`
do
echo $job
    wait $job || let "FAIL+=1"
done

# unpack some files that can't be used when zipped
pigz reference/human_cdna_96.fa.gz
pigz reference/human_dna_96.fa.gz
pigz reference/human_96.gtf.gz

pigz reference/human_cdna_96.fa.gz
pigz reference/mouse_dna_96.fa.gz
pigz reference/mouse_96.gtf.gz


# 2) Start building index using the different methods

# Salmon run create index
time salmon index -p $num_threads -t reference/human_cdna_96.fa -i index/salmon/salmon_human_96
saveIndexS3 "index/salmon/salmon_human_96" "salmon_human_96" "mssm-genecount-combined"

time salmon index -p $num_threads  -t reference/mouse_cdna_96.fa -i index/salmon/salmon_mouse_96
saveIndexS3 "index/salmon/salmon_mouse_96" "salmon_mouse_96" "mssm-genecount-combined"


# kallisto run create index
time kallisto index -i index/kallisto/kallisto_human_96.idx reference/human_cdna_96.fa
saveIndexS3 "index/kallisto/kallisto_human_96.idx" "kallisto_human_96.idx" "mssm-genecount-combined"

time kallisto index -i index/kallisto/kallisto_mouse_96.idx reference/mouse_cdna_96.fa
saveIndexS3 "index/kallisto/kallisto_mouse_96.idx" "kallisto_mouse_96" "mssm-genecount-combined"


# STAR run create index
# we use the genomeSAsparseD flag, it will create a smaller index that will require less memory

GENOME="reference/human_dna_96.fa"
GTF="reference/human_96.gtf"
OUTPUT="index/star/human_96"
STAR \
    --runMode genomeGenerate \
    --genomeDir $OUTPUT \
    --genomeFastaFiles $GENOME \
    --sjdbGTFfile $GTF \
    --runThreadN $num_threads \
    --genomeSAsparseD 2 \
    --genomeSAindexNbases 13
saveIndexS3 "index/star/human_96" "star_human_96" "mssm-genecount-combined"


GENOME="reference/mouse_dna_96.fa"
GTF="reference/mouse_96.gtf"
OUTPUT="index/star/mouse_96"
time STAR \
    --runMode genomeGenerate \
    --genomeDir $OUTPUT \
    --genomeFastaFiles $GENOME \
    --sjdbGTFfile $GTF \
    --runThreadN $num_threads \
    --genomeSAsparseD 2 \
    --genomeSAindexNbases 13
saveIndexS3 "index/star/mouse_96" "star_mouse_96" "mssm-genecount-combined"

# HISAT2 run create index
time hisat2-build -p $num_threads reference/human_dna_96.fa index/hisat2/human_96/human
saveIndexS3 "index/hisat2/human_96" "hisat2_human_96" "mssm-genecount-combined"

time hisat2-build -p $num_threads reference/mouse_dna_96.fa index/hisat2/mouse_96/mouse
saveIndexS3 "index/hisat2/mouse_96" "hisat2_mouse_96" "mssm-genecount-combined"


mkdir index/bwa/human_96
time bwa index -p index/bwa/human_96/human_96 reference/human_dna_96.fa
tar cf - index/bwa/human_96 | pigz > index/bwa/bwa_human_96.tar.gz
aws s3 cp index/bwa/bwa_human_96.tar.gz s3://mssm-genecount-combined/bwa_human_96.tar.gz


#time tools/bwa/bwa-mem2 index -p index/bwa-mem2/mouse_96/mouse_96 reference/mouse_dna_96.fa
#tar cf - index/bwa-mem2/mouse_96 | pigz > index/bwa-mem2/bwa_mouse_96.tar.gz
#aws s3 cp index/bwa-mem2/bwa_mouse_96.tar.gz s3://mssm-genecount-combined/bwa_mouse_96.tar.gz

curl https://s3.amazonaws.com/mssm-genecount-combined/bwa_human_96.tar.gz -o index/bwa-mem2/bwa_human_96.tar.gz
tar -I pigz -xvf index/bwa-mem2/bwa_human_96.tar.gz
rm index/bwa-mem2/bwa_human_96.tar.gz


curl https://s3.amazonaws.com/mssm-genecount-combined/bwa_mouse_96.tar.gz -o index/bwa-mem2/bwa_mouse_96.tar.gz
tar -I pigz -xvf index/bwa-mem2/bwa_mouse_96.tar.gz
rm index/bwa-mem2/bwa_mouse_96.tar.gz


# completed the genome index build

# create mapping of transcript ids to gene symbols
Rscript scripts/buildMapping.r








