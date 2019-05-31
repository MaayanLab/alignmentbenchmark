
curl http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip -o tools/trimm.zip


SRA_FILE="SRR7460337"
time downloadSRA $SRA_FILE "sradata"

java -jar tools/trimmomatic-0.39.jar -threads 16 SRR7460337_1.fastq 


java -jar ~/Trimmomatic-0.36/trimmomatic-0.39.jar PE -phred33 \
    ${f}_1.fastq.gz \
    ${f}_2.fastq.gz \
    ${f}_R1_paired.fq.gz \
    ${f}1_unpaired.fq.gz \
    ${f}_2_paired.fq.gz \
    ${f}2_unpaired.fq.gz \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:35


java -jar tools/trimmomatic-0.39.jar PE -phred33 \
    SRR7460337_1.fastq SRR7460337_2.fastq
