#####################################
#   AUTHOR: ALEXANDER LACHMANN      #
#   DATE: 5/16/2019                 #
#   Mount Sinai School of Medicine  #
#####################################

downloadIndex(){
    mkdir -p index/$1
    curl https://s3.amazonaws.com/$3/$2 -o index/$1/$2
    tar -I pigz -xvf index/$1/$2
    rm index/$1/$2
}

# this script will download all files needed to align reads. Alternatively the script create_index can be run
# this would built indeces from scratch. The indeces for this project are all publically available at
# s3://mssm-genecount-combined/


# 1) Download all data required for alignment or transcript quantification
#    We are using the Ensembl genomes and annotations in this project

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
pigz -d reference/human_96.gtf.gz
pigz -d reference/mouse_96.gtf.gz

# this will download a large amout of data and will take some time, this will still be much faster than building the 
# indeces from scratch
downloadIndex "salmon" "salmon_human_96.tar.gz" "mssm-genecount-combined"
downloadIndex "kallisto" "kallisto_human_96.idx.tar.gz" "mssm-genecount-combined"
downloadIndex "star" "star_human_96.tar.gz" "mssm-genecount-combined"
downloadIndex "hisat2" "hisat2_human_96.tar.gz" "mssm-genecount-combined"
downloadIndex "bwa" "bwa_human_96.tar.gz" "mssm-genecount-combined"






