#####################################
#   AUTHOR: ALEXANDER LACHMANN      #
#   DATE: 5/16/2019                 #
#   Mount Sinai School of Medicine  #
#####################################

mkdir tools

conda config --add channels conda-forge
conda config --add channels bioconda

conda create -y -n alignment
conda activate alignment

# install aligners
conda install -y salmon

conda install -y kallisto

conda install -y hisat2

conda install -y star

conda install -y bwa

# install some other tools
conda install -y sra-tools
conda install -y htseq
conda install -y subread
conda install -y trimmomatic
