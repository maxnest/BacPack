FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update -y
RUN apt-get install -y build-essential
RUN apt-get install -y wget unzip libgl1-mesa-glx libegl1-mesa libxrandr2 libxrandr2 libxss1 libxcursor1 libxcomposite1 libasound2 libxi6 libxtst6 libfindbin-libs-perl libgomp1 libdatetime-perl libxml-simple-perl libdigest-md5-perl git default-jre bioperl vim cpanminus

SHELL ["/bin/bash", "--login", "-c"]

RUN mkdir Soft

# ANACONDA #
RUN wget https://repo.anaconda.com/archive/Anaconda3-2022.10-Linux-x86_64.sh -O ~/anaconda.sh && \
    /bin/bash ~/anaconda.sh -b -p /opt/conda && \
    rm ~/anaconda.sh && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc

# FASTQC #
RUN cd Soft && wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip && \
    unzip fastqc_v0.12.1.zip && rm fastqc_v0.12.1.zip

# FastP #
RUN cd Soft && wget http://opengene.org/fastp/fastp && chmod a+x ./fastp

# SPAdes #
RUN cd Soft && wget https://cab.spbu.ru/files/release3.15.4/SPAdes-3.15.4-Linux.tar.gz && \
    tar -xzvf SPAdes-3.15.4-Linux.tar.gz && rm SPAdes-3.15.4-Linux.tar.gz

# QUAST #
RUN cd Soft && wget https://sourceforge.net/projects/quast/files/quast-5.2.0.tar.gz/download && \
    mv download quast.tar.gz && tar -xzvf quast.tar.gz && rm quast.tar.gz

# FastANI #
RUN cd Soft && wget --quiet https://github.com/ParBLiSS/FastANI/releases/download/v1.33/fastANI-Linux64-v1.33.zip && \
    unzip fastANI-Linux64-v1.33.zip && rm fastANI-Linux64-v1.33.zip

# Prodigal #
RUN cd Soft && wget --quiet https://github.com/hyattpd/Prodigal/releases/download/v2.6.3/prodigal.linux && \
    chmod +x prodigal.linux

# Prokka #
RUN cpanm Test::RequiresInternet && cpanm Bio::Perl
RUN cd Soft && git clone https://github.com/tseemann/prokka.git && /Soft/prokka/bin/prokka --setupdb
RUN echo "### Prokka ###" >> ~/.bashrc && echo 'export PATH="/Soft/prokka:$PATH"' >> ~/.bashrc && \
    echo 'export PATH="/Soft/prokka/bin:$PATH"' >> ~/.bashrc && \
    echo 'export PATH="/Soft/prokka/binaries:$PATH"' >> ~/.bashrc && \
    echo 'export PATH="/Soft/prokka/db:$PATH"' >> ~/.bashrc && \
    echo 'export PATH="/Soft/prokka/binaries/linux:$PATH"' >> ~/.bashrc
RUN cpanm Bio::SearchIO::hmmer3
RUN cd Soft && wget --quiet https://ftp.ncbi.nih.gov/toolbox/ncbi_tools/converters/by_program/tbl2asn/linux64.tbl2asn.gz && \
    gunzip linux64.tbl2asn.gz && mv linux64.tbl2asn tbl2asn && chmod +x tbl2asn && \
    echo "### tbl2asn ###" >> ~/.bashrc && echo 'export PATH="/Soft:$PATH"' >> ~/.bashrc && \
    rm /Soft/prokka/binaries/linux/tbl2asn

# CryProcessor #
RUN cd Soft && git clone https://github.com/lab7arriam/cry_processor && /opt/conda/bin/python -m pip install biopython && \
    echo "### CryProcessor ###" >> ~/.bashrc && echo 'export PATH="/Soft/cry_processor:$PATH"' >> ~/.bashrc

# RabbitQCPlus #
RUN apt-get install libz-dev
RUN cd Soft && git clone https://github.com/RabbitBio/RabbitQCPlus.git && cd RabbitQCPlus && make -j4

# Flye #
RUN cd Soft && git clone https://github.com/fenderglass/Flye && cd Flye && make

# MEDAKA #
RUN apt-get update && apt-get install --no-install-recommends -y make cmake file bzip2 gcc g++ python zlib1g-dev wget curl libffi-dev libncurses5-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev python3-all-dev python3.8-venv python3-virtualenv virtualenv git-lfs && apt-get clean && apt-get autoclean && rm -rf /var/lib/apt/lists/*
RUN git lfs install --skip-repo
RUN cd Soft && git clone https://github.com/nanoporetech/medaka.git && cd medaka && make install

# BWA #
RUN cd Soft && git clone https://github.com/lh3/bwa.git && cd bwa && make

# SAMtools #
RUN cd Soft && wget --quiet https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2 && \
    tar -xvf samtools-1.16.1.tar.bz2 && rm samtools-1.16.1.tar.bz2 && cd samtools-1.16.1 && \
    ./configure --prefix=/Soft/samtools-1.16.1/ && make && make install && \
    echo "### SAMtools (v1.16.1) ###" >> ~/.bashrc && echo 'export PATH="/Soft/samtools-1.16.1:$PATH"' >> ~/.bashrc

# PILON #
RUN cd Soft && wget --quiet https://github.com/broadinstitute/pilon/releases/download/v1.24/pilon-1.24.jar

# MMseq2 #
RUN cd Soft && wget --quiet https://mmseqs.com/latest/mmseqs-linux-avx2.tar.gz && tar xvfz mmseqs-linux-avx2.tar.gz && \
    echo "### MMseqs2 ###" >> ~/.bashrc && echo 'export PATH="/Soft/mmseqs/bin:$PATH"' >> ~/.bashrc

# Snakemake #
RUN conda create -n snakemake_env && conda activate snakemake_env && conda install -c conda-forge mamba && mamba install -c conda-forge -c bioconda snakemake && python -m pip install biopython

# BacPack and programs installed using yaml-files from it #

# BacPack #
RUN cd Soft && git clone https://github.com/maxnest/BacPack

# IPG #
RUN cd Soft/BacPack/resources/IPG && gunzip *.gz && cat *.fasta > IPG_Bt_sequence.fasta

# BUSCO #
RUN /opt/conda/bin/conda init bash
RUN conda env create -f /Soft/BacPack/envs/BUSCO_env.yaml && \
    cd Soft/BacPack/resources/BUSCO_db && tar -xzvf bacilli_odb10.2021-02-23.tar.gz && tar -xzvf bacillales_odb10.2021-02-23.tar.gz

# BtToxin_Digger #
RUN conda env create -f /Soft/BacPack/envs/bttoxin_digger_env.yaml 

# DeepBGC #
RUN conda env create -f /Soft/BacPack/envs/DeepBGC_env.yaml

# antiSMASH #
RUN conda env create -f /Soft/BacPack/envs/antiSMASH_env.yaml
RUN cd Soft && wget --quiet https://dl.secondarymetabolites.org/releases/6.1.1/antismash-6.1.1.tar.gz && \
    tar -zxf antismash-6.1.1.tar.gz && rm antismash-6.1.1.tar.gz && \
    conda activate antismash_env && pip install ./antismash-6.1.1 && \
    conda install meme==4.11.2 && download-antismash-databases

# CheckM #
RUN conda env create -f /Soft/BacPack/envs/CheckM_env.yaml
RUN mkdir /Soft/CheckM && cd /Soft/CheckM && wget --quiet https://zenodo.org/record/7401545/files/checkm_data_2015_01_16.tar.gz?download=1 && \
    mv 'checkm_data_2015_01_16.tar.gz?download=1' checkm_data_2015_01_16.tar.gz && \
    tar -xzvf checkm_data_2015_01_16.tar.gz && rm checkm_data_2015_01_16.tar.gz
RUN conda activate checkm && checkm data setRoot /Soft/CheckM/

