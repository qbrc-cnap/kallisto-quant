From debian:stretch

RUN apt-get update && \
  apt-get install -y build-essential \
    wget \
    zlib1g-dev \
    libcurl4-openssl-dev \
    libncurses-dev \
    libbz2-dev \
    liblzma-dev \
    libxml2-dev \
    libgfortran3 \
    libopenblas-dev \
    libcairo2-dev \
    zip \
    python3-dev \
    python3-pip \
    default-jre \
    pandoc \
    git

# Install all software under /opt/software:
RUN mkdir -p /opt/software 

# Get the miniconda installer script and run it in silent mode:
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /opt/software/miniconda.sh && \
  bash /opt/software/miniconda.sh -b -p /opt/software/miniconda
ENV PATH="/opt/software/miniconda/bin:${PATH}"

# install kallisto and sleuth via conda
RUN conda create --name r36 r=3.6.0
SHELL ["conda", "run", "-n", "r36", "/bin/bash", "-c"]
RUN conda config --add channels defaults && \
  conda config --add channels bioconda && \
  conda config --add channels conda-forge && \
  #conda install -y bioconductor-deseq2 && \
  conda install -y bioconductor-tximport && \
  conda install -y kallisto 
  
  #&& \
  #conda install -y r-sleuth && \
  #conda install -y -c r r-ggplot2 && \
  #conda install -y -c r r-gridextra

# install trimmomatic
RUN cd /opt/software/ && \
    wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip && \
    unzip Trimmomatic-0.39.zip

# Install fastQC:
RUN cd /opt/software && \
  wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.8.zip && \
  unzip fastqc_v0.11.8.zip && \
  chmod +x /opt/software/FastQC/fastqc
ENV PATH="/opt/software/FastQC:${PATH}"

# Install libStatGen and fastQValidator:
RUN cd /opt/software && \
  git clone https://github.com/statgen/libStatGen.git && \
  cd libStatGen && \
  make && \
  cd .. && \
  git clone https://github.com/statgen/fastQValidator.git && \
  cd fastQValidator && \
  make && \
  make install

# Install some Python3 libraries:
ADD requirements.txt /opt/software/
RUN pip3 install -r /opt/software/requirements.txt

# Install multiqc:
RUN pip3 install multiqc

# Install Jinja2 for templating.  Likely installed with other python packages, but just to be sure:
RUN pip3 install Jinja2

# for multiqc to read commandline:
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8

ADD concat_tpm.py /opt/software/

# The scripts for performing pre-workflow checks:
RUN mkdir /opt/software/precheck
ADD perform_precheck.py /opt/software/precheck/
ADD check_fastq.py /opt/software/precheck/

# The script for generating the markdown report:
#ADD generate_report.py /usr/local/bin/
#RUN chmod +x /usr/local/bin/generate_report.py

# The script for organizing files prior to compression:
#ADD organize_report.py /opt/software

# The scripts for preparing files for sleuth:
#ADD move_files.py /opt/software/
#ADD create_sleuth_annotation_file.py /opt/software/
#ADD pivot_normalized_counts.py /opt/software/
#ADD split_pdx_genes.py /opt/software/

# the actual script to run sleuth and tximport:
#ADD sleuth.R /opt/software
COPY run_tximport.R /opt/software

# Add the report template and styling:
#RUN mkdir /opt/report
#ADD report.md /opt/report/
#ADD report.css /opt/report/

ENTRYPOINT ["/bin/bash"]

