# Pirus_PipeGatkSnpIndels
The gatk pipe for snp/indels (using haplotypecaller) implementation for Pirus

# Howto use it in a Pirus container
## Requirement
 * You need LXD on your computer to create it
 * You should read the official doc of Pirus and LXD (be sure that the LXD bridge is well configured. Your container will need to access to internet).

##Instructions
    # create a container
    lxc launch images:ubuntu/xenial GatkHaplotypeCaller
    # configure it
    lxc exec GatkHaplotypeCaller -- /bin/bash

    # following directories are mandatory
    mkdir -p /pipeline/{run,inputs,outputs,logs,db,conda}

    # need curl if you want to notify server with the progress of your run
    apt update
    apt install wget curl jq nano git --fix-missing

    # Install miniconda3
    cd /pipeline/conda
    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
    chmod +x Miniconda3-latest-Linux-x86_64.sh
    # Proceed to the installation ("yes" and default option everywhere)
    ./Miniconda3-latest-Linux-x86_64.sh

    # Set the conda environment for the pipe
    git clone https://github.com/REGOVAR/Pirus_PipeGatkSnpIndels.git
    cd Pirus_PipeGatkSnpIndels
    /root/miniconda3/bin/conda env create -n gatk_snp_indels -f environment.yml
    
    # Activate GATK with your license (free)
    # 1 : Create/Login your account on their website : https://software.broadinstitute.org/gatk/download/licensing.php
    # 2 : Choose in the download section the license you need and download 
    # 3 : save the tar.bz2 from the website on your computer (by example in /tmp/ forlder)
    # 4 : put the file into the container
    # exit the container. The next command shall be run from the host
    exit 
    lxc file push /var/tmp/GenomeAnalysisTK-3.7.tar.bz2 GatkHaplotypeCaller/pipeline/conda/
    lxc exec GatkHaplotypeCaller -- /bin/bash
    cd /pipeline/conda
    # 5 : activate the conda environment
    source activate gatk_snp_indels
    # 6 : activate the license
    gatk-register GenomeAnalysisTK-3.7.tar.bz2

