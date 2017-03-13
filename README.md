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
    mkdir -p /pipeline/{inputs,outputs,logs,db,conda}

    # need curl if you want to notify server with the progress of your run
    apt update
    apt install wget curl jq nano git --fix-missing  # optional : graphviz (if you need to draw snakemake graph)

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

    # exit the container
    exit

    # stop it and create an image
    lxc stop GatkHaplotypeCaller
    lxc publish GatkHaplotypeCaller --alias=PirusGatkHaplotypeCaller
    
    lxc image export PirusGatkHaplotypeCaller
    
    # following command must be done as root to avoid image corruption 
    # (as it will try to create symlink to computer resource in /dev folder by example)
    sudo tar xf <the_name_of_lxc_export_something_like_a8d44d24fcs...8fzef54e5>.tar.gz

    # append folowing informations into the metadata.yaml file (/!\ don't forget the "," in json dictionnary)
    sudo nano metadata.yaml

    # if json
    "pirus":
    {
        "name" : "Gatk HaplotypeCaller",
        "description" : "The GATK haplotypecaller pipeline for SNP and Indels.",
        "version": "1.0.0",
        "pirus_api": "1.0.0",
        "license" : "AGPLv3",
        "developers" : ["Anne-Sophie DENOMME-PICHON", "Jérémie ROQUET", "Olivier GUEUDELOT", "Sacha SCHUTZ"]
        "run" : "/pipeline/conda/Pirus_PipeGatkSnpIndels/run.sh",
        "inputs" : "/pipeline/inputs",
        "outputs" : "/pipeline/outputs",
        "databases" : "/pipeline/db",
        "logs" : "/pipeline/logs",
        "form" : "/pipeline/conda/Pirus_PipeGatkSnpIndels/form.json",
        "icon" : "/pipeline/conda/Pirus_PipeGatkSnpIndels/gatk-logo.png"
    }
    
    # Repackage the image in tar.xz
    sudo tar cfJ PirusGatkHaplotypeCaller.tar.xz metadata.yaml rootfs templates
    sudo rm -fr metadata.yaml rootfs templates
    sudo chown olivier:olivier PirusGatkHaplotypeCaller.tar.xz
