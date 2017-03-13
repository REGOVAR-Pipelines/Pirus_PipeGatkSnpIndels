#!/bin/bash
# coding: utf-8

# Hardcoded parameters
DB="/pipeline/db/"
OUT="/pipeline/outputs/"
IN="/pipeline/inputs/"
LOG="/pipeline/logs/"
RUN="/pipeline/conda/Pirus_PipeGatkSnpIndels"
CORES_MAX=10

# Dynamic parameters loaded from config.json
FILE1=`cat ${IN}config.json | jq ".run.file1"`
DURATION=`cat ${IN}config.json | jq ".run.duration"`
CRASH=`cat ${IN}config.json | jq ".run.crash"`
OUTPUT=`cat ${IN}config.json | jq ".run.outfilename"`
NOTIFY_URL=`cat ${IN}config.json | jq ".pirus.notify_url"`

FILE1=`sed -e 's/^"//' -e 's/"$//' <<<"$FILE1"`
OUTPUT=`sed -e 's/^"//' -e 's/"$//' <<<"$OUTPUT"`
NOTIFY_URL=`sed -e 's/^"//' -e 's/"$//' <<<"$NOTIFY_URL"`


# Init the config file of the pipe
curl -X POST -d '{"progress" : {"min":"0", "max":"5", "value":"1", "label" : "1 / 5"}}' ${NOTIFY_URL}




cd $RUN
source activate gatk_snp_indels


# Create pipe graphique
snakemake all.VQSR.vcf -np --dag | dot -Tsvg > ${OUT}graph.svg
curl -X POST -d '{"output" : {"filename":"pipe.plan.svg"}}' ${NOTIFY_URL}
curl -X POST -d '{"progress" : {"min":"0", "max":"5", "value":"2", "label" : "2 / 5"}}' ${NOTIFY_URL}


# Run the pipe
snakemake all.final.vcf --cores $CORES_MAX

# Finish
curl -X POST -d '{"progress" : {"min":"0", "max":"5", "value":"5", "label" : "5 / 5"}}' ${NOTIFY_URL}
