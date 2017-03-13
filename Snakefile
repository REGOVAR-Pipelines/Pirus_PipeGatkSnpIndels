import glob, os


configfile : "config.yml"
workdir: config["workdir"]




rule clean:
    shell:
        "rm -f sim.*"

# ----------------------------------------------------
# Common
# ----------------------------------------------------


# rule unzip:
#     input:
#         "{filename}.gz"
#     output:
#         "{filename}"
#     shell:
#         "gzip -d {input}"


# rule bgzip:
#     input:
#         "{filename}.vcf"
#     output:
#         "{filename}.vcf.gz"
#     shell:
#         "bgzip {input}"

# rule vcf_index:
#     input:
#         "{filename}.vcf.gz"
#     output:
#         "{filename}.vcf.gz.tbi"
#     shell:
#         "tabix -p vcf {input}"

rule fq_to_fastq:
    input:
        "{filename}.fq"
    output:
        "{filename}.fastq"
    shell:
        "mv {input} {output}"

# ----------------------------------------------------
# Alignment
# ----------------------------------------------------


rule qcreport:
    input: 
        "{filename}.fastq.gz"
    output:
        "html/{filename}_fastqc.html"
    log:
        "html/{filename}.report.log"
    shell:
        "mkdir -p html;"
        "fastqc {input} -o html 2>&1 > {log}"




#TODO : il faut gérer plus proprement les header du sam généré,
#       gatk est très fragile avec ça. ID doit être égale à bwa, tous doivent être renseignés
rule alignment:
    input:
        fwd = "{sample}_R1.fastq",
        rev = "{sample}_R2.fastq"
    output:
        "{sample}.sam"
    threads:128 
    log:
        "{sample}.bwa.log"
    shell:
        "bwa mem -R '@RG\tID:bwa\tSM:{wildcards.sample}\tLB:library1\tPL:ILLUMINA' -t {threads} {config[ref]} {input.fwd} {input.rev} > {output}  2> {log}"



rule sorting_sam:
    input:
        "{sample}.sam"
    output:
        "{sample}.sorted.sam"
    threads:
        128
    log:
        "{sample}.sorted.log"
    shell:
        #"picard SortSam INPUT={input} OUTPUT={output} SORT_ORDER=coordinate 2> {log}"
        "samtools sort {input} --threads {threads} -o {output}"


rule mark_duplicate:
  input:
      "{sample}.sorted.sam"
  output:
      "{sample}.markDuplicate.sam",
      "{sample}.markDuplicate.metrics"
  log:
      "{sample}.markDuplicate.log"
  shell:
      "picard MarkDuplicatesWithMateCigar INPUT={input} OUTPUT={output[0]} M={output[1]} 2> {log}"


rule sam_to_bam:
    input:
        "{filename}.sam"
    output:
        bam = "{filename}.bam"
    shell:
        "samtools view  -Sbh {input} > {output.bam}"


rule index_bam:
    input:
        "{filename}.bam"
    output:
        "{filename}.bam.bai"
    shell:
        "samtools index {input} {output}"



# ----------------------------------------------------
# BQSR : Base Recalibration
# ----------------------------------------------------


# TODO : add 1000G, dbSNP, ...
rule bam_RealignerTargetCreator:
    input:
        "{sample}.markDuplicate.bam",
        "{sample}.markDuplicate.bam.bai"
    output:
        "{sample}.indel.intervals"
    threads:
        128
    shell:
        "gatk -Xms{config[ram]}G -Xmx{config[ram]}G -T RealignerTargetCreator -nt {threads} -R {config[ref]} -I {input[0]} -o {output}"


# TODO : add 1000G, dbSNP, ...
rule bam_IndelRealigner:
    input:
        bam = "{sample}.markDuplicate.bam",
        bai = "{sample}.markDuplicate.bam.bai",
        itv = "{sample}.indel.intervals"
    output:
        "{sample}.realigned.bam"
    shell:
        "gatk -Xms20G -Xmx20G -T IndelRealigner -R {config[ref]} -I {input.bam} -targetIntervals {input.itv} -o {output}"


rule bam_BaseRecalibrator:
    input:
        "{filename}.bam"
    output:
        "{filename}.table"
    threads:
        128
    shell:
        "gatk -Xms{config[ram]}G -Xmx{config[ram]}G -T BaseRecalibrator -R {config[ref]} -nct {threads} -knownSites {config[dbsnp]} -I {input} -o {output}"



rule bam_PrintRead:
    input:
        bam = "{filename}.bam",
        tbl = "{filename}.table"
    output:
        "{filename}.BQSR.bam"
    threads:
        128
    shell:
        "gatk -Xms{config[ram]}G -Xmx{config[ram]}G -T PrintReads -R {config[ref]} -nct {threads} -I {input.bam} --BQSR {input.tbl} -o {output}"



# TODO :
#  - second pass to analyze covariation remaining after recalibration
#  - Generating before/after recalibration plot : https://software.broadinstitute.org/gatk/documentation/article?id=2801









rule fasta_dict_creation:
    input:
        "{reference}.fa"
    output:
        "{reference}.dict"
    shell:
        "picard CreateSequenceDictionary REFERENCE={input} OUTPUT={output}"








# ----------------------------------------------------
# Variant Calling
# ----------------------------------------------------


rule variant_calling_samtools:
    input:
        bam = "{sample}.realigned.BQSR.bam",
        bai = "{sample}.realigned.BQSR.bam.bai"
    output:
        "{sample}.samtools.calling.vcf"
    log:
        "{sample}.samtools.calling.vcf.log"
    shell:
        "samtools mpileup -ugf {config[ref]} {input.bam} | bcftools call -vm -O v -o {output}  2> {log}"


rule variant_calling_HC:
    input:
        bam = "{sample}.realigned.BQSR.bam",
        bai = "{sample}.realigned.BQSR.bam.bai"
    output:
        "{sample}.hc.calling.vcf"
    log:
        "{sample}.hc.calling.vcf.log"
    threads:
        128
    shell:
        "gatk -Xms{config[ram]}G -Xmx{config[ram]}G -T HaplotypeCaller -R {config[ref]} -nct {threads} -I {input.bam} -o {output} \
            --emitRefConfidence GVCF --variant_index_type LINEAR \
            --variant_index_parameter 128000 \
            --dbsnp:vcfinput,VCF {config[dbsnp]} "

# TODO : UnifiedGenotyper variant calling
# TODO : Platipus variant calling



# ----------------------------------------------------
# VCF Merging/joining
# ----------------------------------------------------
# Mergin each sample's gvcf into one gvcf
rule variant_CombineGVCFs:
    # All (and only) gvcf produced by HC can be used by this tool
    input:
        expand("{sample}.hc.calling.vcf", sample=config['samples'])
    output:
        "all.hc.combine.vcf"
    params:
        input_files = " ".join(["--variant {}".format(vcf) for vcf in os.listdir(".") if vcf.endswith(".hc.calling.vcf")])
    shell:
        "gatk -Xms{config[ram]}G -Xmx{config[ram]}G -T CombineGVCFs -R {config[ref]} {params.input_files} -o {output}"


# Performs the multi-sample joint aggregation 
# will combine all spanning records, produce correct genotype likelihoods, re-genotype the newly merged record
# Need gvcf produced by HC with '-ERC GVCF' mode
rule variant_GenotypeGVCFs:
    input:
        "all.hc.combine.vcf"
    output:
        "all.vcf"
    threads:
        128
    shell:
        "gatk -Xms{config[ram]}G -Xmx{config[ram]}G -T GenotypeGVCFs -R {config[ref]} -nt {threads} --variant {input} -o {output}"








# ----------------------------------------------------
# VQSR : Variants Recalibration
# ----------------------------------------------------


rule build_snp_recalibration:
    input:
        "all.vcf"
    output:
        recal  = "all.snp.recal",
        tranch = "all.snp.tranches"
    threads:
        128
    shell:
        "gatk -Xms{config[ram]}G -Xmx{config[ram]}G -T VariantRecalibrator -mode SNP  \
            -R {config[ref]} -nt {threads} -input {input} \
            -resource:hapmap,known=false,training=true,truth=true,prior=15.0 {config[hapmap]} \
            -resource:omni,known=false,training=true,truth=true,prior=12.0 {config[omni]} \
            -resource:1000G,known=false,training=true,truth=false,prior=10.0 {config[_1000g]} \
            -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {config[dbsnp]} \
            -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an InbreedingCoeff \
            -recalFile {output.recal} \
            -tranchesFile {output.tranch} -rscriptFile output.plots.R"
        # 
        # 
        # -resource:1000G,known=false,training=true,truth=false,prior=10.0 ${OTG} \


rule apply_snp_recalibration:
    input:
        vcf   = "all.vcf",
        recal = "all.snp.recal",
        trch  = "all.snp.tranches"
    output:
        "all.snp.recal.vcf"
    threads:
        128
    shell:
        "gatk -Xms{config[ram]}G -Xmx{config[ram]}G -T ApplyRecalibration -mode SNP \
            -R {config[ref]} -nct {threads} -input {input.vcf} \
            --ts_filter_level 99.5 \
            -recalFile {input.recal} \
            -tranchesFile {input.trch} \
            -o {output}"



rule build_indels_recalibration:
    input:
        "all.vcf"
    output:
        recal  = "all.indels.recal",
        tranch = "all.indels.tranches"
    threads:
        128
    shell:
        "gatk -Xms{config[ram]}G -Xmx{config[ram]}G -T VariantRecalibrator -mode INDEL \
            -R {config[ref]} -nct {threads} -input {input} \
            -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {config[dbsnp]} \
            -resource:mills,known=false,training=true,truth=true,prior=12.0 {config[mills]} \
            -an FS -an MQ -an MQRankSum -an QD -an ReadPosRankSum --maxGaussians 4\
            -recalFile {output.recal} \
            -tranchesFile {output.tranch}"


rule apply_indels_recalibration:
    input:
        vcf   = "all.snp.recal.vcf",
        recal = "all.indels.recal",
        trch  = "all.indels.tranches"
    output:
        "all.VQSR.vcf"
    threads:
        128
    shell:
        "gatk -Xms{config[ram]}G -Xmx{config[ram]}G -T ApplyRecalibration -mode INDEL  \
            -R {config[ref]} -nct {threads} -input {input.vcf} \
            --ts_filter_level 99.0 \
            -recalFile {input.recal} \
            -tranchesFile {input.trch} \
            -o {output}"



rule filter_vqsr:
    input:
        vcf   = "all.VQSR.vcf"
    output:
        "all.final.vcf"
    threads:
        128
    shell:
        "gatk -Xms{config[ram]}G -Xmx{config[ram]}G -T SelectVariants -env -R {config[ref]} -nct {threads} -V {input.vcf} -o {output}"


# ----------------------------------------------------
# Annotations
# ----------------------------------------------------


rule annotation_snpeff:
    input:
        "{sample}.vcf"
    output:
        "{sample}.snpeff.vcf"
    log:
        "{sample}.snpeff.log"
    shell:
        "snpEff -Xmx{config[ram]}g -v {config[snpEff_ref]} {input} > {output} "




