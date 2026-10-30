import os
STRAINSEEKER_DB = os.path.abspath(os.environ.get("PBP_STRAINSEEKER_DB", "refs/strainseeker_db"))
import pandas as pd
MyList = pd.read_csv("list/tpa.ss14.sampleList.txt", header=None)
SAMPLES = list(MyList[0])


rule all:
    input:
        expand("VCF/ss14/{sample}_HaploCaller.ss14.pass.vcf.gz",          sample=SAMPLES),
        expand("VCF/ss14/{sample}_HaploCaller.ss14.pass.norm.vcf.gz",     sample=SAMPLES),
        expand("VCF/ss14/{sample}_HaploCaller.ss14.pass.ann.vcf.gz",      sample=SAMPLES),
        expand("VCF/ss14/{sample}.ann.tsv",                               sample=SAMPLES),
        expand("consensus/ss14/{sample}.ss14.consensus.fasta",            sample=SAMPLES),
        expand("consensus/ss14/{sample}.ss14.snp_only.consensus.fasta",   sample=SAMPLES)

# 最终输出说明:
# {sample}_HaploCaller.ss14.pass.vcf.gz
#   └── SNP:   QD<2, MQ<40, FS>60, SOR>4, MQRankSum<-12.5, ReadPosRankSum<-8, DP<3
#   └── INDEL: QD<2, FS>200, ReadPosRankSum<-20, DP<3
# {sample}.ss14.consensus.fasta
#   └── 包含 SNP+INDEL，DP<3 区域用 N 遮蔽


#######################################################
##### T. Pallidum subsp SS14 reference ################
#######################################################

REF             = "refs/TPA_SS14.fasta"
HOST_REF_HUMAN = os.environ.get("PBP_HOST_REF_HUMAN", "refs/host_human_bbmap")
HOST_REF_RABBIT = os.environ.get("PBP_HOST_REF_RABBIT", "refs/host_rabbit_bbmap")
BBMAP_CACHE_SCRIPT = "script/p0_00_bbmap_node_cache.sh"
BBMAP_GUARD_SCRIPT = "script/p0_00_run_bbmap_guarded.sh"
BBMAP_ENV = os.environ.get("PBP_BBMAP_ENV", "pbp_bbmap_runtime")
BBMAP_STARTUP_TIMEOUT = os.environ.get("PBP_BBMAP_STARTUP_TIMEOUT", "900")
TRIMMOMATIC_ADAPTERS = os.environ.get("PBP_TRIMMOMATIC_ADAPTERS", "refs/TruSeq3-PE.fa")
PICARD_TMP      = "tmp"

#######################################################
##### 1. Pre-alignment QC #############################
#######################################################

rule adapter_trim:
    input:
        fq_1 = "raw_fastq/ss14/{sample}_1.fastq.gz",
        fq_2 = "raw_fastq/ss14/{sample}_2.fastq.gz"
    output:
        fq_1          = "processed_fastq/ss14/adapter_trimmed/{sample}_R1.paired_trim.fq.gz",
        fq_2          = "processed_fastq/ss14/adapter_trimmed/{sample}_R2.paired_trim.fq.gz",
        fq_1_unpaired = "processed_fastq/ss14/adapter_trimmed/{sample}_R1.unpaired.fq.gz",
        fq_2_unpaired = "processed_fastq/ss14/adapter_trimmed/{sample}_R2.unpaired.fq.gz",
        log_file      = "processed_fastq/ss14/adapter_trimmed/{sample}_trim_log.txt"
    threads: 8
    resources: mem_mb=16000
    params:
        adapters = TRIMMOMATIC_ADAPTERS
    shell:
        r"""
        conda run -n tpa_main_env \
            trimmomatic PE -threads {threads} -phred33 \
                -Xmx12g \
                -trimlog {output.log_file} \
                {input.fq_1} {input.fq_2} \
                {output.fq_1} {output.fq_1_unpaired} \
                {output.fq_2} {output.fq_2_unpaired} \
                ILLUMINACLIP:{params.adapters}:2:30:10:8:TRUE LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
        """

rule fq_rm_host_human:
    input:
        fq_1 = "processed_fastq/ss14/adapter_trimmed/{sample}_R1.paired_trim.fq.gz",
        fq_2 = "processed_fastq/ss14/adapter_trimmed/{sample}_R2.paired_trim.fq.gz"
    output:
        fq_clean  = temp("processed_fastq/ss14/host_filtered/{sample}.rm_human_R1.fq.gz"),
        fq_clean2 = temp("processed_fastq/ss14/host_filtered/{sample}.rm_human_R2.fq.gz"),
        fq_cnt    = "processed_fastq/ss14/host_filtered/{sample}.human_contamination.fq.gz"
    log:
        "logs/bbmap/ss14/{sample}.human.log"
    threads: 16
    resources: mem_mb=100000, big_mem_slots=1, sge_slots=64
    params:
        ref = HOST_REF_HUMAN,
        cache_key = "hg19_build1",
        cache_script = BBMAP_CACHE_SCRIPT,
        guard_script = BBMAP_GUARD_SCRIPT,
        bbmap_env = BBMAP_ENV,
        startup_timeout = BBMAP_STARTUP_TIMEOUT
    shell:
        r"""
        host_ref=$(bash {params.cache_script} prepare "{params.ref}" "{params.cache_key}")
        bash {params.guard_script} "{log}" "{params.startup_timeout}" \
          conda run --no-capture-output -n "{params.bbmap_env}" \
            bbmap.sh build=1 threads={threads} -Xmx64g \
                minid=0.95 maxindel=3 bandwidthratio=0.16 bandwidth=12 \
                quickmatch fast minhits=2 pigz unpigz \
                path="$host_ref" \
                in={input.fq_1} in2={input.fq_2} \
                outm={output.fq_cnt} \
                outu1={output.fq_clean} outu2={output.fq_clean2}
        """

rule fq_rm_host_rabbit:
    input:
        fq_1 = "processed_fastq/ss14/host_filtered/{sample}.rm_human_R1.fq.gz",
        fq_2 = "processed_fastq/ss14/host_filtered/{sample}.rm_human_R2.fq.gz"
    output:
        fq_clean  = temp("processed_fastq/ss14/host_filtered/{sample}.rm_rabbit_R1.fq.gz"),
        fq_clean2 = temp("processed_fastq/ss14/host_filtered/{sample}.rm_rabbit_R2.fq.gz"),
        fq_cnt    = "processed_fastq/ss14/host_filtered/{sample}.rabbit_contamination.fq.gz"
    log:
        "logs/bbmap/ss14/{sample}.rabbit.log"
    threads: 16
    resources: mem_mb=100000, big_mem_slots=1, sge_slots=64
    params:
        ref = HOST_REF_RABBIT,
        cache_key = "oryCun2_build1",
        cache_script = BBMAP_CACHE_SCRIPT,
        guard_script = BBMAP_GUARD_SCRIPT,
        bbmap_env = BBMAP_ENV,
        startup_timeout = BBMAP_STARTUP_TIMEOUT
    shell:
        r"""
        host_ref=$(bash {params.cache_script} prepare "{params.ref}" "{params.cache_key}")
        bash {params.guard_script} "{log}" "{params.startup_timeout}" \
          conda run --no-capture-output -n "{params.bbmap_env}" \
            bbmap.sh build=1 threads={threads} -Xmx64g \
                minid=0.95 maxindel=3 bandwidthratio=0.16 bandwidth=12 \
                quickmatch fast minhits=2 pigz unpigz \
                path="$host_ref" \
                in={input.fq_1} in2={input.fq_2} \
                outm={output.fq_cnt} \
                outu1={output.fq_clean} outu2={output.fq_clean2}
        """

rule fq_rm_broken:
    input:
        fq_1 = "processed_fastq/ss14/host_filtered/{sample}.rm_rabbit_R1.fq.gz",
        fq_2 = "processed_fastq/ss14/host_filtered/{sample}.rm_rabbit_R2.fq.gz"
    output:
        fq_1 = "processed_fastq/ss14/host_filtered/{sample}_R1.host_filtered.fq.gz",
        fq_2 = "processed_fastq/ss14/host_filtered/{sample}_R2.host_filtered.fq.gz"
    threads: 8
    resources: mem_mb=40000
    shell:
        r"""
        conda run -n tpa_main_env \
            repair.sh tossbrokenreads -Xmx32g \
                in={input.fq_1} in2={input.fq_2} \
                out={output.fq_1} out2={output.fq_2}
        """

rule fq_seqtk_sub:
    input:
        fq_1 = "processed_fastq/ss14/host_filtered/{sample}_R1.host_filtered.fq.gz",
        fq_2 = "processed_fastq/ss14/host_filtered/{sample}_R2.host_filtered.fq.gz"
    output:
        temp("processed_fastq/ss14/strainseeker/{sample}.sub")
    threads: 2
    resources: mem_mb=8000
    shell:
        r"""
        conda run -n tpa_main_env \
            bash -c 'seqtk sample -s100 {input.fq_1} 100000 > {output}'
        """

rule strainseeker:
    input:
        "processed_fastq/ss14/strainseeker/{sample}.sub"
    output:
        "processed_fastq/ss14/strainseeker/{sample}.seeker.txt"
    threads: 2
    resources: mem_mb=16000
    params: root = os.getcwd()
    shell:
        r"""
        mkdir -p tmp/strainseeker/{wildcards.sample}
        cd tmp/strainseeker/{wildcards.sample} && \
        conda run -n tpa_main_env \
            seeker.pl \
                -i {params.root}/processed_fastq/ss14/strainseeker/{wildcards.sample}.sub \
                -d {STRAINSEEKER_DB:q} \
                -o {params.root}/processed_fastq/ss14/strainseeker/{wildcards.sample}.seeker.txt
        """

#######################################################
##### 2. Alignment ####################################
#######################################################

rule bwa_align:
    input:
        fa           = REF,
        strainseeker = "processed_fastq/ss14/strainseeker/{sample}.seeker.txt",
        fq_1         = "processed_fastq/ss14/host_filtered/{sample}_R1.host_filtered.fq.gz",
        fq_2         = "processed_fastq/ss14/host_filtered/{sample}_R2.host_filtered.fq.gz"
    output:
        temp("bam/ss14/aln/{sample}.ss14.raw.bam")
    threads: 4
    resources: mem_mb=12000
    params:
        rg = r"@RG\tID:{sample}\tSM:{sample}\tLB:{sample}_L"
    shell:
        r"""
        conda run -n tpa_main_env \
            bash -c 'bwa mem -R '"'"'{params.rg}'"'"' -t 2 {input.fa} {input.fq_1} {input.fq_2} \
                | samtools view -@ 2 -b -o {output}'
        """

rule bam_fixmate:
    input:
        "bam/ss14/aln/{sample}.ss14.raw.bam"
    output:
        temp("bam/ss14/aln/{sample}.ss14.fixmate.bam")
    threads: 4
    resources: mem_mb=8000
    shell:
        r"""
        conda run -n tpa_main_env \
            samtools fixmate -@ {threads} -O bam {input} {output}
        """

rule bam_sort:
    input:
        "bam/ss14/aln/{sample}.ss14.fixmate.bam"
    output:
        temp("bam/ss14/aln/{sample}.ss14.sorted.bam")
    threads: 4
    resources: mem_mb=12000
    shell:
        r"""
        mkdir -p tmp
        conda run -n tpa_main_env \
            samtools sort -T tmp/{wildcards.sample} -@ {threads} -O bam -o {output} {input}
        """

rule bam_index:
    input:
        "bam/ss14/aln/{sample}.ss14.sorted.bam"
    output:
        temp("bam/ss14/aln/{sample}.ss14.sorted.bam.bai")
    threads: 4
    resources: mem_mb=8000
    shell:
        r"""
        conda run -n tpa_main_env \
            samtools index -@ {threads} {input}
        """

rule flagstat:
    input:
        bam = "bam/ss14/aln/{sample}.ss14.sorted.bam",
        bai = "bam/ss14/aln/{sample}.ss14.sorted.bam.bai"
    output:
        temp("bam/ss14/aln/{sample}.ss14.flagstat")
    threads: 4
    resources: mem_mb=8000
    shell:
        r"""
        conda run -n tpa_main_env \
            bash -c 'samtools flagstat -@ {threads} {input.bam} > {output}'
        """

rule filter_flags:
    input:
        bam      = "bam/ss14/aln/{sample}.ss14.sorted.bam",
        bai      = "bam/ss14/aln/{sample}.ss14.sorted.bam.bai",
        flagstat = "bam/ss14/aln/{sample}.ss14.flagstat"
    output:
        temp("bam/ss14/aln/{sample}.ss14.flag_filt.bam")
    threads: 4
    resources: mem_mb=8000
    shell:
        r"""
        conda run -n tpa_main_env \
            samtools view -@ {threads} -h -F 4 -F 256 -f 2 -b -o {output} {input.bam}
        """

rule index_flag_filter:
    input:
        "bam/ss14/aln/{sample}.ss14.flag_filt.bam"
    output:
        temp("bam/ss14/aln/{sample}.ss14.flag_filt.bam.bai")
    threads: 4
    resources: mem_mb=8000
    shell:
        r"""
        conda run -n tpa_main_env \
            samtools index -@ {threads} {input}
        """

rule mark_dups:
    input:
        bam = "bam/ss14/aln/{sample}.ss14.flag_filt.bam",
        bai = "bam/ss14/aln/{sample}.ss14.flag_filt.bam.bai"
    output:
        bam     = "bam/ss14/aln/{sample}.ss14.dedup.bam",
        metrics = "bam/ss14/aln/{sample}.dedup.metrics"
    threads: 4
    resources: mem_mb=16000
    log: "bam/ss14/aln/{sample}.ss14.markdup.log.txt"
    params:
        tmp = PICARD_TMP
    shell:
        r"""
        mkdir -p {params.tmp}
        conda run -n tpa_main_env \
            bash -c 'picard -Xmx12g -XX:ParallelGCThreads={threads} \
                -Djava.io.tmpdir={params.tmp} \
                MarkDuplicates \
                INPUT={input.bam} \
                OUTPUT={output.bam} \
                METRICS_FILE={output.metrics} \
                TMP_DIR={params.tmp} \
                2> {log}'
        """

rule index_dedup:
    input:
        "bam/ss14/aln/{sample}.ss14.dedup.bam"
    output:
        "bam/ss14/aln/{sample}.ss14.dedup.bam.bai"
    threads: 4
    resources: mem_mb=8000
    shell:
        r"""
        conda run -n tpa_main_env \
            samtools index -@ {threads} {input}
        """

rule bam_fixmate_b:
    input:
        bam = "bam/ss14/aln/{sample}.ss14.dedup.bam",
        bai = "bam/ss14/aln/{sample}.ss14.dedup.bam.bai"
    output:
        bam = temp("bam/ss14/aln/{sample}.ss14.matefixed.bam"),
        bai = temp("bam/ss14/aln/{sample}.ss14.matefixed.bai")
    threads: 4
    resources: mem_mb=16000
    log: "bam/ss14/aln/{sample}.ss14.matefixed.log.txt"
    params:
        tmp = PICARD_TMP
    shell:
        r"""
        mkdir -p {params.tmp}
        conda run -n tpa_main_env \
            bash -c 'picard -Xmx12g -XX:ParallelGCThreads={threads} \
                -Djava.io.tmpdir={params.tmp} \
                FixMateInformation \
                INPUT={input.bam} \
                OUTPUT={output.bam} \
                CREATE_INDEX=true \
                TMP_DIR={params.tmp} \
                2> {log}'
        """

#######################################################
##### 3. Post-alignment QC ############################
#######################################################

rule remove_mismatch:
    input:
        bam = "bam/ss14/aln/{sample}.ss14.matefixed.bam"
    output:
        bam    = temp("bam/ss14/postfilter/{sample}.ss14.rm_mismatch.bam"),
        failed = "bam/ss14/postfilter/{sample}.ss14.mismatch_reads.txt"
    threads: 2
    resources: mem_mb=16000
    shell:
        r"""
        conda run -n tpa_ngsutils_env \
            bamutils filter {input.bam} {output.bam} -failed {output.failed} -mismatch 5
        """

rule short_alignment:
    input:
        "bam/ss14/postfilter/{sample}.ss14.rm_mismatch.bam"
    output:
        temp("bam/ss14/postfilter/{sample}.ss14.minlen.bam")
    threads: 2
    resources: mem_mb=16000
    shell:
        r"""
        conda run -n tpa_ngsutils_env \
            bamutils filter {input} {output} -minlen 35
        """

rule soft_clip_a:
    input:
        "bam/ss14/postfilter/{sample}.ss14.minlen.bam"
    output:
        temp("bam/ss14/postfilter/{sample}.ss14.softclips.tmp.bam")
    threads: 2
    resources: mem_mb=32000
    shell:
        r"""
        conda run -n tpa_ngsutils_env \
            bamutils removeclipping {input} {output}
        """

rule soft_clip_b:
    input:
        "bam/ss14/postfilter/{sample}.ss14.softclips.tmp.bam"
    output:
        "bam/ss14/postfilter/{sample}.read_names_to_remove_highSoftClip.txt"
    threads: 6
    resources: mem_mb=8000
    shell:
        r"""
        conda run -n tpa_main_env \
            bash -c 'samtools view -@ 4 {input} | grep "ZC:f:" | \
            awk "{{ for (i = 1; i <= NF; i++) {{ if (\$i ~ /ZC:f:/) {{ print \$1, \$i }} }} }}" | \
            sed "s/ZC:f://" | \
            awk -v MAX_SOFTCLIP="0.05" "{{ if (\$2 > MAX_SOFTCLIP) {{ print \$1 }} }}" | \
            cut -f 1 | sort | uniq > {output}'
        """

rule soft_clip_c:
    input:
        bam      = "bam/ss14/postfilter/{sample}.ss14.minlen.bam",
        read_out = "bam/ss14/postfilter/{sample}.read_names_to_remove_highSoftClip.txt"
    output:
        bam_sc       = temp("bam/ss14/postfilter/{sample}.ss14.softclip_filtOut_reads.bam"),
        bam_filtered = temp("bam/ss14/postfilter/{sample}.ss14.rm_softclips.bam")
    threads: 2
    resources: mem_mb=40000
    params:
        tmp = PICARD_TMP
    shell:
        r"""
        mkdir -p {params.tmp}
        conda run -n tpa_main_env \
            picard -Xmx32g \
                -Djava.io.tmpdir={params.tmp} \
                FilterSamReads \
                INPUT={input.bam} OUTPUT={output.bam_sc} \
                READ_LIST_FILE={input.read_out} FILTER=includeReadList \
                TMP_DIR={params.tmp}
        conda run -n tpa_main_env \
            picard -Xmx32g \
                -Djava.io.tmpdir={params.tmp} \
                FilterSamReads \
                INPUT={input.bam} OUTPUT={output.bam_filtered} \
                READ_LIST_FILE={input.read_out} FILTER=excludeReadList \
                TMP_DIR={params.tmp}
        """

rule chimeric_aln:
    input:
        "bam/ss14/postfilter/{sample}.ss14.rm_softclips.bam"
    output:
        bam               = temp("bam/ss14/postfilter/{sample}.ss14.rm_chim.bam"),
        bam_filtOut_reads = temp("bam/ss14/postfilter/{sample}.ss14.chim_reads.bam")
    threads: 2
    resources: mem_mb=8000
    shell:
        r"""
        conda run -n tpa_main_env \
            samtools view -@ 2 -h -f 2048 -b -o {output.bam_filtOut_reads} {input}
        conda run -n tpa_main_env \
            samtools view -@ 2 -h -F 2048 -b -o {output.bam} {input}
        """

rule hard_clip:
    input:
        "bam/ss14/postfilter/{sample}.ss14.rm_chim.bam"
    output:
        bam = temp("bam/ss14/postfilter/{sample}.ss14.rm_hardclips.bam")
    threads: 4
    resources: mem_mb=8000
    shell:
        r"""
        conda run -n tpa_main_env \
            bash -c 'samtools view -@ 2 -h {input} | awk '"'"'$6 !~ /H/'"'"' | samtools view -bh -o {output.bam} -'
        """

rule MAPQ_filt:
    input:
        "bam/ss14/postfilter/{sample}.ss14.rm_hardclips.bam"
    output:
        bam = temp("bam/ss14/postfilter/{sample}.ss14.MAPQ40.bam")
    threads: 4
    resources: mem_mb=8000
    shell:
        r"""
        conda run -n tpa_main_env \
            bash -c 'samtools view -@ 2 -h -bq 40 {input} > {output.bam}'
        """

rule rm_singletone:
    input:
        "bam/ss14/postfilter/{sample}.ss14.MAPQ40.bam"
    output:
        "bam/ss14/postfilter/{sample}.ss14.final.bam"
    threads: 4
    resources: mem_mb=8000
    shell:
        r"""
        conda run -n tpa_main_env \
            samtools view -@ 2 -h -F 8 -b -o {output} {input}
        """

rule index_final_bam:
    input:
        "bam/ss14/postfilter/{sample}.ss14.final.bam"
    output:
        "bam/ss14/postfilter/{sample}.ss14.final.bam.bai"
    threads: 2
    resources: mem_mb=8000
    shell:
        r"""
        conda run -n tpa_main_env \
            samtools index -@ {threads} {input} {output}
        """

#######################################################
##### 4. Variant calling ##############################
#######################################################

rule make_gvcf:
    input:
        fa  = REF,
        bam = "bam/ss14/postfilter/{sample}.ss14.final.bam",
        bai = "bam/ss14/postfilter/{sample}.ss14.final.bam.bai"
    output:
        "VCF/ss14/{sample}_HaploCaller.ss14.raw.g.vcf.gz"
    threads: 4
    resources: mem_mb=28000
    params:
        java_opts = "-Xmx20g",
        gatk_opts = "-ploidy 1 -stand-call-conf 30 -ERC GVCF"
    shell:
        r"""
        conda run -n tpa_gatk_env \
            gatk --java-options "{params.java_opts}" HaplotypeCaller \
                -R {input.fa} {params.gatk_opts} \
                -I {input.bam} -O {output}
        """

rule make_raw_vcf:
    input:
        fa   = REF,
        gvcf = "VCF/ss14/{sample}_HaploCaller.ss14.raw.g.vcf.gz"
    output:
        "VCF/ss14/{sample}_HaploCaller.ss14.raw.vcf.gz"
    threads: 1
    resources: mem_mb=12000
    params:
        java_opts = "-Xmx8g"
    shell:
        r"""
        conda run -n tpa_gatk_env \
            gatk --java-options "{params.java_opts}" GenotypeGVCFs \
                -R {input.fa} -V {input.gvcf} -O {output}
        """

rule select_snp:
    input:
        fa  = REF,
        vcf = "VCF/ss14/{sample}_HaploCaller.ss14.raw.vcf.gz"
    output:
        "VCF/ss14/{sample}_HaploCaller.ss14.raw.snp.vcf.gz"
    threads: 1
    resources: mem_mb=12000
    params:
        java_opts = "-Xmx8g"
    shell:
        r"""
        conda run -n tpa_gatk_env \
            gatk --java-options "{params.java_opts}" SelectVariants \
                -R {input.fa} -V {input.vcf} \
                --select-type-to-include SNP \
                -O {output}
        """

rule select_indel:
    input:
        fa  = REF,
        vcf = "VCF/ss14/{sample}_HaploCaller.ss14.raw.vcf.gz"
    output:
        temp("VCF/ss14/{sample}_HaploCaller.ss14.raw.indel.vcf.gz")
    threads: 1
    resources: mem_mb=12000
    params:
        java_opts = "-Xmx8g"
    shell:
        r"""
        conda run -n tpa_gatk_env \
            gatk --java-options "{params.java_opts}" SelectVariants \
                -R {input.fa} -V {input.vcf} \
                --select-type-to-include INDEL \
                -O {output}
        """

rule filter_snp:
    input:
        fa  = REF,
        vcf = "VCF/ss14/{sample}_HaploCaller.ss14.raw.snp.vcf.gz"
    output:
        temp("VCF/ss14/{sample}_HaploCaller.ss14.filt.snp.vcf.gz")
    threads: 1
    resources: mem_mb=12000
    params:
        java_opts = "-Xmx8g"
    shell:
        r"""
        conda run -n tpa_gatk_env \
            bash -c 'gatk --java-options "{params.java_opts}" VariantFiltration \
                -R {input.fa} \
                -V {input.vcf} \
                -O {output} \
                --filter-expression "QD < 2.0"              --filter-name "SNP_QD2" \
                --filter-expression "MQ < 40.0"             --filter-name "SNP_MQ40" \
                --filter-expression "FS > 60.0"             --filter-name "SNP_FS60" \
                --filter-expression "SOR > 4.0"             --filter-name "SNP_SOR4" \
                --filter-expression "MQRankSum < -12.5"     --filter-name "SNP_MQRankSum-12.5" \
                --filter-expression "ReadPosRankSum < -8.0" --filter-name "SNP_ReadPosRankSum-8" \
                --filter-expression "DP < 3.0"              --filter-name "SNP_LowCoverage"'
        """

rule filter_indel:
    input:
        fa  = REF,
        vcf = "VCF/ss14/{sample}_HaploCaller.ss14.raw.indel.vcf.gz"
    output:
        temp("VCF/ss14/{sample}_HaploCaller.ss14.filt.indel.vcf.gz")
    threads: 1
    resources: mem_mb=12000
    params:
        java_opts = "-Xmx8g"
    shell:
        r"""
        conda run -n tpa_gatk_env \
            bash -c 'gatk --java-options "{params.java_opts}" VariantFiltration \
                -R {input.fa} \
                -V {input.vcf} \
                -O {output} \
                --filter-expression "QD < 2.0"               --filter-name "INDEL_QD2" \
                --filter-expression "FS > 200.0"             --filter-name "INDEL_FS200" \
                --filter-expression "ReadPosRankSum < -20.0" --filter-name "INDEL_ReadPosRankSum-20" \
                --filter-expression "DP < 3.0"               --filter-name "INDEL_LowCoverage"'
        """

rule select_pass_snp:
    input:
        fa  = REF,
        vcf = "VCF/ss14/{sample}_HaploCaller.ss14.filt.snp.vcf.gz"
    output:
        "VCF/ss14/{sample}_HaploCaller.ss14.pass.snp.vcf.gz"
    threads: 1
    resources: mem_mb=12000
    params:
        java_opts = "-Xmx8g"
    shell:
        r"""
        conda run -n tpa_gatk_env \
            gatk --java-options "{params.java_opts}" SelectVariants \
                -R {input.fa} -V {input.vcf} -O {output} \
                -select "vc.isNotFiltered()"
        """

rule select_pass_indel:
    input:
        fa  = REF,
        vcf = "VCF/ss14/{sample}_HaploCaller.ss14.filt.indel.vcf.gz"
    output:
        temp("VCF/ss14/{sample}_HaploCaller.ss14.pass.indel.vcf.gz")
    threads: 1
    resources: mem_mb=12000
    params:
        java_opts = "-Xmx8g"
    shell:
        r"""
        conda run -n tpa_gatk_env \
            gatk --java-options "{params.java_opts}" SelectVariants \
                -R {input.fa} -V {input.vcf} -O {output} \
                -select "vc.isNotFiltered()"
        """

rule merge_pass_vcf:
    input:
        fa    = REF,
        snp   = "VCF/ss14/{sample}_HaploCaller.ss14.pass.snp.vcf.gz",
        indel = "VCF/ss14/{sample}_HaploCaller.ss14.pass.indel.vcf.gz"
    output:
        "VCF/ss14/{sample}_HaploCaller.ss14.pass.vcf.gz"
    threads: 1
    resources: mem_mb=12000
    params:
        java_opts = "-Xmx8g"
    shell:
        r"""
        conda run -n tpa_gatk_env \
            gatk --java-options "{params.java_opts}" MergeVcfs \
                -I {input.snp} \
                -I {input.indel} \
                -O {output}
        """

#######################################################
##### 4b. Normalize → Annotate → Extract #############
#######################################################

rule normalize_vcf:
    input:
        fa  = REF,
        vcf = "VCF/ss14/{sample}_HaploCaller.ss14.pass.vcf.gz"
    output:
        vcf = "VCF/ss14/{sample}_HaploCaller.ss14.pass.norm.vcf.gz",
        tbi = "VCF/ss14/{sample}_HaploCaller.ss14.pass.norm.vcf.gz.tbi"
    threads: 1
    resources: mem_mb=8000
    shell:
        r"""
        conda run -n tpa_main_env \
            bash -c 'bcftools norm -f {input.fa} -m -any {input.vcf} | bgzip -c > {output.vcf}'
        conda run -n tpa_main_env \
            bcftools index -t {output.vcf}
        """

rule annotate_vcf:
    input:
        vcf = "VCF/ss14/{sample}_HaploCaller.ss14.pass.norm.vcf.gz",
        tbi = "VCF/ss14/{sample}_HaploCaller.ss14.pass.norm.vcf.gz.tbi"
    output:
        vcf  = "VCF/ss14/{sample}_HaploCaller.ss14.pass.ann.vcf.gz",
        tbi  = "VCF/ss14/{sample}_HaploCaller.ss14.pass.ann.vcf.gz.tbi",
        html = "VCF/ss14/{sample}_snpEff_summary.html",
        csv  = "VCF/ss14/{sample}_snpEff_summary.csv"
    threads: 2
    resources: mem_mb=12000
    params:
        db       = "TPA_SS14",
        java_mem = "-Xmx8g"
    shell:
        r"""
        conda run -n tpa_main_env \
            bash -c 'snpEff {params.java_mem} \
                -v {params.db} \
                -stats {output.html} \
                -csvStats {output.csv} \
                {input.vcf} | bgzip -c > {output.vcf}'
        conda run -n tpa_main_env \
            bcftools index -t {output.vcf}
        """

rule extract_ann_fields:
    input:
        vcf = "VCF/ss14/{sample}_HaploCaller.ss14.pass.ann.vcf.gz"
    output:
        tsv = "VCF/ss14/{sample}.ann.tsv"
    threads: 1
    resources: mem_mb=8000
    shell:
        r"""
        conda run -n tpa_main_env \
            bash -c 'SnpSift extractFields -s "," -e "." \
                {input.vcf} \
                CHROM POS REF ALT \
                "ANN[*].GENE" \
                "ANN[*].GENEID" \
                "ANN[*].EFFECT" \
                "ANN[*].IMPACT" \
                "ANN[*].HGVS_C" \
                "ANN[*].HGVS_P" \
                "ANN[*].RANK" \
                GEN[0].DP GEN[0].GT \
            > {output.tsv}'
        """

#######################################################
##### 5. Consensus FASTA ##############################
#######################################################

rule make_depth_mask:
    input:
        bam = "bam/ss14/postfilter/{sample}.ss14.final.bam",
        bai = "bam/ss14/postfilter/{sample}.ss14.final.bam.bai"
    output:
        depth = temp("consensus/ss14/{sample}.ss14.depth.txt"),
        mask  = "consensus/ss14/{sample}.ss14.lowcov_mask.bed"
    threads: 2
    resources: mem_mb=8000
    shell:
        r"""
        mkdir -p consensus/ss14
        conda run -n tpa_main_env \
            bash -c 'samtools depth -a -q 0 -Q 0 {input.bam} > {output.depth} && \
            awk "\$3 < 3 {{print \$1\"\t\"(\$2-1)\"\t\"\$2}}" {output.depth} | \
            bedtools merge -i - > {output.mask}'
        """

rule make_consensus:
    input:
        fa   = REF,
        vcf  = "VCF/ss14/{sample}_HaploCaller.ss14.pass.vcf.gz",
        mask = "consensus/ss14/{sample}.ss14.lowcov_mask.bed"
    output:
        "consensus/ss14/{sample}.ss14.consensus.fasta"
    threads: 1
    resources: mem_mb=8000
    shell:
        r"""
        conda run -n tpa_main_env bcftools index -f {input.vcf}
        conda run -n tpa_main_env \
            bash -c 'bcftools consensus \
                -f {input.fa} \
                --mask {input.mask} \
                --mask-with N \
                --sample {wildcards.sample} \
                {input.vcf} | \
            sed "s/^>.*/>{wildcards.sample}/" > {output}'
        """

rule make_snp_consensus:
    input:
        fa   = REF,
        vcf  = "VCF/ss14/{sample}_HaploCaller.ss14.pass.snp.vcf.gz",
        mask = "consensus/ss14/{sample}.ss14.lowcov_mask.bed"
    output:
        "consensus/ss14/{sample}.ss14.snp_only.consensus.fasta"
    threads: 1
    resources: mem_mb=8000
    shell:
        r"""
        conda run -n tpa_main_env bcftools index -f {input.vcf}
        conda run -n tpa_main_env \
            bash -c 'bcftools consensus \
                -f {input.fa} \
                --mask {input.mask} \
                --mask-with N \
                --sample {wildcards.sample} \
                {input.vcf} | \
            sed "s/^>.*/>{wildcards.sample}/" > {output}'
        """
