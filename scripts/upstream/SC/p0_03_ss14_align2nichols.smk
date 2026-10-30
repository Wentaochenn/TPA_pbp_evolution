import pandas as pd
MyList = pd.read_csv("list/tpa.ss14.sampleList.txt", header=None)
SAMPLES = list(MyList[0])


rule all:
    input:
        expand("VCF/ss14_align2_nichols/{sample}_HaploCaller.ss14_a2n.pass.vcf.gz",          sample=SAMPLES),
        expand("VCF/ss14_align2_nichols/{sample}_HaploCaller.ss14_a2n.pass.norm.vcf.gz",     sample=SAMPLES),
        expand("VCF/ss14_align2_nichols/{sample}_HaploCaller.ss14_a2n.pass.ann.vcf.gz",      sample=SAMPLES),
        expand("VCF/ss14_align2_nichols/{sample}.ann.tsv",                               sample=SAMPLES),
        expand("consensus/ss14_align2_nichols/{sample}.ss14_a2n.consensus.fasta",            sample=SAMPLES),
        expand("consensus/ss14_align2_nichols/{sample}.ss14_a2n.snp_only.consensus.fasta",   sample=SAMPLES)

# 最终输出说明:
# {sample}_HaploCaller.ss14_a2n.pass.vcf.gz
#   └── SNP:   QD<2, MQ<40, FS>60, SOR>4, MQRankSum<-12.5, ReadPosRankSum<-8, DP<3
#   └── INDEL: QD<2, FS>200, ReadPosRankSum<-20, DP<3
# {sample}.ss14_a2n.consensus.fasta
#   └── 包含 SNP+INDEL，DP<3 区域用 N 遮蔽


#######################################################
##### T. Pallidum subsp SS14 reference ################
#######################################################

REF             = "refs/TPA_Nichols.fasta"
HOST_REF_HUMAN = os.environ.get("PBP_HOST_REF_HUMAN", "refs/host_human_bbmap")
HOST_REF_RABBIT = os.environ.get("PBP_HOST_REF_RABBIT", "refs/host_rabbit_bbmap")
TRIMMOMATIC_ADAPTERS = os.environ.get("PBP_TRIMMOMATIC_ADAPTERS", "refs/TruSeq3-PE.fa")
PICARD_TMP      = "tmp"

#######################################################
##### 1. Pre-alignment QC #############################
#######################################################
'''
rule adapter_trim:
    input:
        fq_1 = "raw_fastq/ss14/{sample}_1.fastq.gz",
        fq_2 = "raw_fastq/ss14/{sample}_2.fastq.gz"
    output:
        fq_1          = "processed_fastq/ss14_align2_nichols/adapter_trimmed/{sample}_R1.paired_trim.fq.gz",
        fq_2          = "processed_fastq/ss14_align2_nichols/adapter_trimmed/{sample}_R2.paired_trim.fq.gz",
        fq_1_unpaired = "processed_fastq/ss14_align2_nichols/adapter_trimmed/{sample}_R1.unpaired.fq.gz",
        fq_2_unpaired = "processed_fastq/ss14_align2_nichols/adapter_trimmed/{sample}_R2.unpaired.fq.gz",
        log_file      = "processed_fastq/ss14_align2_nichols/adapter_trimmed/{sample}_trim_log.txt"
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
        fq_1 = "processed_fastq/ss14_align2_nichols/adapter_trimmed/{sample}_R1.paired_trim.fq.gz",
        fq_2 = "processed_fastq/ss14_align2_nichols/adapter_trimmed/{sample}_R2.paired_trim.fq.gz"
    output:
        fq_clean  = temp("processed_fastq/ss14_align2_nichols/host_filtered/{sample}.rm_human_R1.fq.gz"),
        fq_clean2 = temp("processed_fastq/ss14_align2_nichols/host_filtered/{sample}.rm_human_R2.fq.gz"),
        fq_cnt    = "processed_fastq/ss14_align2_nichols/host_filtered/{sample}.human_contamination.fq.gz"
    threads: 16
    resources: mem_mb=80000, big_mem_slots=1
    params:
        ref = HOST_REF_HUMAN
    shell:
        r"""
        conda run -n tpa_main_env \
            bbmap.sh build=1 threads={threads} -Xmx64g \
                minid=0.95 maxindel=3 bandwidthratio=0.16 bandwidth=12 \
                quickmatch fast minhits=2 pigz unpigz \
                path={params.ref} \
                in={input.fq_1} in2={input.fq_2} \
                outm={output.fq_cnt} \
                outu1={output.fq_clean} outu2={output.fq_clean2}
        """

rule fq_rm_host_rabbit:
    input:
        fq_1 = "processed_fastq/ss14_align2_nichols/host_filtered/{sample}.rm_human_R1.fq.gz",
        fq_2 = "processed_fastq/ss14_align2_nichols/host_filtered/{sample}.rm_human_R2.fq.gz"
    output:
        fq_clean  = temp("processed_fastq/ss14_align2_nichols/host_filtered/{sample}.rm_rabbit_R1.fq.gz"),
        fq_clean2 = temp("processed_fastq/ss14_align2_nichols/host_filtered/{sample}.rm_rabbit_R2.fq.gz"),
        fq_cnt    = "processed_fastq/ss14_align2_nichols/host_filtered/{sample}.rabbit_contamination.fq.gz"
    threads: 16
    resources: mem_mb=80000, big_mem_slots=1
    params:
        ref = HOST_REF_RABBIT
    shell:
        r"""
        conda run -n tpa_main_env \
            bbmap.sh build=1 threads={threads} -Xmx64g \
                minid=0.95 maxindel=3 bandwidthratio=0.16 bandwidth=12 \
                quickmatch fast minhits=2 pigz unpigz \
                path={params.ref} \
                in={input.fq_1} in2={input.fq_2} \
                outm={output.fq_cnt} \
                outu1={output.fq_clean} outu2={output.fq_clean2}
        """

rule fq_rm_broken:
    input:
        fq_1 = "processed_fastq/ss14_align2_nichols/host_filtered/{sample}.rm_rabbit_R1.fq.gz",
        fq_2 = "processed_fastq/ss14_align2_nichols/host_filtered/{sample}.rm_rabbit_R2.fq.gz"
    output:
        fq_1 = "processed_fastq/ss14_align2_nichols/host_filtered/{sample}_R1.host_filtered.fq.gz",
        fq_2 = "processed_fastq/ss14_align2_nichols/host_filtered/{sample}_R2.host_filtered.fq.gz"
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
        fq_1 = "processed_fastq/ss14_align2_nichols/host_filtered/{sample}_R1.host_filtered.fq.gz",
        fq_2 = "processed_fastq/ss14_align2_nichols/host_filtered/{sample}_R2.host_filtered.fq.gz"
    output:
        temp("processed_fastq/ss14_align2_nichols/strainseeker/{sample}.sub")
    threads: 2
    resources: mem_mb=8000
    shell:
        r"""
        conda run -n tpa_main_env \
            bash -c 'seqtk sample -s100 {input.fq_1} 100000 > {output}'
        """

rule strainseeker:
    input:
        "processed_fastq/ss14_align2_nichols/strainseeker/{sample}.sub"
    output:
        "processed_fastq/ss14_align2_nichols/strainseeker/{sample}.seeker.txt"
    threads: 2
    resources: mem_mb=16000
    shell:
        r"""
        mkdir -p tmp/strainseeker/{wildcards.sample}
        cd tmp/strainseeker/{wildcards.sample} && \
        conda run -n tpa_main_env \
            seeker.pl \
                -i {workflow.basedir}/processed_fastq/ss14_align2_nichols/strainseeker/{wildcards.sample}.sub \
                -d {STRAINSEEKER_DB:q} \
                -o {workflow.basedir}/processed_fastq/ss14_align2_nichols/strainseeker/{wildcards.sample}.seeker.txt
        """
'''
#######################################################
##### 2. Alignment ####################################
#######################################################

rule bwa_align:
    input:
        fa           = REF,
        fq_1         = "processed_fastq/ss14/host_filtered/{sample}_R1.host_filtered.fq.gz",
        fq_2         = "processed_fastq/ss14/host_filtered/{sample}_R2.host_filtered.fq.gz"
    output:
        temp("bam/ss14_align2_nichols/aln/{sample}.ss14_a2n.raw.bam")
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
        "bam/ss14_align2_nichols/aln/{sample}.ss14_a2n.raw.bam"
    output:
        temp("bam/ss14_align2_nichols/aln/{sample}.ss14_a2n.fixmate.bam")
    threads: 4
    resources: mem_mb=8000
    shell:
        r"""
        conda run -n tpa_main_env \
            samtools fixmate -@ {threads} -O bam {input} {output}
        """

rule bam_sort:
    input:
        "bam/ss14_align2_nichols/aln/{sample}.ss14_a2n.fixmate.bam"
    output:
        temp("bam/ss14_align2_nichols/aln/{sample}.ss14_a2n.sorted.bam")
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
        "bam/ss14_align2_nichols/aln/{sample}.ss14_a2n.sorted.bam"
    output:
        temp("bam/ss14_align2_nichols/aln/{sample}.ss14_a2n.sorted.bam.bai")
    threads: 4
    resources: mem_mb=8000
    shell:
        r"""
        conda run -n tpa_main_env \
            samtools index -@ {threads} {input}
        """

rule flagstat:
    input:
        bam = "bam/ss14_align2_nichols/aln/{sample}.ss14_a2n.sorted.bam",
        bai = "bam/ss14_align2_nichols/aln/{sample}.ss14_a2n.sorted.bam.bai"
    output:
        temp("bam/ss14_align2_nichols/aln/{sample}.ss14_a2n.flagstat")
    threads: 4
    resources: mem_mb=8000
    shell:
        r"""
        conda run -n tpa_main_env \
            bash -c 'samtools flagstat -@ {threads} {input.bam} > {output}'
        """

rule filter_flags:
    input:
        bam      = "bam/ss14_align2_nichols/aln/{sample}.ss14_a2n.sorted.bam",
        bai      = "bam/ss14_align2_nichols/aln/{sample}.ss14_a2n.sorted.bam.bai",
        flagstat = "bam/ss14_align2_nichols/aln/{sample}.ss14_a2n.flagstat"
    output:
        temp("bam/ss14_align2_nichols/aln/{sample}.ss14_a2n.flag_filt.bam")
    threads: 4
    resources: mem_mb=8000
    shell:
        r"""
        conda run -n tpa_main_env \
            samtools view -@ {threads} -h -F 4 -F 256 -f 2 -b -o {output} {input.bam}
        """

rule index_flag_filter:
    input:
        "bam/ss14_align2_nichols/aln/{sample}.ss14_a2n.flag_filt.bam"
    output:
        temp("bam/ss14_align2_nichols/aln/{sample}.ss14_a2n.flag_filt.bam.bai")
    threads: 4
    resources: mem_mb=8000
    shell:
        r"""
        conda run -n tpa_main_env \
            samtools index -@ {threads} {input}
        """

rule mark_dups:
    input:
        bam = "bam/ss14_align2_nichols/aln/{sample}.ss14_a2n.flag_filt.bam",
        bai = "bam/ss14_align2_nichols/aln/{sample}.ss14_a2n.flag_filt.bam.bai"
    output:
        bam     = "bam/ss14_align2_nichols/aln/{sample}.ss14_a2n.dedup.bam",
        metrics = "bam/ss14_align2_nichols/aln/{sample}.dedup.metrics"
    threads: 4
    resources: mem_mb=16000
    log: "bam/ss14_align2_nichols/aln/{sample}.ss14_a2n.markdup.log.txt"
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
        "bam/ss14_align2_nichols/aln/{sample}.ss14_a2n.dedup.bam"
    output:
        "bam/ss14_align2_nichols/aln/{sample}.ss14_a2n.dedup.bam.bai"
    threads: 4
    resources: mem_mb=8000
    shell:
        r"""
        conda run -n tpa_main_env \
            samtools index -@ {threads} {input}
        """

rule bam_fixmate_b:
    input:
        bam = "bam/ss14_align2_nichols/aln/{sample}.ss14_a2n.dedup.bam",
        bai = "bam/ss14_align2_nichols/aln/{sample}.ss14_a2n.dedup.bam.bai"
    output:
        bam = temp("bam/ss14_align2_nichols/aln/{sample}.ss14_a2n.matefixed.bam"),
        bai = temp("bam/ss14_align2_nichols/aln/{sample}.ss14_a2n.matefixed.bai")
    threads: 4
    resources: mem_mb=16000
    log: "bam/ss14_align2_nichols/aln/{sample}.ss14_a2n.matefixed.log.txt"
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
        bam = "bam/ss14_align2_nichols/aln/{sample}.ss14_a2n.matefixed.bam"
    output:
        bam    = temp("bam/ss14_align2_nichols/postfilter/{sample}.ss14_a2n.rm_mismatch.bam"),
        failed = "bam/ss14_align2_nichols/postfilter/{sample}.ss14_a2n.mismatch_reads.txt"
    threads: 2
    resources: mem_mb=16000
    shell:
        r"""
        conda run -n tpa_ngsutils_env \
            bamutils filter {input.bam} {output.bam} -failed {output.failed} -mismatch 5
        """

rule short_alignment:
    input:
        "bam/ss14_align2_nichols/postfilter/{sample}.ss14_a2n.rm_mismatch.bam"
    output:
        temp("bam/ss14_align2_nichols/postfilter/{sample}.ss14_a2n.minlen.bam")
    threads: 2
    resources: mem_mb=16000
    shell:
        r"""
        conda run -n tpa_ngsutils_env \
            bamutils filter {input} {output} -minlen 35
        """

rule soft_clip_a:
    input:
        "bam/ss14_align2_nichols/postfilter/{sample}.ss14_a2n.minlen.bam"
    output:
        temp("bam/ss14_align2_nichols/postfilter/{sample}.ss14_a2n.softclips.tmp.bam")
    threads: 2
    resources: mem_mb=32000
    shell:
        r"""
        conda run -n tpa_ngsutils_env \
            bamutils removeclipping {input} {output}
        """

rule soft_clip_b:
    input:
        "bam/ss14_align2_nichols/postfilter/{sample}.ss14_a2n.softclips.tmp.bam"
    output:
        "bam/ss14_align2_nichols/postfilter/{sample}.read_names_to_remove_highSoftClip.txt"
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
        bam      = "bam/ss14_align2_nichols/postfilter/{sample}.ss14_a2n.minlen.bam",
        read_out = "bam/ss14_align2_nichols/postfilter/{sample}.read_names_to_remove_highSoftClip.txt"
    output:
        bam_sc       = temp("bam/ss14_align2_nichols/postfilter/{sample}.ss14_a2n.softclip_filtOut_reads.bam"),
        bam_filtered = temp("bam/ss14_align2_nichols/postfilter/{sample}.ss14_a2n.rm_softclips.bam")
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
        "bam/ss14_align2_nichols/postfilter/{sample}.ss14_a2n.rm_softclips.bam"
    output:
        bam               = temp("bam/ss14_align2_nichols/postfilter/{sample}.ss14_a2n.rm_chim.bam"),
        bam_filtOut_reads = temp("bam/ss14_align2_nichols/postfilter/{sample}.ss14_a2n.chim_reads.bam")
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
        "bam/ss14_align2_nichols/postfilter/{sample}.ss14_a2n.rm_chim.bam"
    output:
        bam = temp("bam/ss14_align2_nichols/postfilter/{sample}.ss14_a2n.rm_hardclips.bam")
    threads: 4
    resources: mem_mb=8000
    shell:
        r"""
        conda run -n tpa_main_env \
            bash -c 'samtools view -@ 2 -h {input} | awk '"'"'$6 !~ /H/'"'"' | samtools view -bh -o {output.bam} -'
        """

rule MAPQ_filt:
    input:
        "bam/ss14_align2_nichols/postfilter/{sample}.ss14_a2n.rm_hardclips.bam"
    output:
        bam = temp("bam/ss14_align2_nichols/postfilter/{sample}.ss14_a2n.MAPQ40.bam")
    threads: 4
    resources: mem_mb=8000
    shell:
        r"""
        conda run -n tpa_main_env \
            bash -c 'samtools view -@ 2 -h -bq 40 {input} > {output.bam}'
        """

rule rm_singletone:
    input:
        "bam/ss14_align2_nichols/postfilter/{sample}.ss14_a2n.MAPQ40.bam"
    output:
        "bam/ss14_align2_nichols/postfilter/{sample}.ss14_a2n.final.bam"
    threads: 4
    resources: mem_mb=8000
    shell:
        r"""
        conda run -n tpa_main_env \
            samtools view -@ 2 -h -F 8 -b -o {output} {input}
        """

rule index_final_bam:
    input:
        "bam/ss14_align2_nichols/postfilter/{sample}.ss14_a2n.final.bam"
    output:
        "bam/ss14_align2_nichols/postfilter/{sample}.ss14_a2n.final.bam.bai"
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
        bam = "bam/ss14_align2_nichols/postfilter/{sample}.ss14_a2n.final.bam",
        bai = "bam/ss14_align2_nichols/postfilter/{sample}.ss14_a2n.final.bam.bai"
    output:
        "VCF/ss14_align2_nichols/{sample}_HaploCaller.ss14_a2n.raw.g.vcf.gz"
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
        gvcf = "VCF/ss14_align2_nichols/{sample}_HaploCaller.ss14_a2n.raw.g.vcf.gz"
    output:
        "VCF/ss14_align2_nichols/{sample}_HaploCaller.ss14_a2n.raw.vcf.gz"
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
        vcf = "VCF/ss14_align2_nichols/{sample}_HaploCaller.ss14_a2n.raw.vcf.gz"
    output:
        "VCF/ss14_align2_nichols/{sample}_HaploCaller.ss14_a2n.raw.snp.vcf.gz"
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
        vcf = "VCF/ss14_align2_nichols/{sample}_HaploCaller.ss14_a2n.raw.vcf.gz"
    output:
        temp("VCF/ss14_align2_nichols/{sample}_HaploCaller.ss14_a2n.raw.indel.vcf.gz")
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
        vcf = "VCF/ss14_align2_nichols/{sample}_HaploCaller.ss14_a2n.raw.snp.vcf.gz"
    output:
        temp("VCF/ss14_align2_nichols/{sample}_HaploCaller.ss14_a2n.filt.snp.vcf.gz")
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
        vcf = "VCF/ss14_align2_nichols/{sample}_HaploCaller.ss14_a2n.raw.indel.vcf.gz"
    output:
        temp("VCF/ss14_align2_nichols/{sample}_HaploCaller.ss14_a2n.filt.indel.vcf.gz")
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
        vcf = "VCF/ss14_align2_nichols/{sample}_HaploCaller.ss14_a2n.filt.snp.vcf.gz"
    output:
        "VCF/ss14_align2_nichols/{sample}_HaploCaller.ss14_a2n.pass.snp.vcf.gz"
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
        vcf = "VCF/ss14_align2_nichols/{sample}_HaploCaller.ss14_a2n.filt.indel.vcf.gz"
    output:
        temp("VCF/ss14_align2_nichols/{sample}_HaploCaller.ss14_a2n.pass.indel.vcf.gz")
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
        snp   = "VCF/ss14_align2_nichols/{sample}_HaploCaller.ss14_a2n.pass.snp.vcf.gz",
        indel = "VCF/ss14_align2_nichols/{sample}_HaploCaller.ss14_a2n.pass.indel.vcf.gz"
    output:
        "VCF/ss14_align2_nichols/{sample}_HaploCaller.ss14_a2n.pass.vcf.gz"
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
        vcf = "VCF/ss14_align2_nichols/{sample}_HaploCaller.ss14_a2n.pass.vcf.gz"
    output:
        vcf = "VCF/ss14_align2_nichols/{sample}_HaploCaller.ss14_a2n.pass.norm.vcf.gz",
        tbi = "VCF/ss14_align2_nichols/{sample}_HaploCaller.ss14_a2n.pass.norm.vcf.gz.tbi"
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
        vcf = "VCF/ss14_align2_nichols/{sample}_HaploCaller.ss14_a2n.pass.norm.vcf.gz",
        tbi = "VCF/ss14_align2_nichols/{sample}_HaploCaller.ss14_a2n.pass.norm.vcf.gz.tbi"
    output:
        vcf  = "VCF/ss14_align2_nichols/{sample}_HaploCaller.ss14_a2n.pass.ann.vcf.gz",
        tbi  = "VCF/ss14_align2_nichols/{sample}_HaploCaller.ss14_a2n.pass.ann.vcf.gz.tbi",
        html = "VCF/ss14_align2_nichols/{sample}_snpEff_summary.html",
        csv  = "VCF/ss14_align2_nichols/{sample}_snpEff_summary.csv"
    threads: 2
    resources: mem_mb=12000
    params:
        db       = "TPA_Nichols",
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
        vcf = "VCF/ss14_align2_nichols/{sample}_HaploCaller.ss14_a2n.pass.ann.vcf.gz"
    output:
        tsv = "VCF/ss14_align2_nichols/{sample}.ann.tsv"
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
        bam = "bam/ss14_align2_nichols/postfilter/{sample}.ss14_a2n.final.bam",
        bai = "bam/ss14_align2_nichols/postfilter/{sample}.ss14_a2n.final.bam.bai"
    output:
        depth = temp("consensus/ss14_align2_nichols/{sample}.ss14_a2n.depth.txt"),
        mask  = "consensus/ss14_align2_nichols/{sample}.ss14_a2n.lowcov_mask.bed"
    threads: 2
    resources: mem_mb=8000
    shell:
        r"""
        mkdir -p consensus/ss14_align2_nichols
        conda run -n tpa_main_env \
            bash -c 'samtools depth -a -q 0 -Q 0 {input.bam} > {output.depth} && \
            awk "\$3 < 3 {{print \$1\"\t\"(\$2-1)\"\t\"\$2}}" {output.depth} | \
            bedtools merge -i - > {output.mask}'
        """

rule make_consensus:
    input:
        fa   = REF,
        vcf  = "VCF/ss14_align2_nichols/{sample}_HaploCaller.ss14_a2n.pass.vcf.gz",
        mask = "consensus/ss14_align2_nichols/{sample}.ss14_a2n.lowcov_mask.bed"
    output:
        "consensus/ss14_align2_nichols/{sample}.ss14_a2n.consensus.fasta"
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
        vcf  = "VCF/ss14_align2_nichols/{sample}_HaploCaller.ss14_a2n.pass.snp.vcf.gz",
        mask = "consensus/ss14_align2_nichols/{sample}.ss14_a2n.lowcov_mask.bed"
    output:
        "consensus/ss14_align2_nichols/{sample}.ss14_a2n.snp_only.consensus.fasta"
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
