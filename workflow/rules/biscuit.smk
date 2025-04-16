###-----------------------------------------------------------------------------------------------------------------###
# Notes on where variables are defined, if not defined in biscuit.smk
#
# output_directory - workflow/Snakefile
# config           - workflow/Snakefile
# BISCUIT_INDEX_FORMATS - workflow/Snakefile
#
###-----------------------------------------------------------------------------------------------------------------###

if config['build_ref_with_methylation_controls']:
    rule build_ref_with_methylation_controls:
        input:
            config['ref']['fasta'],
        output:
            ref = 'merged_reference/merged.fa.gz',
            refdir = directory('merged_reference/'),
            idxfil = expand('merged_reference/merged.fa.gz.{ext}', ext=BISCUIT_INDEX_FORMATS),
        log:
            f'{output_directory}/logs/build_merged_reference_index.log',
        benchmark:
            f'{output_directory}/benchmarks/build_merged_reference_index.txt',
        threads: 2
        resources:
            mem_gb=config['hpc_parameters']['small_memory_gb'],
            time = config['runtime']['medium'],
        conda:
            '../envs/biscuit.yaml'
        envmodules:
            config['envmodules']['biscuit'],
            config['envmodules']['samtools'],
        shell:
            """
            mkdir -p {output.refdir}

            if (file {input} | grep -q "extra field"); then
                cat <(bgzip -d {input}) <(zcat bin/puc19.fa.gz) <(zcat bin/lambda.fa.gz) | bgzip > {output.ref}
            elif (file {input} | grep -q "gzip compressed data, was"); then
                cat <(zcat {input}) <(zcat bin/puc19.fa.gz) <(zcat bin/lambda.fa.gz) | bgzip > {output.ref}
            else
                cat {input} <(zcat bin/puc19.fa.gz) <(zcat bin/lambda.fa.gz) | bgzip > {output.ref}
            fi

            biscuit index {output.ref} 2> {log}
            samtools faidx {output.ref} 2>> {log}
            """

def get_biscuit_reference(wildcards):
    if config['build_ref_with_methylation_controls']: # Currently not set up to be generated
        return 'merged_reference/merged.fa.gz'
    else:
        return config['ref']['fasta']

def get_biscuit_index(wildcards):
    if config['build_ref_with_methylation_controls']: # Currently not set up to be generated
        return 'merged_reference/merged.fa.gz'
    else:
        return config['ref']['index']

def get_rename_fastq_output_R1(wildcards):
    cp_output = checkpoints.rename_fastq_files.get().output.symlink_dir

    if config['trim_galore']['trim_before_biscuit']:
        return output_directory + '/analysis/trim_reads/' + wildcards.sample + '-R1_val_1_merged.fq.gz'
    else:
        IDX, = glob_wildcards(cp_output + '/' + wildcards.sample + '-{id}-R1.fastq.gz')
        files = list(expand(cp_output + '/' + wildcards.sample + '-{idx}-R1.fastq.gz', idx = IDX))
        files.sort()
        return files
        
def get_rename_fastq_output_R2(wildcards):
    cp_output = checkpoints.rename_fastq_files.get().output.symlink_dir

    if config['trim_galore']['trim_before_biscuit']:
        return output_directory + '/analysis/trim_reads/' + wildcards.sample + '-R2_val_2_merged.fq.gz'
    else:
        IDX, = glob_wildcards(cp_output + '/' + wildcards.sample + '-{id}-R2.fastq.gz')
        files = list(expand(cp_output + '/' + wildcards.sample + '-{idx}-R2.fastq.gz', idx = IDX))
        files.sort()
        return files

rule biscuit_sifter:
    input:
        reference = get_biscuit_reference,
        index = get_biscuit_index,
        R1 = get_rename_fastq_output_R1,
        R2 = get_rename_fastq_output_R2,
    output:
        bam = f'{output_directory}/analysis/align/{{sample}}.sorted.markdup.bam',
        bai = f'{output_directory}/analysis/align/{{sample}}.sorted.markdup.bam.bai',
        dup = f'{output_directory}/analysis/align/{{sample}}.dupsifter.stat',
    params:
        args_list = config['biscuit']['args_align'],
        LB = config['sam_header']['LB'],
        PL = config['sam_header']['PL'],
        PU = config['sam_header']['PU'],
        SM = '{sample}', # also used for ID
        al_threads = config['hpc_parameters']['biscuit_sifter_threads'],
        st_threads = config['hpc_parameters']['samtools_index_threads'],
    log:
        biscuit = f'{output_directory}/logs/biscuit/biscuit_sifter.{{sample}}.log',
        dupsifter = f'{output_directory}/logs/biscuit/dupsifter.{{sample}}.log',
        samtools_sort = f'{output_directory}/logs/biscuit/samtools_sort.{{sample}}.log',
        samtools_index = f'{output_directory}/logs/biscuit/samtools_index.{{sample}}.log',
    benchmark:
        f'{output_directory}/benchmarks/biscuit_sifter/{{sample}}.txt'
    threads: config['hpc_parameters']['biscuit_sifter_threads'] + config['hpc_parameters']['samtools_index_threads']
    resources:
        mem_gb = config['hpc_parameters']['max_memory_gb'],
        time = config['runtime']['long'],
    conda:
        '../envs/biscuit.yaml'
    envmodules:
        config['envmodules']['biscuit'],
        config['envmodules']['dupsifter'],
        config['envmodules']['samtools'],
        config['envmodules']['htslib'],
    shell:
        """
        # biscuitSifter pipeline
        biscuit align \
            {params.args_list} \
            -@ {params.al_threads} \
            -R '@RG\tLB:{params.LB}\tID:{params.SM}\tPL:{params.PL}\tPU:{params.PU}\tSM:{params.SM}' \
            {input.index} \
            <(zcat {input.R1}) \
            <(zcat {input.R2}) 2> {log.biscuit} | \
        dupsifter --stats-output {output.dup} {input.reference} 2> {log.dupsifter} | \
        samtools sort -@ {params.st_threads} -m 5G -o {output.bam} -O BAM - 2> {log.samtools_sort}

        samtools index -@ {params.st_threads} {output.bam} 2> {log.samtools_index}
        """

rule biscuit_pileup:
    input:
        ref = get_biscuit_reference,
        bam = f'{output_directory}/analysis/align/{{sample}}.sorted.markdup.bam',
    output:
        vcf_gz = f'{output_directory}/analysis/pileup/{{sample}}.vcf.gz',
        vcf_tabix = f'{output_directory}/analysis/pileup/{{sample}}.vcf.gz.tbi',
        meth = f'{output_directory}/analysis/pileup/{{sample}}.vcf_meth_average.tsv',
        bed_gz = f'{output_directory}/analysis/pileup/{{sample}}.bed.gz',
        bed_tbi = f'{output_directory}/analysis/pileup/{{sample}}.bed.gz.tbi',
    params:
        args_pileup = config['biscuit']['args_pileup'],
        args_vcf2bed_cg = config['biscuit']['args_vcf2bed_cg'],
        vcf = f'{output_directory}/analysis/pileup/{{sample}}.vcf',
        bed = f'{output_directory}/analysis/pileup/{{sample}}.bed',
    log:
        pileup = f'{output_directory}/logs/biscuit_pileup/{{sample}}.pileup.log',
        vcf_gz = f'{output_directory}/logs/biscuit_pileup/{{sample}}.vcf_gz.log',
        vcf_tbi = f'{output_directory}/logs/biscuit_pileup/{{sample}}.vcf_tbi.log',
        vcf2bed = f'{output_directory}/logs/biscuit_pileup/{{sample}}.vcf2bed.log',
        bed_gz = f'{output_directory}/logs/biscuit_pileup/{{sample}}.bed_gz.log',
        bed_tbi = f'{output_directory}/logs/biscuit_pileup/{{sample}}.bed_tabix.log',
    benchmark:
        f'{output_directory}/benchmarks/biscuit_pileup/{{sample}}.txt',
    threads: config['hpc_parameters']['pileup_threads']
    resources:
        mem_gb = config['hpc_parameters']['intermediate_memory_gb'],
        time = config['runtime']['medium'],
    wildcard_constraints:
        sample = '.*[^(_mergecg)]',
    conda:
        '../envs/biscuit.yaml'
    envmodules:
        config['envmodules']['biscuit'],
        config['envmodules']['htslib'],
    shell:
        """
        biscuit pileup {params.args_pileup} -@ {threads} -o {params.vcf} {input.ref} {input.bam} 2> {log.pileup}

        bgzip -@ {threads} {params.vcf} 2> {log.vcf_gz}
        tabix -p vcf {output.vcf_gz} 2> {log.vcf_tbi}

        biscuit vcf2bed {params.args_vcf2bed_cg} -t cg {output.vcf_gz} 1> {params.bed} 2> {log.vcf2bed}
        bgzip -@ {threads} {params.bed} 2> {log.bed_gz}
        tabix -p bed {output.bed_gz} 2> {log.bed_tbi}
        """

rule biscuit_mergecg:
    input:
        ref = get_biscuit_reference,
        bed = f'{output_directory}/analysis/pileup/{{sample}}.bed.gz',
    output:
        mergecg_gz = f'{output_directory}/analysis/pileup/{{sample}}_mergecg.bed.gz',
        mergecg_tbi = f'{output_directory}/analysis/pileup/{{sample}}_mergecg.bed.gz.tbi',
    params:
        mergecg = f'{output_directory}/analysis/pileup/{{sample}}_mergecg.bed',
        args_mergecg = config['biscuit']['args_mergecg'],
    log:
        mergecg = f'{output_directory}/logs/biscuit_pileup/mergecg.{{sample}}.log',
        mergecg_gz = f'{output_directory}/logs/biscuit_pileup/mergecg_gz.{{sample}}.log',
        mergecg_tbi = f'{output_directory}/logs/biscuit_pileup/mergecg_tabix.{{sample}}.log',
    benchmark:
        f'{output_directory}/benchmarks/biscuit_mergecg/{{sample}}.txt',
    threads: 8
    resources:
        mem_gb = config['hpc_parameters']['intermediate_memory_gb'],
        time = config['runtime']['medium'],
    wildcard_constraints:
        sample = '.*[^(_mergecg)]'
    conda:
        '../envs/biscuit.yaml'
    envmodules:
        config['envmodules']['biscuit'],
        config['envmodules']['htslib'],
    shell:
        """
        biscuit mergecg {params.args_mergecg} {input.ref} {input.bed} 1> {params.mergecg} 2> {log.mergecg}

        bgzip -@ {threads} {params.mergecg} 2> {log.mergecg_gz}
        tabix -p bed {output.mergecg_gz} 2> {log.mergecg_tbi}
        """

rule biscuit_ch_vcf2bed:
    input:
        vcf_gz = f'{output_directory}/analysis/pileup/{{sample}}.vcf.gz',
        vcf_tabix = f'{output_directory}/analysis/pileup/{{sample}}.vcf.gz.tbi',
    output:
        bed = temp(f'{output_directory}/analysis/biscuit_ch_vcf2bed/{{sample}}_ch.bed'),
        bed_gz = f'{output_directory}/analysis/biscuit_ch_vcf2bed/{{sample}}_ch.bed.gz',
        bed_tbi = f'{output_directory}/analysis/biscuit_ch_vcf2bed/{{sample}}_ch.bed.gz.tbi',
    log:
        vcf2bed = f'{output_directory}/logs/biscuit_ch_vcf2bed/{{sample}}_ch.vcf2bed.log',
        bed_gz = f'{output_directory}/logs/biscuit_ch_vcf2bed/{{sample}}_ch.bed_gz.log',
        bed_tbi = f'{output_directory}/logs/biscuit_ch_vcf2bed/{{sample}}_ch.bed_tabix.log',
    benchmark:
        f'{output_directory}/benchmarks/biscuit_ch_vcf2bed/{{sample}}_ch.txt',
    threads: config['hpc_parameters']['pileup_threads']
    resources:
        mem_gb = config['hpc_parameters']['intermediate_memory_gb'],
        time = config['runtime']['medium'],
    conda:
        '../envs/biscuit.yaml'
    envmodules:
        config['envmodules']['biscuit'],
        config['envmodules']['htslib'],
    shell:
        """
        biscuit vcf2bed -t ch -e {input.vcf_gz} 1> {output.bed} 2> {log.vcf2bed}
        bgzip -c -@ {threads} {output.bed} 1> {output.bed_gz} 2> {log.bed_gz}
        tabix -p bed {output.bed_gz} 2> {log.bed_tbi}
        """

rule split_ch_bed:
    input:
        bed_gz = f'{output_directory}/analysis/biscuit_ch_vcf2bed/{{sample}}_ch.bed.gz',
        bed_tbi = f'{output_directory}/analysis/biscuit_ch_vcf2bed/{{sample}}_ch.bed.gz.tbi',
    output:
        bed_gz = f'{output_directory}/analysis/split_ch_bed/{{cx}}/{{sample}}_{{cx}}.bed.gz',
        bed_tbi = f'{output_directory}/analysis/split_ch_bed/{{cx}}/{{sample}}_{{cx}}.bed.gz.tbi',
    log:
        bed_gz = f'{output_directory}/logs/split_ch_bed/{{cx}}/{{sample}}_{{cx}}.bed_gz.log',
        bed_tbi = f'{output_directory}/logs/split_ch_bed/{{cx}}/{{sample}}_{{cx}}.bed_tabix.log',
    benchmark:
        f'{output_directory}/benchmarks/split_ch_bed/{{cx}}/{{sample}}_{{cx}}.txt',
    threads: config['hpc_parameters']['pileup_threads']
    resources:
        mem_gb = config['hpc_parameters']['intermediate_memory_gb'],
        time = config['runtime']['medium'],
    conda:
        '../envs/biscuit.yaml'
    envmodules:
        config['envmodules']['biscuit'],
        config['envmodules']['htslib'],
    shell:
        """
        # Filter for requested CH methylation type.
        # Rearrange columns so that beta is fourth column and coverage is fifth.
        # Delete columns 3-6 to make it compatible with readBiscuit(merged=FALSE)
        zcat {input.bed_gz} | \
        perl -F"\\t" -lane 'next unless $F[5] =~ /^{wildcards.cx}$/; die unless scalar(@F)==9; print join("\\t", @F[0,1,2,7,8])' | \
        bgzip -@ {threads} -c 1> {output.bed_gz} 2> {log.bed_gz}

        tabix -p bed {output.bed_gz} 2> {log.bed_tbi}
        """

def get_bed_for_bigwigs(wildcards):
    if wildcards.cx == "CG":
        return f'{output_directory}/analysis/pileup/{wildcards.sample}.bed.gz'
    if wildcards.cx == "mergecg":
        return f'{output_directory}/analysis/pileup/{wildcards.sample}_mergecg.bed.gz'
    if wildcards.cx in ["CA", "CC", "CT"]:
        return f'{output_directory}/analysis/split_ch_bed/{wildcards.cx}/{wildcards.sample}_{wildcards.cx}.bed.gz'

rule beta_bigwigs:
    input:
        chrom_sizes=f'{output_directory}/analysis/prep_chromsizes_file/chrom_sizes.tsv',
        bed_gz = get_bed_for_bigwigs
    output:
        bed = temp(f'{output_directory}/analysis/beta_bigwigs/{{sample}}_{{cx}}.bed'),
        bigwig=f'{output_directory}/analysis/beta_bigwigs/{{sample}}_{{cx}}.bw',
    params:
        min_cov = config['make_bigwigs']['min_depth']
    log:
        f'{output_directory}/logs/beta_bigwigs/{{sample}}_{{cx}}.log',
    benchmark:
        f'{output_directory}/benchmarks/beta_bigwigs/{{sample}}_{{cx}}.txt'
    threads: 1 
    resources:
        mem_gb = config['hpc_parameters']['intermediate_memory_gb'],
        time = config['runtime']['medium'],
    conda:
        '../envs/ucsc.yaml'
    envmodules:
        config['envmodules']['ucsc'],
        config['envmodules']['htslib'],
    shell:
        """
        # Make BED for bedGraphToBigWig. Filter for minimum coverage.
        zcat {input.bed_gz} | \
        perl -F"\\t" -lane 'print $_ if $F[4] >= {params.min_cov}' | \
        cut -f1-4 > {output.bed}

        bedGraphToBigWig {output.bed} {input.chrom_sizes} {output.bigwig}
        """

rule biscuit_snps:
    input:
        bam = f'{output_directory}/analysis/align/{{sample}}.sorted.markdup.bam',
        vcf_gz = f'{output_directory}/analysis/pileup/{{sample}}.vcf.gz',
    output:
        snp_bed_gz = f'{output_directory}/analysis/snps/{{sample}}.snp.bed.gz',
        snp_bed_gz_tbi = f'{output_directory}/analysis/snps/{{sample}}.snp.bed.gz.tbi',
    params:
        snp_bed = f'{output_directory}/analysis/snps/{{sample}}.snp.bed',
        args_vcf2bed_snp = config['biscuit']['args_vcf2bed_snp'],
    log:
        f'{output_directory}/logs/snps/snps.{{sample}}.log',
    benchmark:
        f'{output_directory}/benchmarks/biscuit_snps/{{sample}}.txt',
    threads: 8
    resources:
        mem_gb = config['hpc_parameters']['intermediate_memory_gb'],
        time = config['runtime']['medium'],
    conda:
        '../envs/biscuit.yaml'
    envmodules:
        config['envmodules']['biscuit'],
        config['envmodules']['htslib'],
    shell:
        """
        biscuit vcf2bed {params.args_vcf2bed_snp} -t snp {input.vcf_gz} > {params.snp_bed} 2> {log}
        bgzip -@ {threads} {params.snp_bed} 2>> {log}
        tabix -p bed {output.snp_bed_gz} 2>> {log}
        """

rule biscuit_epiread:
    input:
        ref = get_biscuit_reference,
        bam = f'{output_directory}/analysis/align/{{sample}}.sorted.markdup.bam',
        snps = f'{output_directory}/analysis/snps/{{sample}}.snp.bed.gz',
        snps_tbi = f'{output_directory}/analysis/snps/{{sample}}.snp.bed.gz.tbi',
    output:
        epibed_gz = f'{output_directory}/analysis/epiread/{{sample}}.epibed.gz',
        epibed_gz_tbi = f'{output_directory}/analysis/epiread/{{sample}}.epibed.gz.tbi',
    params:
        epibed = f'{output_directory}/analysis/epiread/{{sample}}.epibed',
        args_epiread = config['biscuit']['args_epiread'],
    log:
        f'{output_directory}/logs/epiread/epiread.{{sample}}.log',
    benchmark:
        f'{output_directory}/benchmarks/biscuit_epiread/{{sample}}.txt'
    threads: config['hpc_parameters']['pileup_threads']
    resources:
        mem_gb = config['hpc_parameters']['intermediate_memory_gb'],
        time = config['runtime']['medium'],
    conda:
        '../envs/biscuit.yaml'
    envmodules:
        config['envmodules']['biscuit'],
        config['envmodules']['htslib'],
    shell:
        """
        biscuit epiread \
            {params.args_epiread} \
            -@ {threads} \
            -B {input.snps} \
            {input.ref} \
            {input.bam} | \
        sort -k1,1 -k2,2n > {params.epibed} 2> {log}

        bgzip -@ {threads} {params.epibed} 2>> {log}
        tabix -p bed {output.epibed_gz} 2>> {log}
        """
