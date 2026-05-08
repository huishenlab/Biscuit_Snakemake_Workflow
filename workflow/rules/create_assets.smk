###-----------------------------------------------------------------------------------------------------------------###
# Notes on where variables are defined, if not defined in obs_exp.smk
#
# output_directory      - workflow/Snakefile
# config                - workflow/Snakefile
# get_biscuit_reference - workflow/rules/biscuit.smk
#
###-----------------------------------------------------------------------------------------------------------------###

# Verify the provided genome is one of the available genomes
def check_genome(wildcards):
    if config['create_assets']['run'] and config['create_assets']['genome'] not in ['hg38', 'hg19', 'mm10', 'mm9']:
        print('hg38, hg19, mm10, and mm9 are the only available genomes for creating asset files')
        exit(1)
    elif config['create_assets']['run'] and config['create_assets']['genome'] in ['hg38', 'hg19', 'mm10', 'mm9']:
        return config['create_assets']['genome']

# Determine which GENCODE annotation file to download
# TODO: Can add mm39 as an option (release_M38 as of 8 May 2026)
def get_annotation_file(wildcards):
    base = 'http://ftp.ebi.ac.uk/pub/databases/gencode'
    gen = check_genome(wildcards)
    if gen == 'hg38':
        return f'{base}/Gencode_human/release_38/gencode.v38.annotation.gtf.gz'
    elif gen == 'hg19':
        return f'{base}/Gencode_human/release_19/gencode.v19.annotation.gtf.gz'
    elif gen == 'mm9':
        return f'{base}/Gencode_mouse/release_M1/gencode.vM1.annotation.gtf.gz'
    elif gen == 'mm10':
        return f'{base}/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz'

rule general_assets:
    input:
        ref = get_biscuit_reference,
    output:
        cpg = f'assets/{config["create_assets"]["genome"]}/cpg.bed.gz',
        bot = f'assets/{config["create_assets"]["genome"]}/windows100bp.gc_content.bot10p.bed.gz',
        top = f'assets/{config["create_assets"]["genome"]}/windows100bp.gc_content.top10p.bed.gz',
        cgi = f'assets/{config["create_assets"]["genome"]}/cpg_islands.bed.gz',
        opn = f'assets/{config["create_assets"]["genome"]}/cpg_open_seas.bed.gz',
        shl = f'assets/{config["create_assets"]["genome"]}/cpg_shelves.bed.gz',
        sho = f'assets/{config["create_assets"]["genome"]}/cpg_shores.bed.gz',
        exn = f'assets/{config["create_assets"]["genome"]}/exon.bed.gz',
        gen = f'assets/{config["create_assets"]["genome"]}/genic.bed.gz',
        inr = f'assets/{config["create_assets"]["genome"]}/intergenic.bed.gz',
        msk = f'assets/{config["create_assets"]["genome"]}/rmsk.bed.gz',
    params:
        annt = get_annotation_file,
        gen = check_genome,
        dir = f'assets/{config["create_assets"]["genome"]}',
        tmp_sho = f'assets/{config["create_assets"]["genome"]}/shores.tmp.bed',
        tmp_shl = f'assets/{config["create_assets"]["genome"]}/shelves.tmp.bed',
        tmp_gen = f'assets/{config["create_assets"]["genome"]}/gene.tmp.bed',
        tmp_exn = f'assets/{config["create_assets"]["genome"]}/exon.tmp.bed',
    log:
        f'{output_directory}/logs/general_assets.log',
    benchmark:
        f'{output_directory}/benchmarks/general_assets.txt',
    threads: 1
    resources:
        mem_gb = config['hpc_parameters']['intermediate_memory_gb'],
        time = config['runtime']['medium']
    conda:
        '../envs/biscuit.yaml'
    envmodules:
        config['envmodules']['bedtools'],
        config['envmodules']['biscuit'],
    shell:
        """
        set +o pipefail;

        echo "Create output directory" 1> {log}
        mkdir -p {params.dir}
        
        echo "Set reference location" 1> {log}
        # Reference location
        REFLOC={input.ref}.fai

        echo "Retrieve CpG island file" 1> {log}
        # CpG Islands from UCSC on only the canonical chromosomes
        wget -qO- http://hgdownload.cse.ucsc.edu/goldenpath/{params.gen}/database/cpgIslandExt.txt.gz | \
        gunzip -c | \
        awk 'BEGIN{{ OFS="\\t"; }}{{ print $2, $3, $4; }}' | \
        awk '{{ if ($1 ~ /^chr[1234567890XYM]{{1,2}}$/) {{ print }} }}' | \
        bedtools sort -i - | \
        gzip -c > {output.cgi}
        echo "CpG island file retrieved" 1> {log}

        # Create middles for finding locations
        bedtools slop \
            -i {output.cgi} \
            -g ${{REFLOC}} \
            -b 2000 | \
        bedtools merge -i - \
        > {params.tmp_sho}

        bedtools slop \
            -i {output.cgi} \
            -g ${{REFLOC}} \
            -b 4000 | \
        bedtools merge -i - \
        > {params.tmp_shl}

        # CpG Open Seas (intervening locations)
        sort -k1,1 -k2,2n ${{REFLOC}} | \
        bedtools complement -L -i {params.tmp_shl} -g - | \
        gzip -c > {output.opn}

        # CpG Shelves (shores +/- 2kb)
        bedtools subtract \
            -a {params.tmp_shl} \
            -b {params.tmp_sho} | \
        gzip -c > {output.shl}

        # CpG Shores (island +/- 2kb)
        bedtools subtract \
            -a {params.tmp_sho} \
            -b {output.cgi} | \
        gzip -c > {output.sho}

        # Clean up CpG temporary files
        rm -f {params.tmp_shl} {params.tmp_sho}
        echo "CpG seascape created and cleaned up" 1> {log}

        # BISCUIT assets (includes CpG and top/bottom 10% GC-content files)
        build_biscuit_QC_assets.pl --outdir {params.dir} --ref {input.ref} 1> {log}

        # Repeat-masked file
        if [[ {params.gen} != "mm9" ]]; then
            wget -qO- http://hgdownload.cse.ucsc.edu/goldenpath/{params.gen}/database/rmsk.txt.gz | \
            gunzip -c | \
            awk 'BEGIN{{ OFS="\\t"; }}{{ print $6, $7, $8; }}' | \
            awk '{{ if ($1 ~ /^chr[1234567890XYM]{{1,2}}$/) {{ print }} }}' | \
            bedtools sort -i - | \
            bedtools merge -i - | \
            gzip -c > {output.msk}
        else
            wget -q -nd -r -nH -np --accept-regex 'chr[1234567890XYM]*_rmsk.txt.gz$' \
                --regex-type=posix http://hgdownload.cse.ucsc.edu/goldenpath/mm9/database/

            cat chr*_rmsk.txt.gz | \
            gunzip -c | \
            awk 'BEGIN{{ OFS="\\t"; }}{{ print $6, $7, $8; }}' | \
            bedtools sort -i - | \
            bedtools merge -i - | \
            gzip -c > {output.msk}

            if [[ -f index.html ]]; then rm index.html; fi
            if [[ -f robots.txt ]]; then rm robots.txt; fi
            rm chr*_rmsk.txt.gz
        fi

        # Genic, intergenic, and exon files
        wget -qO- {params.annt} | \
        gunzip -c | \
        awk 'BEGIN{{ OFS="\\t"; }}{{ if ($3 == "gene") {{ print $1, $4, $5 > "{params.tmp_gen}"; }} else if ($3 == "exon") {{ print $1, $4, $5 > "{params.tmp_exn}"; }} }}'

        bedtools sort -i {params.tmp_gen} | \
        bedtools merge -i - | \
        gzip -c > {output.gen}

        bedtools sort -i {params.tmp_exn} | \
        bedtools merge -i - | \
        gzip -c > {output.exn}

        sort -k1,1 -k2,2n {input.ref}.fai | \
        bedtools complement -L -i {output.gen} -g - | \
        gzip -c > {output.inr}

        # Clean up
        rm -f {params.tmp_gen} {params.tmp_exn}
        """       

rule bismap_assets:
    input:
        cpg = f'assets/{config["create_assets"]["genome"]}/cpg.bed.gz',
        cgi = f'assets/{config["create_assets"]["genome"]}/cpg_islands.bed.gz',
        exn = f'assets/{config["create_assets"]["genome"]}/exon.bed.gz',
        gen = f'assets/{config["create_assets"]["genome"]}/genic.bed.gz',
        inr = f'assets/{config["create_assets"]["genome"]}/intergenic.bed.gz',
        msk = f'assets/{config["create_assets"]["genome"]}/rmsk.bed.gz',
        ref = get_biscuit_reference,
    output:
        bis = f'assets/{config["create_assets"]["genome"]}/k100.bismap.bedgraph.gz',
        sor = f'assets/{config["create_assets"]["genome"]}/k100.bismap.sorted.bed.gz',
        wgt = f'assets/{config["create_assets"]["genome"]}/k100.bismap.10kb_avg.bed.gz',
        cpg = f'assets/{config["create_assets"]["genome"]}/cpg_bismap.bed.gz',
        cgi = f'assets/{config["create_assets"]["genome"]}/cgi_bismap.bed.gz',
        exn = f'assets/{config["create_assets"]["genome"]}/exon_bismap.bed.gz',
        gen = f'assets/{config["create_assets"]["genome"]}/genic_bismap.bed.gz',
        inr = f'assets/{config["create_assets"]["genome"]}/intergenic_bismap.bed.gz',
        msk = f'assets/{config["create_assets"]["genome"]}/rmsk_bismap.bed.gz',
        byn = f'assets/{config["create_assets"]["genome"]}/genome_window_10kb.bed.gz',
    params:
        gen = check_genome,
        dir = f'assets/{config["create_assets"]["genome"]}',
    log:
        f'{output_directory}/logs/bismap_assets.log',
    benchmark:
        f'{output_directory}/benchmarks/bismap_assets.txt',
    threads: 8
    resources:
        mem_gb = config['hpc_parameters']['intermediate_memory_gb'],
        time = config['runtime']['medium']
    conda:
        '../envs/biscuit.yaml'
    envmodules:
        config['envmodules']['bedtools'],
    shell:
        """
        set +o pipefail;

        wget -P {params.dir} --no-check-certificate -q https://bismap.hoffmanlab.org/raw/{params.gen}/k100.bismap.bedgraph.gz

        # Create average mappability scores in 10kb windows
        zcat {output.bis} | sort --parallel={threads} -k1,1 -k2,2n -T {params.dir} | gzip > {output.sor}

        bedtools makewindows -w 10000 -g {input.ref}.fai | gzip 1> {output.byn} 2> {log}
        bedtools map -a {output.byn} -b {output.sor} -c 4 -o mean | gzip 1> {output.wgt} 2>> {log}

        bedtools intersect -a {output.bis} -b {input.cpg} | gzip -c 1> {output.cpg} 2>> {log}
        bedtools intersect -a {output.bis} -b {input.cgi} | gzip -c 1> {output.cgi} 2>> {log}
        bedtools intersect -a {output.bis} -b {input.exn} | gzip -c 1> {output.exn} 2>> {log}
        bedtools intersect -a {output.bis} -b {input.gen} | gzip -c 1> {output.gen} 2>> {log}
        bedtools intersect -a {output.bis} -b {input.inr} | gzip -c 1> {output.inr} 2>> {log}
        bedtools intersect -a {output.bis} -b {input.msk} | gzip -c 1> {output.msk} 2>> {log}
        """

rule wcgw_assets:
    input:
        cpg = f'assets/{config["create_assets"]["genome"]}/cpg.bed.gz',
        ref = get_biscuit_reference,
    output:
        expand(f'assets/{config["create_assets"]["genome"]}/{{ctxt}}_{{ncpgs}}_neighbors.bed.gz',
            ctxt=['scgs', 'scgw', 'wcgw'],
            ncpgs=['0', '1', '2', '3p']
        ),
    params:
        dir = f'assets/{config["create_assets"]["genome"]}',
        cpg = os.getcwd(),
    log:
        f'{output_directory}/logs/wcgw_assets.log',
    benchmark:
        f'{output_directory}/benchmarks/wcgw_assets.txt',
    threads: config['hpc_parameters']['wcgw_threads']
    resources:
        mem_gb = config['hpc_parameters']['small_memory_gb'],
        time = config['runtime']['short'],
    conda:
        '../envs/biscuit.yaml'
    envmodules:
        config['envmodules']['bedtools'],
        config['envmodules']['parallel'],
    shell:
        """
        bash workflow/scripts/create_context_beds.sh \
            -t {threads} \
            -o {params.dir} \
            {input.ref} \
            {params.cpg}/{input.cpg} 2> {log}
        """
