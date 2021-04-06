from os.path import join, abspath, expanduser, exists
localrules: bwa_index_setup, postprocess, label_bins, extract_DAStool, combine_final_reports

def get_sample_assemblies_reads(sample_file):
    sample_reads = {}
    sample_assemblies = {}
    with open(sample_file) as sf:
        for l in sf.readlines():
            s = l.strip().split("\t")
            if len(s) == 1 or s[0] == 'Sample' or s[0] == '#Sample' or s[0].startswith('#'):
                continue
            if len(s) < 3:
                sys.exit('Badly formatted sample_file')
            sample = s[0]
            assembly = s[1]
            reads_split = s[2].split(',')
            if sample in sample_reads:
                print(sample)
                raise ValueError("Non-unique sample encountered!")
            # get read pairs and singles from read specification
            if (len(reads_split) == 3) or (len(reads_split) == 2):
                sample_reads[sample] = reads_split[0:2]
            elif len(reads_split)==1:
                sample_reads[sample] = reads_split[0]
                # sys.exit('must be paired end reads')
            sample_assemblies[sample] = assembly
    return sample_reads, sample_assemblies

# Read in sample and outdir from config file
sample_file = config['sample_file']
outdir = config['outdir_base']
sample_reads, sample_assemblies = get_sample_assemblies_reads(sample_file)
sample_list = list(sample_reads.keys())


# are we using a non-standard (non ncbi) taxonomy
if 'custom_taxonomy' in config:
    custom_taxonomy = config['custom_taxonomy']
else:
    custom_taxonomy = False

# convert outdir to absolute path
if outdir[0] == '~':
    outdir = expanduser(outdir)
outdir = abspath(outdir)


# Determine if long reads
if 'long_read' in config and config['long_read']:
    long_read = True
else:
    long_read = False


# to speedup execution - remove samples that are done completely
# only if this new option specified in config
if 'skip_finished' in config:
    if config['skip_finished']:
        sample_list_new = [s for s in sample_list if not exists(join(outdir, s, "final", s + ".tsv")) ]
        sample_list = sample_list_new

# params from vamb 
INDEX_SIZE = config['index_size']
MM_MEM = config['minimap_mem']
MM_PPN = config['minimap_ppn']
VAMB_MEM = config['vamb_mem']
VAMB_PPN = config['vamb_ppn']
CHECKM_MEM = config['checkm_mem']
CHECKM_PPN = config['checkm_ppn']
SAMPLE_DATA = config['sample_data']
CONTIGS = config['contigs']
VAMB_PARAMS = config['vamb_params']
VAMB_PRELOAD = config['vamb_preload']

# list of per-sample assemblies for vamb
contigs_list = []
fh_in = open(CONTIGS, 'r')
for line in fh_in:
    line = line.rstrip()
    contigs_list.append(line)

print('##################################################################')
print(' SAMPLE LIST ')
print(sample_list)
print('##################################################################')

def get_metabat_bins(wildcards):
    outputs = checkpoints.metabat.get(**wildcards).output[0]
    return glob_wildcards(join(outputs, "{metabat_bin}.fa")).metabat_bin

def get_maxbin_bins(wildcards):
    outputs = checkpoints.maxbin.get(**wildcards).output[0]
    return glob_wildcards(join(outputs, "{maxbin_bin}.fasta")).maxbin_bin

def get_vamb_bins(wildcards):
    outputs = checkpoints.vamb.get(**wildcards).output[0]
    return glob_wildcards(join(outputs, "{vamb_bin}.fasta")).vamb_bin

def get_DAStool_bins(wildcards):
    outputs = checkpoints.extract_DAStool.get(**wildcards).output[0]
    # outputs = join(outdir, wildcards.sample, "DAS_tool_bins")
    return glob_wildcards(join(outputs, "{bin}.fa")).bin

rule all:
    input:
        ## miscellaneous
        config['kraken2db'],
        expand(join(outdir, "{sample}/idx/{sample}.fa.amb"), sample = sample_list),
        ## vamb specific
        join(outdir,jgi_matrix/jgi.abundance.dat"),
        join(outdir,"vamb/clusters.tsv"),
        join(outdir,"vamb/checkm.results"),
        #expand(join(outdir, "{sample}/vamb/checkm/checkm.results"), sample = config['sample']),
        ## Post-processing
        expand(join(outdir, "{sample}/classify/bin_species_calls.tsv"), sample = sample_list),
        expand(join(outdir, "{sample}/final/{sample}.tsv"), sample = sample_list),
        expand(join(outdir, "{sample}/final/{sample}_simple.tsv"), sample = sample_list),
        join(outdir, "binning_table_all_full.tsv"),


###############################################
####### Prepping input for vamb binning #######
###############################################

rule cat_contigs:
    input:
        contigs_list
    output:
        join(outdir,"/vamb/contigs.flt.fna.gz")
    params:
        walltime="864000", nodes="1", ppn="1", mem="10gb"
        outfolder= join(outdir,"/vamb/")
        
    threads:
        int(1)
    log:
        join(outdir,"/vamb/log/contigs/catcontigs.log")
    #conda:
    #    "envs/vamb.yaml"
    shell:"""
         mkdir -p {params.outfolder}
         /vol/sci/bio/data/moran.yassour/lab/Tools/vamb/bin/concatenate.py {output} {input} -m 2000
      """

rule index:
    input:
        contigs = join(outdir,"/vamb/contigs.flt.fna.gz")
    output:
        mmi = join(outdir,"/vamb/contigs.flt.mmi")
    params:
        walltime="864000", nodes="1", ppn="1", mem="90gb"
    threads:
        int(1)
    log:
        join(outdir,"/vamb/log/contigs/index.log")
    conda: 
        "/vol/sci/bio/data/moran.yassour/lab/Tools/vamb/workflow/envs/minimap2.yaml"
    shell:
        "minimap2 -I {INDEX_SIZE} -d {output} {input} 2> {log}"

rule dict:
    input:
        contigs = "contigs.flt.fna.gz"
    output:
        dict = "contigs.flt.dict"
    params:
        walltime="864000", nodes="1", ppn="1", mem="10gb"
    threads:
        int(1)
    log:
        join(outdir,"/vamb/log/contigs/dict.log")
    conda:
        "/vol/sci/bio/data/moran.yassour/lab/Tools/vamb/workflow/envs/samtools.yaml"
    shell:
        "samtools dict {input} | cut -f1-3 > {output} 2> {log}"

rule minimap:
    input:
        fq = lambda wildcards: sample2path[wildcards.sample],
        mmi = "contigs.flt.mmi",
        dict = "contigs.flt.dict"
    output:
        bam = temp("mapped/{sample}.bam")
    params:
        walltime="864000", nodes="1", ppn=MM_PPN, mem=MM_MEM
    threads:
        int(MM_PPN)
    log:
        join(outdir,"/vamb/log/map/{sample}.minimap.log")
    conda:
        "/vol/sci/bio/data/moran.yassour/lab/Tools/vamb/workflow/envs/minimap2.yaml"
    shell:
        '''minimap2 -t {threads} -ax sr {input.mmi} {input.fq} | grep -v "^@" | cat {input.dict} - | samtools view -F 3584 -b - > {output.bam} 2>{log}'''

rule sort:
    input:
        "mapped/{sample}.bam"
    output:
        temp("mapped/{sample}.sort.bam")
    params:
        walltime="864000", nodes="1", ppn="2", mem="15gb",
        prefix="mapped/tmp.{sample}"
    threads:
        int(2)
    log:
        join(outdir,"/vamb/log/map/{sample}.sort.log")
    conda:
        "/vol/sci/bio/data/moran.yassour/lab/Tools/vamb/workflow/envs/samtools.yaml"
    shell:
        "samtools sort {input} -T {params.prefix} --threads 1 -m 3G -o {output} 2>{log}"

rule jgi:
    input:
        bam = "mapped/{sample}.sort.bam"
    output:
        jgi = temp("jgi/{sample}.raw.jgi")
    params:
        walltime="864000", nodes="1", ppn="1", mem="10gb"
    threads:
        int(1)
    log:
        join(outdir,"/vamb/log/jgi/{sample}.jgi")
    conda:
        "/vol/sci/bio/data/moran.yassour/lab/Tools/vamb/workflow/envs/metabat2.yaml"
    shell:
        "jgi_summarize_bam_contig_depths --noIntraDepthVariance --outputDepth {output} {input} 2>{log}"

rule cut_column1to3: 
    input:
        "jgi/%s.raw.jgi" % IDS[0] 
    output:
        "jgi/jgi.column1to3"
    params:
        walltime="86400", nodes="1", ppn="1", mem="1gb"
    log:
        join(outdir,"/vamb/log/jgi/column1to3")
    shell: 
        "cut -f1-3 {input} > {output} 2>{log}"

rule cut_column4to5:
    input:
        "jgi/{sample}.raw.jgi"
    output:
        "jgi/{sample}.cut.jgi"
    params:
        walltime="86400", nodes="1", ppn="1", mem="1gb"
    log:
        join(outdir,"/vamb/log/jgi/{sample}.cut.log")
    shell: 
        "cut -f1-3 --complement {input} > {output} 2>{log}"

rule paste_abundances:
    input:
        column1to3="jgi/jgi.column1to3",
        data=expand("jgi/{sample}.cut.jgi", sample=IDS)
    output:
        "jgi_matrix/jgi.abundance.dat" 
    params:
        walltime="86400", nodes="1", ppn="1", mem="1gb"
    log:
        join(outdir,"/vamb/log/jgi/paste_abundances.log")
    shell: 
        "paste {input.column1to3} {input.data} > {output} 2>{log}" 

rule vamb:
    input:
        jgi = join(outdir,"/vamb/jgi_matrix/jgi.abundance.dat"),
        contigs = join(outdir,"/vamb/contigs.flt.fna.gz")
    output:
        join(outdir,"/vamb/clusters.tsv"),
        join(outdir,"/vamb/latent.npz"),
        join(outdir,"/vamb/lengths.npz"),
        join(outdir,"/vamb/log.txt"),
        join(outdir,"/vamb/model.pt"),
        join(outdir,"/vamb/mask.npz"),
        join(outdir,"/vamb/tnf.npz")
    params:
        vamb_folder = join(outdir,"/vamb")
        walltime="86400", nodes="1", ppn=VAMB_PPN, mem=VAMB_MEM
    log:
        join(outdir,"/vamb/log/vamb.log")
    threads:
        int(VAMB_threads)
    #conda:
        #"/vol/sci/bio/data/moran.yassour/lab/Tools/vamb/workflow/envs/vamb.yaml"
    shell:"""
        {VAMB_PRELOAD}
        #rm -rf vamb
        /vol/sci/bio/data/moran.yassour/lab/Tools/vamb/bin/vamb --outdir {params.vamb} --fasta {input.contigs} --jgi {input.jgi} {VAMB_PARAMS} 2>{log}
        """

#rule split_sample_bins:
#essential for running DAS_Tool as we are used to

##### to revise

rule checkm_vamb:
    input:
        join(outdir,"/vamb/clusters.tsv")
    output:
        join(outdir,"/vamb/checkm.tsv")
    params:
        bins = join(outdir,"/vamb/bins"),
        outdir = join(outdir,"/vamb/checkm")
    log:
        join(outdir, "/vamb/log/checkm.log")
    resources:
        mem = 128,
        time = 12
    threads: 4
    singularity:
        "shub://bsiranosian/bin_genomes:checkm"
    shell:"""
        checkm lineage_wf -f {output} --tab_table -t {threads} -x fna {params.bins} {params.outdir} 2>{log}
        """
        

#############################################################
####### Prepping input for metabat and maxbin binning #######
#############################################################

# copy assembly fasta file to working directory
rule bwa_index_setup:
    input:
        lambda wildcards: sample_assemblies[wildcards.sample]
    output:
        join(outdir, "{sample}/idx/{sample}.fa")
    resources:
        mem=1,
        time=1
    threads: 1
    shell: """
        cp {input} {output}
        """

# index assembly file
rule bwa_index:
    input:
        join(outdir, "{sample}/idx/{sample}.fa")
    output:
        join(outdir, "{sample}/idx/{sample}.fa.amb"),
        join(outdir, "{sample}/idx/{sample}.fa.ann"),
        join(outdir, "{sample}/idx/{sample}.fa.bwt"),
        join(outdir, "{sample}/idx/{sample}.fa.pac"),
        join(outdir, "{sample}/idx/{sample}.fa.sa")
    log:
        join(outdir, "{sample}/logs/bwa_index.log")
    conda: "envs/align.yaml"
    resources:
        mem=8,
        time=2
    threads: 1
    shell: """
        bwa index {input}
        """

# Align reads to assembly
rule bwa_align:
    input:
        reads = lambda wildcards: sample_reads[wildcards.sample],
        asm = join(outdir, "{sample}/idx/{sample}.fa"),
        amb = join(outdir, "{sample}/idx/{sample}.fa.amb"),
        ann = join(outdir, "{sample}/idx/{sample}.fa.ann"),
        bwt = join(outdir, "{sample}/idx/{sample}.fa.bwt"),
        pac = join(outdir, "{sample}/idx/{sample}.fa.pac"),
        sa = join(outdir, "{sample}/idx/{sample}.fa.sa")
    output:
        join(outdir, "{sample}/{sample}.bam")
    log:
        join(outdir, "{sample}/logs/bwa_mem.log")
    conda: "envs/align.yaml"
    resources:
        mem=16,
        time=12
    threads: 8
    shell: """
        bwa mem -t {threads} {input.asm} {input.reads} | samtools sort --threads {threads} > {output}
        """

# Index bam file to prepare for bamidx
rule bam_idx:
    input:
        join(outdir, "{sample}/{sample}_lr.bam") if long_read else join(outdir, "{sample}/{sample}.bam") #choose a long read alignment or short read alignment
    output:
        join(outdir, "{sample}/{sample}_lr.bam.bai") if long_read else join(outdir, "{sample}/{sample}.bam.bai") #choose a long read alignment or short read alignment
    log:
        join(outdir, "{sample}/logs/bamidx.log")
    conda: "envs/align.yaml"
    threads: 8
    resources:
        mem = 4,
        time = 1
    shell:
        "samtools index -@ {threads} {input}"
        
rule align_lr:
    input:
        join(outdir, "{sample}/idx/{sample}.fa"),
        reads = lambda wildcards: sample_reads[wildcards.sample]
    log:
        join(outdir, "{sample}/logs/align_lr.log")
    output:
        join(outdir, "{sample}/{sample}_lr.bam")
    conda: "envs/align.yaml"
    resources:
        mem=48,
        time=6
    threads: 16
    shell: """
        minimap2 -t {threads} -ax map-ont {input} | samtools sort --threads {threads} > {output}
        """

# Generate a depth file from BAM file for MetaBat input
rule metabat_pre:
    input:
        join(outdir, "{sample}/{sample}_lr.bam") if long_read else join(outdir, "{sample}/{sample}.bam") #choose a long read alignment or short read alignment
    output:
        single = join(outdir, "{sample}/{sample}.fa.depth.txt"),
    params:
        paired = join(outdir, "{sample}/{sample}.fa.paired.txt"),
    singularity:
        "docker://quay.io/biocontainers/metabat2:2.15--h137b6e9_0"
    shell: """
        jgi_summarize_bam_contig_depths --outputDepth {output.single} --pairedContigs {params.paired} --minContigLength 1000 --minContigDepth 1  {input} --percentIdentity 50
        """

#####################################################################
################ Binning with metabat and maxbin ####################
#####################################################################
checkpoint metabat:
    input:
        asm = join(outdir, "{sample}/idx/{sample}.fa"),
        depth = join(outdir, "{sample}/{sample}.fa.depth.txt"),
    output:
        directory(join(outdir, "{sample}/metabat/bins/"))
    singularity:
        "docker://quay.io/biocontainers/metabat2:2.15--h137b6e9_0"
    resources:
        mem=64,
        time=24
    threads: 4
    params:
        outstring = join(outdir, "{sample}/metabat/bins/bin")
    shell: """
        metabat2 --seed 1 -t {threads} --unbinned --inFile {input.asm} --outFile {params.outstring} --abdFile {input.depth}
        # if no bins produced, copy contigs to bin.unbinned
        if [ $(ls {output} | wc -l ) == "0" ]; then
            cp {input.asm} {output}/bin.unbinned.fa
        fi

        # check for bin.tooShort.fa thats empty
        if [ -f {output}/bin.tooShort.fa ]; then
            echo "Found bin.tooShort.fa"
            if [ $(cat {output}/bin.tooShort.fa | wc -l ) == "0" ]; then
                echo "Removing bin.tooShort.fa"
                rm {output}/bin.tooShort.fa
            fi
        fi
        """

checkpoint maxbin:
    input:
        contigs = join(outdir, "{sample}/idx/{sample}.fa"),
        depth = rules.metabat_pre.output.single
    output:
        directory(join(outdir, "{sample}/maxbin/bins"))
    singularity:
        "docker://quay.io/biocontainers/maxbin2:2.2.7--he1b5a44_1"
    params:
        outfolder = join(outdir, "{sample}/maxbin/"),
        logfile = join(outdir, "{sample}/maxbin/maxbin.log"),
        abundance_file = join(outdir, "{sample}/maxbin/{sample}_maxbin_depth.txt"),
    resources:
        time=lambda wildcards, attempt: attempt * 2
    threads: 8
    shell: """
        if [ -d {params.outfolder} ]; then rm -r {params.outfolder}; fi
        mkdir -p {params.outfolder}
        cd {params.outfolder}

        # create abundance file, which we already have from the metabat rule
        cut -f 1,4 {input.depth} | tail -n +2 > {params.abundance_file}
        run_MaxBin.pl -contig {input.contigs} -out maxbin \
        -abund {params.abundance_file} -thread {threads} || true
        
        # look in log  if dataset cant be binned, just copy contigs over
        if $(grep -q "cannot be binned" {params.logfile}); then
            echo 'DATASET CANNOT BE BINNED'
            mkdir {params.outfolder}/bins/
            cp {input.contigs} {params.outfolder}/bins/maxbin.unbinned.fasta
        elif ls {params.outfolder}/*.fasta 1> /dev/null 2>&1; then
            echo "FOUND BINS"
            mkdir {params.outfolder}/bins/
            mv {params.outfolder}/*.fasta {params.outfolder}/bins/
        else 
            echo "PROGRAM MUST HAVE FAILED"
            echo $(ls {params.outfolder})
            exit 1
        fi
        """

# checkm for metabat output
rule checkm_metabat:
    input:
        lambda wildcards: expand(join(outdir, "{sample}/metabat/bins/{metabat_bin}.fa"), metabat_bin = get_metabat_bins(wildcards), sample = wildcards.sample)
    output:
        join(outdir, "{sample}/metabat/checkm/checkm.tsv")
    log:
        join(outdir, "{sample}/metabat/logs/checkm.log")
    singularity:
        "shub://bsiranosian/bin_genomes:checkm"
    resources:
        mem = 128,
        time = 12
    threads: 4
    params:
        binfolder = join(outdir, "{sample}/metabat/bins/"),
        checkmfolder = join(outdir, "{sample}/metabat/checkm/"),
    shell: """
        rm -rf {sample}/metabat/checkm/*
        checkm lineage_wf -t {threads} -x fa --tab_table -f {output} {params.binfolder} {params.checkmfolder}
        """

# checkm for maxbin output
rule checkm_maxbin:
    input:
        lambda wildcards: expand(join(outdir, "{sample}/maxbin/bins/{maxbin_bin}.fasta"), maxbin_bin = get_maxbin_bins(wildcards), sample = wildcards.sample)
    output:
        join(outdir, "{sample}/maxbin/checkm/checkm.tsv")
    log:
        join(outdir, "{sample}/maxbin/logs/checkm.log")
    singularity:
        "shub://bsiranosian/bin_genomes:checkm"
    resources:
        mem = 128,
        time = 12
    threads: 4
    params:
        binfolder = join(outdir, "{sample}/maxbin/bins"),
        checkmfolder = join(outdir, "{sample}/maxbin/checkm/")
    shell: """
        rm -rf {sample}/maxbin/checkm/*
        checkm lineage_wf -t {threads} -x fasta --tab_table -f {output} {params.binfolder} {params.checkmfolder}
        """


# Aggregate binning results using DAStool
rule DAStool:
    input:
        lambda wildcards: expand(join(outdir, "{sample}/metabat/bins/{metabat_bin}.fa"), metabat_bin = get_metabat_bins(wildcards), sample = wildcards.sample),
        lambda wildcards: expand(join(outdir, "{sample}/maxbin/bins/{maxbin_bin}.fasta"), maxbin_bin = get_maxbin_bins(wildcards), sample = wildcards.sample),
        lambda wildcards: expand(join(outdir,"{sample}/vamb/bins/{vamb_bin}.fna"), vamb_bin = get_vamb_bins(wildcards), sample = wildcards.sample),
        contigs = join(outdir, "{sample}/idx/{sample}.fa"),
    output: 
        join(outdir, "{sample}/DAS_tool/completed.txt")
    conda: "envs/das_tool.yaml"
    params:
        outfolder = join(outdir, "{sample}/DAS_tool/"),
        metabat_dir = join(outdir, "{sample}/metabat/bins/"),
        maxbin_dir = join(outdir, "{sample}/maxbin/bins/"),
        vamb_dir = join(outdir, "{sample}/vamb/bins/"),
        metabat_tsv = join(outdir, "{sample}/DAS_tool/metabat_scaffold2bin.tsv"),
        maxbin_tsv = join(outdir, "{sample}/DAS_tool/maxbin_scaffold2bin.tsv"),
        vamb_tsv = join(outdir, "{sample}/DAS_tool/vamb_scaffold2bin.tsv"),
        logfile = join(outdir, "{sample}/DAS_tool/_DASTool.log"),
    threads: 8
    resources:
        time=lambda wildcards, attempt: attempt * 6
    shell: """
        # Prepare scaffold2bin file for each set of bins
        Fasta_to_Scaffolds2Bin.sh -e fa -i {params.metabat_dir} > {params.metabat_tsv}
        Fasta_to_Scaffolds2Bin.sh -e fasta -i {params.maxbin_dir} > {params.maxbin_tsv}
        Fasta_to_Scaffolds2Bin.sh -e fna -i {params.vamb_dir} > {params.vamb_tsv}
        # Run DAS_Tool
        DAS_Tool -i {params.metabat_tsv},{params.maxbin_tsv},{params.vamb_tsv} \
        -l metabat,maxbin,vamb -c {input.contigs} -o {params.outfolder} \
        --search_engine diamond --threads {threads} --write_bins 1 --write_unbinned 1 || true

        if $(grep -q "No bins with bin-score >0.5 found" {params.logfile}) || $(grep -q "single copy gene prediction using diamond failed" {params.logfile}); then
            echo 'DATASET CANNOT BE BINNED'
            mkdir {params.outfolder}/_DASTool_bins
            cp {input.contigs} {params.outfolder}/_DASTool_bins/unbinned.fa
            touch {output}
        elif ls {params.outfolder}/_DASTool_bins/*.fa 1> /dev/null 2>&1; then
            echo "FOUND BINS"
            touch {output}
        else 
            echo "PROGRAM MUST HAVE FAILED IN SOME OTHER WAY"
            exit 1
        fi
        """

# extract the bins as a separate rule to speed up execution
checkpoint extract_DAStool:
    input: rules.DAStool.output
    output:
        directory(join(outdir, "{sample}/DAS_tool_bins"))
    params:
        old_binfolder = join(outdir, "{sample}/DAS_tool/_DASTool_bins"),
        new_binfolder = join(outdir, "{sample}/DAS_tool_bins"),
    shell: """
        mkdir -p {params.new_binfolder}
        cp {params.old_binfolder}/*.fa {params.new_binfolder}
    """


# checkm for DAStool output
rule checkm_DAStool:
    input:
        lambda wildcards: expand(join(outdir, "{sample}/DAS_tool_bins/{bin}.fa"), bin = get_DAStool_bins(wildcards), sample = wildcards.sample)
    output:
        join(outdir, "{sample}/DAS_tool/checkm/checkm.tsv")
    log:
        join(outdir, "{sample}/DAS_tool/logs/checkm.log")
    singularity:
        "shub://bsiranosian/bin_genomes:checkm"
    resources:
        mem = 128,
        time = 12
    threads: 4
    params:
        binfolder = join(outdir, "{sample}/DAS_tool_bins"),
        checkmfolder = join(outdir, "{sample}/DAS_tool/checkm/"),
        bin_ex = ".fa"
    shell: """
        rm -rf {params.checkmfolder}
        checkm lineage_wf -t {threads} -x {params.bin_ex} --tab_table -f {output} {params.binfolder} {params.checkmfolder}
        """

#####################################################
############ ANALYSIS OF DAS_TOOL BINS ##############
#####################################################

# Use aragorn to detect tRNA and tmRNA genes
rule aragorn:
    input:
        join(outdir, "{sample}/DAS_tool_bins/{bin}.fa")
    output:
        join(outdir, "{sample}/rna/trna/{bin}.fa.txt")
    log:
        join(outdir, "{sample}/logs/aragorn_{bin}.log")
    singularity:
        "docker://quay.io/biocontainers/prokka:1.14.6--pl526_0"
    resources:
        mem = 8,
        time = 1
    shell: """
        aragorn -t {input} -o {output}
        """

# Use barrnap to predict ribosomal RNA
rule barrnap:
    input:
        join(outdir, "{sample}/DAS_tool_bins/{bin}.fa")
    output:
        join(outdir, "{sample}/rna/rrna/{bin}.fa.txt")
    log:
        join(outdir, "{sample}/logs/barrnap_{bin}.log")
    singularity:
        "docker://quay.io/biocontainers/prokka:1.14.6--pl526_0"
    resources:
        mem = 8,
        time = 1
    shell: """
        barrnap {input} > {output}
        """

# Use quast to evaluate genome assemblies
rule quast:
    input:
        join(outdir, "{sample}/DAS_tool_bins/{bin}.fa")
    output:
        join(outdir, "{sample}/quast/{bin}.fa/report.tsv")
    log:
        join(outdir, "{sample}/logs/quast_{bin}.log")
    resources:
        mem=8,
        time=1
    params:
        quastfolder = join(outdir, "{sample}/quast/{bin}/"),
        thresholds = "0,10000,50000,100000,250000,500000,1000000,2000000,3000000"
    shell: """
        mkdir -p {params.quastfolder}
        python /vol/sci/bio/data/moran.yassour/lab/Tools/quast-5.0.2/quast.py -o {params.quastfolder} {input} --contig-thresholds {params.thresholds} --fast
        """

rule pull_prokka:
    input: join(outdir, "{sample}/idx/{sample}.fa"),
    output: join(outdir, "{sample}/prokka/pulled.txt") 
    params:
        prokka_folder = join(outdir, "{wildcards.sample}/prokka")
    conda:
         "/vol/sci/bio/data/moran.yassour/lab/Tools/Anaconda/envs/envs_2/prokka.yml"
    shell: """
        mkdir -p {params.prokka_folder}
        touch {output}
    """

# Use prokka to annotate metagenomes
rule prokka:
    input:
        join(outdir, "{sample}/DAS_tool_bins/{bin}.fa"),
        rules.pull_prokka.output
    output:
        join(outdir, "{sample}/prokka/{bin}.fa/{sample}_{bin}.fa.gff")
    log:
        join(outdir, "{sample}/logs/prokka_{bin}.log")
    conda: 
        "/vol/sci/bio/data/moran.yassour/lab/Tools/Anaconda/envs/envs_2/prokka.yml"
    resources:
        mem = 48,
        time = lambda wildcards, attempt: 4 * attempt,
    threads: 8
    params:
        prokka_folder = join(outdir, "{sample}/prokka/{bin}"),
        prefix = "{sample}_{bin}.fa"
    shell: """
        # don't run this on unbinned contigs, takes forever
        if [[ {wildcards.bin} =~ "unbinned" ]]; then
            touch {output}
            touch {params.prokka_folder}/prokka_skipped.out
        else
            prokka {input} --outdir {params.prokka_folder} --prefix {params.prefix} \
            --addgenes --mincontiglen 200 --force --cpus {threads} --noanno
            #--centre X --compliant 
        fi
        """

# Retrieve stats on mapping from whole metagenome sample
rule bam_idxstats:
    input:
        join(outdir, "{sample}/{sample}_lr.bam") if long_read else join(outdir, "{sample}/{sample}.bam"), #choose a long read alignment or short read alignment,
        join(outdir, "{sample}/{sample}_lr.bam.bai") if long_read else join(outdir, "{sample}/{sample}.bam.bai"), #choose a long read alignment or short read alignment,
    output:
        join(outdir, "{sample}/{sample}_lr.bam.bai.tsv") if long_read else join(outdir, "{sample}/{sample}.bam.bai.tsv"), #choose a long read alignment or short read alignment,
    log:
        join(outdir, "{sample}/logs/bamidxstats.log")
    conda: "envs/align.yaml"
    resources:
        mem = 2,
        time = 2
    shell:
        "samtools idxstats {input[0]} > {output}"

# Get mapping stats for each bin
rule bin_idxstats:
    input:
        join(outdir, "{sample}/DAS_tool_bins/{bin}.fa"),
        join(outdir, "{sample}/{sample}_lr.bam.bai.tsv") if long_read else join(outdir, "{sample}/{sample}.bam.bai.tsv"), #choose a long read alignment or short read alignment,
    output:
        join(outdir, "{sample}/coverage/raw/{bin}.tsv")
    log:
        join(outdir, "{sample}/logs/coverage_idxstats_{bin}.log")
    resources:
        mem = 2,
        time = 6
    shell: """
        grep '>' {input[0]} | tr -d '>' | cut -f 1 -d " " | xargs -I foo -n 1 grep -P 'foo\t' {input[1]} > {output}
    """

# Determine coverage for each bin
rule bin_coverage:
    input:
        rules.bin_idxstats.output
    output:
        join(outdir, "{sample}/coverage/{bin}.txt")
    log:
        join(outdir, "{sample}/logs/coverage_{bin}.log")
    conda: "envs/assign_species.yaml"
    resources:
        mem = 2,
        time = 1
    params:
        read_length = config['read_length'],
        sample = lambda wildcards: wildcards.sample
    script:
        "scripts/bin_coverage.py"

# index DAStool bins
rule fasta_index:
    input:
        join(outdir, "{sample}/DAS_tool_bins/{bin}.fa")
    output:
        join(outdir, "{sample}/DAS_tool_bins/{bin}.fa.fai")
    log:
        join(outdir, "{sample}/logs/faidx_{bin}.log")
    conda: "envs/align.yaml"
    resources:
        mem = 8,
        time = 1
    threads: 1
    shell:
        "samtools faidx {input}"

# Use kraken2 to assign taxonomic labels
rule kraken2:
    input:
        join(outdir, "{sample}/idx/{sample}.fa")
    output:
        krak = join(outdir, "{sample}/classify/{sample}.krak"),
        krak_report = join(outdir, "{sample}/classify/{sample}.krak.report")
    log:
        join(outdir, "{sample}/logs/kraken_class.log")
    
    params:
        db = config['kraken2db']
    resources:
        mem = 256,
        time = 6
    threads: 4
    shell: """
        /vol/sci/bio/data/moran.yassour/lab/Tools/kraken2/kraken2 --db {params.db} --db {params.db} --threads {threads} \
        --output {output.krak} --report {output.krak_report} {input}
        """

rule label_bins:
    input:
        krak = join(outdir, "{sample}/classify/{sample}.krak"),
        bins = lambda wildcards: expand(join(outdir, "{sample}/DAS_tool_bins/{bin}.fa.fai"),
                                        bin = get_DAStool_bins(wildcards), sample = wildcards.sample)
    output:
        join(outdir, "{sample}/classify/bin_species_calls.tsv")
    log:
        join(outdir, "{sample}/logs/assign_species.log")
    conda: "envs/assign_species.yaml"
    params:
        binfolder = join(outdir, "{sample}/DAS_tool_bins/"),
        custom_taxonomy = custom_taxonomy
    script:
        "scripts/assign_species.py"

rule postprocess:
    input:
        prokka = lambda wildcards: expand(join(outdir, "{sample}/prokka/{bin}.fa/{sample}_{bin}.fa.gff"), bin = get_DAStool_bins(wildcards), sample = wildcards.sample),
        quast = lambda wildcards: expand(join(outdir, "{sample}/quast/{bin}.fa/report.tsv"), bin = get_DAStool_bins(wildcards), sample = wildcards.sample),
        checkm = join(outdir, "{sample}/DAS_tool/checkm/checkm.tsv"),
        trna = lambda wildcards: expand(join(outdir, "{sample}/rna/trna/{bin}.fa.txt"), bin = get_DAStool_bins(wildcards), sample = wildcards.sample),
        rrna = lambda wildcards: expand(join(outdir, "{sample}/rna/rrna/{bin}.fa.txt"), bin = get_DAStool_bins(wildcards), sample = wildcards.sample),
        classify = rules.label_bins.output,
        coverage = lambda wildcards: expand(join(outdir, "{sample}/coverage/{bin}.txt"), bin = get_DAStool_bins(wildcards), sample = wildcards.sample),
    output:
        full = join(outdir, "{sample}/final/{sample}.tsv"),
        simple = join(outdir, "{sample}/final/{sample}_simple.tsv")
    singularity: "docker://r-base"
    params:
        bins = lambda wildcards: get_DAStool_bins(wildcards),
        # bins = glob_wildcards(join(outdir, "{sample}/DAS_tool_bins", "{bin}.fa")).bin,
        sample = lambda wildcards: wildcards.sample
    script: 
        "scripts/postprocess.R"

rule combine_final_reports:
    input:
        all_full = expand(join(outdir, "{sample}/final/{sample}.tsv"), sample=sample_list),
        single_full = expand(join(outdir, "{sample}/final/{sample}.tsv"), sample=sample_list[0]),
        all_simple = expand(join(outdir, "{sample}/final/{sample}_simple.tsv"), sample=sample_list),
        single_simple = expand(join(outdir, "{sample}/final/{sample}_simple.tsv"), sample=sample_list[0]),
    output:
        full = join(outdir, "binning_table_all_full.tsv"),
        simple = join(outdir, "binning_table_all_simple.tsv"),
    shell: """
        head -n 1 {input.single_full} > {output.full}
        tail -n +2 -q {input.all_full} >> {output.full}
        head -n 1 {input.single_simple} > {output.simple}
        tail -n +2 -q {input.all_simple} >> {output.simple}
    """
