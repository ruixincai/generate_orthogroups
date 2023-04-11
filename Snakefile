#!/usr/bin/env python3

from pathlib import Path

proteomes_directory = 'data/proteomes'

bbmap = 'docker://quay.io/biocontainers/bbmap:39.01--h5c4e2a8_0'

pggb = 'docker://ghcr.io/pangenome/pggb:20230113201558a9a04c'

vg = 'docker://quay.io/vgteam/vg:v1.46.0'

ref_proteomes = {
    'aalb': Path(proteomes_directory,
                 'GCF_006496715.1_Aalbo_primary.1_protein.faa'),
    'agam': Path(proteomes_directory,
                 'GCF_000005575.2_AgamP3_protein.faa'),
    'cpip': Path(proteomes_directory,
                 'GCF_016801865.1_TS_Cpip_V1_protein.faa')
}

gff_directory = 'data/gff'

ref_gff = {
    'aalb': Path(gff_directory,
        'GCF_006496715.1_Aalbo_primary.1_genomic.gff'),
    'agam': Path(gff_directory,
        'GCF_000005575.2_AgamP3_genomic.gff'),
    'cpip': Path(gff_directory,
        'GCF_016801865.1_TS_Cpip_V1_genomic.gff')
}

transcript_directory = 'data/transcripts'

ref_transcripts = {
    'aalb': Path(transcript_directory,
        'GCF_006496715.1_Aalbo_primary.1_rna.fna'),
    'agam': Path(transcript_directory,
        'GCF_000005575.2_AgamP3_rna.fna'),
    'cpip': Path(transcript_directory,
        'GCF_016801865.1_TS_Cpip_V1_rna.fna')
}

# leave this as the first rule
rule target:
    input:
        'output/orthoresult/Results_Results/Orthogroups/Orthogroups.tsv',
        expand('output/separated_transcript_protein_tables/{species}.csv',
               species=list(ref_gff.keys())),
        'output/combined_table.csv',
        'output/ortho_match.csv',
        'output/separated_orthogroups',
        "data/transcripts/combined_transcripts.fa"



# Run the OrthoFinder on the reference protein sequences of three mosquito species in FASTA format
rule orthofinder:
    input:
        input_files = ref_proteomes.values() # Path()
    output:
        orthoresult = 'output/orthoresult/Results_Results/Orthogroups/Orthogroups.tsv'
    log:
        'output/logs/orthofinder.log'
    params:
        input_dir = proteomes_directory,
        output_dir = 'output/orthoresult',
        output = 'Results'
    threads:
        workflow.cores
    resources:
        mem_mb = '20480',
        time = '0-12:0:00' # easier to use integer minutes (i.e. 720)
    container:
        'docker://quay.io/biocontainers/orthofinder:2.5.4--hdfd78af_0'
    shell:
        'rm -r {params.output_dir} && ' # change directory to output_dir 
        'orthofinder '
        '-f {params.input_dir} '
        '-o {params.output_dir} '
        '-n {params.output} '
        '-t {threads} '
        '&> {log}'


# the protein accession number and the transcript accession number was paired separately 
# in a single table for each species of mosquitoes
rule separated_transcript_protein_tables:
    input:
        gff = lambda wildcards: ref_gff[wildcards.species]
    output:
        output_tables = 'output/separated_transcript_protein_tables/{species}.csv'
    log:
        'output/logs/separated_transcript_protein_tables_{species}.log'
    threads:
        1
    resources:
        time = '0-0:20:00'
    container:
        'docker://ghcr.io/tomharrop/r-containers:bioconductor_3.17'
    script:
        'src/separated_transcript_protein.R'


# 3 CSV files generated before was combined into a single table database
# a combined tabular database include all transcript accession numbers, protein accession numbers for 
# 3 mosquito species is generated
rule combined_transcript_protein_table:
    input:
        input_files = expand('output/separated_transcript_protein_tables/{species}.csv', 
            species=list(ref_gff.keys()))
    output:
        output_file = 'output/combined_table.csv'
    log:
        'output/logs/combined_transcript_protein_table.log'
    threads:
        1
    resources:
        time = '0-0:20:00'
    container:
        'docker://ghcr.io/tomharrop/r-containers:bioconductor_3.17' # r-bundle-bioconductor/3.12-r-4.0.4
    script:
        'src/combined_transcript_protein.R'


# Match the combined tabular database to the orthogroups generated from OrthoFinder
rule ortho_match:
    input:
        input_orthofinder = 'output/orthoresult/Results_Results/Orthogroups/Orthogroups.tsv',
        input_combinedtable = 'output/combined_table.csv'
    output:
        output_file = 'output/ortho_match.csv'
    log:
        'output/logs/ortho_match.log'
    threads:
        1
    resources:
        time = '0-0:30:00'
    container:
        'docker://ghcr.io/tomharrop/r-containers:bioconductor_3.17' # r-bundle-bioconductor/3.12-r-4.0.4
    script:
        'src/ortho_match.R'


# an number of file is created in a defined directory
checkpoint separated_orthogroups:
    input:
        input_file = 'output/ortho_match.csv'
    output:
        output_dir = directory('output/separated_orthogroups')
    log:
        'output/logs/separated_orthogroups.log' 
    threads:
        1
    resources:
        time = '0-0:20:00'
    container:
        'docker://ghcr.io/tomharrop/r-containers:bioconductor_3.17' # r-bundle-bioconductor/3.12-r-4.0.4
    script:
        'src/separate_ortho.R'


def aggregate(wildcards):
    # look up the ouptputs for the checkpoint
    checkpoint_output = checkpoints.separated_orthogroups.get(**wildcards).output['output_dir']
    # glob the output
    og_txt_path = 'output/separated_orthogroups/{og}.txt'
    all_og_ids = glob_wildcards(og_txt_path).og
    og_seq_output = "output/og_seq/{og}.fa"
    # we MUST return a list of files we want
    return(expand(og_seq_output,
        og = all_og_ids))


# combine 3 reference transcript sequences from three mosquito species into a single FASTA file
rule transcripts_merge:
    input:
        ref_transcripts.values()
    output:
        "data/transcripts/combined_transcripts.fa"
    shell:
        "cat {input} > {output}"

# the checkpoint that shall trigger re-evaluation of the DAG
rule og_seq:
    input:
        og_txt = "output/separated_orthogroups/{og}.txt",
        ref_trans = "data/transcripts/combined_transcripts.fa"
    output:
        "output/og_seq/{og}.fa"
    log:
        'output/logs/og_seq/og_seq.{og}.log'
    resources:
        time = '0-0:1:00'
    container:
        bbmap
    shell:
        'filterbyname.sh '
        'include=true '
        'prefix=true '
        'in={input.ref_trans} ' 
        'names={input.og_txt} ' 
        'out={output} '
        '&>{log}'



rule pggb_index:
    input:
        og_seq_fa = 'output/og_seq/{og}.fa'
    output:
        og_seq_index = 'output/og_seq/{og}.fa.gz',
        index = 'output/og_seq/{og}.fa.gz.fai'
    log:
        'output/logs/pggb_index/{og}.log'
    resources:
        time = '0-0:1:00'
    container: 
        pggb
    shell: 
        'bgzip -c {input.og_seq_fa}  > {output.og_seq_index} 2> {log} '
        '&& '
        'samtools faidx {output.og_seq_index} &>> {log}'


# pggb: Construct graphs for each orthogroup from transcript sequences of three mosquito species
checkpoint pggb:
    input:
        og_seq_index = 'output/og_seq/{og}.fa.gz',
        index = 'output/og_seq/{og}.fa.gz.fai'
    output:
        output_dir = directory('output/pggb/{og}.{identity}.{segment}')
    log:
        'output/logs/pggb/pggb.{og}.{identity}.{segment}.log'
    resources:
        time = '0-0:1:00'
    container: 
        pggb
    threads:
        1
    shell:
        'pggb '
        '-i {input.og_seq_index} '
        '-o {output} '
        # number of haplotypes
        '-n $( zgrep -c "^>" {input.og_seq_index} ) '
        # number of threads
        '-t {threads} '
        # segment length for scaffolding the graph
        '-s {wildcards.segment} '
        # pairwise identity 
        '-p {wildcards.identity} '
        # pruning matches shorter than a given threshold from the initial graph model
        # '-k '
        '&>{log}'

# target orthogroups
list_of_ogs = ['OG0012812','OG0009228','OG0007096','OG0009197','OG0010256','OG0012242','OG0004074','OG0000133', 'OG0008785','OG0001722']

rule pggb_test:
    input:
        expand('output/pggb/{og}.{identity}.{segment}',
            og=list_of_ogs,
            identity=[85, 90, 95, 60, 70, 80],
            segment=[100, 300, 3000]) # minimum segment length is required to be >= 100 bp

# rule vg_convert:
    # input:
        # gfa_input
    # output:
        # 'output/vg/{og}.{identity}.{segment}'

def gfa_input(wildcards):
    # pggb_dir = f'output/pggb/{wildcards.og}.{wildcards.identity}.{wildcards.segment}/{{filename}}.gfa'
    checkpoint_output = checkpoints.pggb.get(**wildcards).output['output_dir']
    gfa_file_path = f'{checkpoint_output}/{{filename}}.gfa'
    gfa_file = glob_wildcards(gfa_file_path).filename
    # gfa_output = str("output/pggb/{wildcards.og}.{wildcards.identity}.{wildcards.segment}/{{filename}}.gfa")
    return(expand(gfa_file_path,filename = gfa_file))


rule vg_test: 
    input:
        expand('output/vg/index/pruned/{og}.{identity}.{segment}.pruned.vg',
            og=list_of_ogs,
            identity=[85, 90, 95, 60, 70, 80],
            segment=[100, 300, 3000])

rule gfa2vg:
    input:
        gfa_input
    output:
        'output/vg/vg_version/{og}.{identity}.{segment}.vg'
    log:
        'output/logs/vg/vg_version.{og}.{identity}.{segment}.log'
    resources:
        time = '0-0:1:00'
    container: 
        vg
    shell:
        'vg convert '
        '-g {input} '
        '-v '
        '> {output} '
        '2> {log}'

# sm1.log
rule vg_index_xg:
    input:
        'output/vg/vg_version/{og}.{identity}.{segment}.vg'
    output:
        'output/vg/index/xg/{og}.{identity}.{segment}.xg'
    log:
        'output/logs/vg/xg_log/vg_index.{og}.{identity}.{segment}.log'
    resources:
        time = '0-0:1:00'
    container: 
        vg
    shell:
        'vg index '
        '-x {output} '
        '{input} '
        '2> {log}'


rule vg_index_pruned:
    input:
        'output/vg/vg_version/{og}.{identity}.{segment}.vg'
    output:
        'output/vg/index/pruned/{og}.{identity}.{segment}.pruned.vg'
    log:
        xg_pruned_log = 'output/logs/vg/xg_pruned_log/vg_index.{og}.{identity}.{segment}.log'
    resources:
        time = '0-0:1:00'
    container: 
        vg
    shell:
        'vg prune -r {input} > {output} '
        '2> {log}'


rule vg_index_gcsa:
    input:
        'output/vg/vg_version/{og}.{identity}.{segment}.vg'
    output:
        'output/vg/index/gcsa/{og}.{identity}.{segment}.gcsa'
    log:
        gcsa_log = 'output/logs/vg/gcsa_log/vg_index.{og}.{identity}.{segment}.log'
    resources:
        time = '0-0:1:00'
    container: 
        vg
    shell:
        'vg index '
        '-x {output.xg} '
        '{input} '
        '2> {log.xg_log}'
        '&& '
        'vg prune -r {input} > {output.xg_pruned} '
        '2> {log.xg_pruned_log}'
        '&& '
        'vg index '
        '-g {output.gcsa} '
        '{output.xg_pruned} '
        '2> {log.gcsa_log}'


        
    #pggb_dir = str('output/pggb/{wildcards.og}.{wildcards.identity}.{wildcards.segment}/{{filename}}.gfa')
    #gfa_file = glob_wildcards(pggb_dir).filename
    #gfa_file_path = 'output/pggb/{wildcards.og}.{wildcards.identity}.{wildcards.segment}/{{filename}}.gfa'
    #return(gfa_file)

# rule vg:
    #input:

    
    #output:
        #directory('output/vg/vg_result/{og}.{identity}.{segment}')
    #log:
        #'output/logs/vg/vg.{og}.{identity}.{segment}.log'
    #resources:
        #time = '0-0:1:00'
    #container: 
        #vg
    #threads:
        #4
    #shell:
        #'vg mpmap '
        #'-n rna'
        #'-t 4 '
        

rule aggregate:
    input:
        aggregate




