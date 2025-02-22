#!/usr/bin/env python3

from pathlib import Path

proteomes_directory = 'data/proteomes'

bbmap = 'docker://quay.io/biocontainers/bbmap:39.01--h5c4e2a8_0'

pggb = 'docker://ghcr.io/pangenome/pggb:20230113201558a9a04c'

vg = 'docker://quay.io/vgteam/vg:v1.47.0'

samtools = 'docker://quay.io/biocontainers/samtools:1.3.1--h0cf4675_11'

rpvg = 'docker://quay.io/jonassibbesen/rpvg'

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
# a combined tabular database include all transcript accession numbers, protein accession numbers,  
# mosquito species is generated
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


# TXT files including transcript IDs for orthogroups
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
# transcript sequences match to transcript IDs
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




# pggb index
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


# construct graphs for each orthogroup from transcript sequences
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


# 10 target orthogroups
list_of_ogs = ['OG0012812','OG0009228','OG0007096','OG0009197','OG0010256','OG0012242','OG0004074','OG0000133', 'OG0008785','OG0001722']


# test pggb on 10 orthogroups with different identity and segment length
rule pggb_test:
    input:
        expand('output/pggb/{og}.{identity}.{segment}',
            og=list_of_ogs,
            identity=[85, 90, 95, 60, 70, 80],
            segment=[100, 300, 3000])  # minimum segment length is required to be >= 100 bp


# get GFA graphs from pggb
# 10 * 6 * 3 = 180 GFA graphs
def gfa_input(wildcards):
    checkpoint_output = checkpoints.pggb.get(**wildcards).output['output_dir']
    gfa_file_path = f'{checkpoint_output}/{{filename}}.gfa'
    gfa_file = glob_wildcards(gfa_file_path).filename
    return(expand(gfa_file_path,filename = gfa_file))




# test vg on 10 orthogroups with different identity and segment length
rule vg_test: 
    input:
        expand('output/vg/vg_bam/SRR1645076/{identity}.{segment}.map.bam',
            identity=[85, 90, 95, 60, 70, 80],
            segment=[100, 300, 3000])


# analysis of GFA graphs
rule vg_stat:
    input:
        gfa_input
    output:
        'output/matrix/{og}.{identity}.{segment}.txt'
    log:
        'output/logs/matrix/{og}.{identity}.{segment}.log'
    resources:
        time = '0-0:1:00'
    container: 
        vg
    shell:
        'vg stats '
        '-z '  # size of graph
        '-N '  # number of nodes in graph
        '-E '  # number of edges in graph
        '-s '  # describe subgraphs of graph
        '{input} '
        '> {output} '
        '2> {log}'


# convert GFA files to vg files
rule gfa2vg:
    input:
        gfa_input
    output:
        'output/vg/vg_format/{og}.{identity}.{segment}.vg'
    log:
        'output/logs/vg/vg_format/{og}.{identity}.{segment}.log'
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


# merge vg files of 10 orthogroups into a single vg file 
# with different identity and segment length
# 6 * 3 = 18 merged graphs
rule merge_vg:
    input:
        vg_files = expand('output/vg/vg_format/{og}.{{identity}}.{{segment}}.vg', 
            og=list_of_ogs)
    output:
        merged_graph = 'output/vg/merged_graph/{identity}.{segment}.merged_graph.vg'
    log:
        'output/logs/vg/merged_graph/{identity}.{segment}.merged_graph.log'
    resources:
        time = '0-0:1:00'
    container: 
        vg
    shell:
        'vg combine '
        '{input} '
        '> {output} '
        '2> {log}'



rule vg_index_xg:
    input:
        'output/vg/merged_graph/{identity}.{segment}.merged_graph.vg'
    output:
        'output/vg/index/xg/{identity}.{segment}.xg'
    log:
        'output/logs/vg/xg_log/{identity}.{segment}.log'
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
        'output/vg/merged_graph/{identity}.{segment}.merged_graph.vg'
    output:
        'output/vg/index/pruned/{identity}.{segment}.pruned.vg'
    log:
        'output/logs/vg/pruned_log/{identity}.{segment}.log'
    resources:
        time = '0-0:1:00'
    container: 
        vg
    shell:
        'vg prune -r {input} > {output} '
        '2> {log}'


rule vg_index_mod_pruned:
    input:
        'output/vg/index/pruned/{identity}.{segment}.pruned.vg'
    output:
        'output/vg/index/mod_pruned/{identity}.{segment}.mod.pruned.vg'
    log:
        'output/logs/vg/mod_pruned_log/{identity}.{segment}.log'
    resources:
        time = '0-0:1:00'
    container: 
        vg
    shell:
        'vg mod -X 256 {input} > {output} '
        '2> {log}'


rule vg_index_gcsa:
    input:
        'output/vg/index/mod_pruned/{identity}.{segment}.mod.pruned.vg'
    output:
        'output/vg/index/gcsa/{identity}.{segment}.gcsa'
    log:
        'output/logs/vg/gcsa_log/{identity}.{segment}.log'
    resources:
        time = '0-0:1:00'
    container: 
        vg
    shell:
        'vg index '
        '-g {output} '
        '{input} '
        '2> {log}'


# rule vg_index_snarls:
    # input:
        # 'output/vg/index/xg/{og}.{identity}.{segment}.xg'
    # output:
        # 'output/vg/index/snarls/{og}.{identity}.{segment}.trivial.snarls' # sm5
    # log:
        # 'output/logs/vg/snarls_log/vg_index.{og}.{identity}.{segment}.log'
    # resources:
        # time = '0-0:1:00'
    # container: 
        # vg
    # shell:
        # 'vg snarls '
        # '-T {input} '
        # '> {output} '
        # '2> {log}'
# -s 
# invalid option -- 's' in vg index
# create the distance index no longer need the snarls files


rule vg_index_dist:
    input:
        xg = 'output/vg/index/xg/{identity}.{segment}.xg'
        # snarl = 'output/vg/index/snarls/{og}.{identity}.{segment}.trivial.snarls'
    output:
        'output/vg/index/dist/{identity}.{segment}.dist'
    log:
        'output/logs/vg/dist_log/{identity}.{segment}.log'
    resources:
        time = '0-0:1:00'
    container: 
        vg
    shell:
        'vg index '
        '-x {input.xg} '
        # '-s {input.snarl} '
        '-j {output} '
        '2> {log}'

        
    #pggb_dir = str('output/pggb/{wildcards.og}.{wildcards.identity}.{wildcards.segment}/{{filename}}.gfa')
    #gfa_file = glob_wildcards(pggb_dir).filename
    #gfa_file_path = 'output/pggb/{wildcards.og}.{wildcards.identity}.{wildcards.segment}/{{filename}}.gfa'
    #return(gfa_file)

rule vg_mpmap:
    input:
        xg = 'output/vg/autoindex/{identity}.{segment}.merged_graph.xg',
        gcsa = 'output/vg/autoindex/{identity}.{segment}.merged_graph.gcsa',
        dist = 'output/vg/autoindex/{identity}.{segment}.merged_graph.dist',
        fastq = 'data/RNA_seq/SRR520427.fq'
    output:
        'output/vg/vg_mpmap/{identity}.{segment}.gamp'
    log:
        'output/logs/vg/vg_mpmap/{identity}.{segment}.log'
    resources:
        time = '0-0:10:00'
    container: 
        vg
    shell:
        'vg mpmap '
        '-x {input.xg} '
        '-g {input.gcsa} '
        '-d {input.dist} ' 
        '-n RNA '
        '-f {input.fastq} '
        '> {output} '
        '2> {log}'



# vg map
rule vg_sim:
    input:
        xg = 'output/vg/index/xg/{identity}.{segment}.xg'
    output:
        'output/vg/vg_sim/{identity}.{segment}.sim.txt'
    log:
        'output/logs/vg/vg_sim/{identity}.{segment}.log'
    resources:
        time = '0-0:1:00'
    container: 
        vg
    shell:
        'vg sim '
        '-n 1000 '
        '-l 150 '
        '-x {input.xg} '
        '> {output} '
        '2> {log}'


# convert merged graphs to merged GFA graphs
rule vg_vg2gfa:
    input:
        merged_graph = 'output/vg/merged_graph/{identity}.{segment}.merged_graph.vg'
    output:
        'output/vg/merged_graph_gfa/{identity}.{segment}.merged_graph.gfa'
    log:
        'output/logs/vg/merged_graph_gfa/{identity}.{segment}.merged_graph.log'
    resources:
        time = '0-0:5:00'
    container: 
        vg
    shell:
        'vg view '
        '{input.merged_graph} > {output} '
        '2> {log}'


# analysis of merged GFA graphs
rule vg_stat_merged:
    input:
        'output/vg/merged_graph_gfa/{identity}.{segment}.merged_graph.gfa'
    output:
        'output/matrix_merged/{identity}.{segment}.txt'
    log:
        'output/logs/matrix_merged/{identity}.{segment}.log'
    resources:
        time = '0-0:1:00'
    container: 
        vg
    shell:
        'vg stats '
        '-N '  # number of nodes in graph
        '-E '  # number of edges in graph
        '-s '  # describe subgraphs of graph
        '{input} '
        '> {output} '
        '2> {log}'


# index merged GFA graphs
rule vg_autoindex:
    input:
        'output/vg/merged_graph_gfa/{identity}.{segment}.merged_graph.gfa'
    output:
        multiext('output/vg/autoindex/{identity}.{segment}.merged_graph', '.xg', '.gcsa', '.gcsa.lcp')
    params:
        prefix = 'output/vg/autoindex/{identity}.{segment}.merged_graph',
        workflow ='map'
    log:
        'output/logs/vg/autoindex/{identity}.{segment}.vg_autoindex.log'
    resources:
        time = '0-0:5:00'
    container: 
        vg
    shell:
        'vg autoindex '
        '--workflow {params.workflow} '
        '--prefix {params.prefix} '
        '--gfa {input} '
        '2> {log}'


# vg map
rule vg_map:
    input:
        rules.vg_autoindex.output,
        fastq = 'data/RNA_seq/SRR1645076.fq'
    output:
        'output/vg/vg_map/SRR1645076/{identity}.{segment}.map.gam'
    params:
        prefix = 'output/vg/autoindex/{identity}.{segment}.merged_graph'
    log:
        'output/logs/vg/vg_map/{identity}.{segment}.vg_map.log'
    threads:    
        lambda wildcards, attempt: 10 * attempt
    resources:  
        time = lambda wildcards, attempt: 10 * attempt
    container: 
        vg
    shell:
        'vg map '
        '-t {threads} '
        '--log-time '
        '-d {params.prefix} '
        '-f {input.fastq} '
        '--interleaved '
        '> {output} '
        '2> {log}'


# gbwt index
rule vg_gbwt:
    input:
        xg = 'output/vg/autoindex/{identity}.{segment}.merged_graph.xg',
        gfa = 'output/vg/merged_graph_gfa/{identity}.{segment}.merged_graph.gfa'
    output:
        'output/vg/gbwt/{identity}.{segment}.gbwt'
    log:
        'output/logs/vg/gbwt/{identity}.{segment}.gbwt.log'
    resources:
        time = '0-0:5:00'
    container: 
        vg
    shell:
        'vg gbwt '
        '-G {input.gfa} '
        '-o {output} '
        '2> {log}'


# analysis of gam graphs
rule vg_stat_gam:
    input:
        'output/vg/vg_map/{identity}.{segment}.map.gam'
    output:
        'output/matrix_gam/{identity}.{segment}.txt'
    log:
        'output/logs/matrix_gam/{identity}.{segment}.log'
    resources:
        time = '0-0:5:00'
    container: 
        vg
    shell:
        'vg stats '
        '-a '  # number of nodes in graph
        '{input} '
        '> {output} '
        '2> {log}'


# infers path posterior probabilities and abundances from variation graph read alignments
# rule rpvg:
    # input:
        # xg = 'output/vg/autoindex/{identity}.{segment}.merged_graph.xg',
        # gbwt = 'output/vg/gbwt/{identity}.{segment}.gbwt',
        # gam = 'output/vg/vg_map/{identity}.{segment}.map.gam'
    # output:
        # txt = 'output/rpvg/{identity}.{segment}.txt'
        # merged = 'output/vg/merged_graph/{identity}.{segment}.merged_graph.vg'
    # log:
        # 'output/logs/rpvg/{identity}.{segment}.rpvg.log'
    # resources:
        # time = '0-0:5:00'
    # container: 
        # rpvg
    # shell:
        # 'rpvg '
        # '-g {input.xg} '
        # '-p {input.gbwt} '
        # '-a {input.gam} '
        # '-i '
        # '-o {output} '
        # '2> {log}'


rule rpvg:
    input:
        xg = 'output/vg/autoindex/{identity}.{segment}.merged_graph.xg',
        gbwt = 'output/vg/gbwt/{identity}.{segment}.gbwt',
        gam = 'output/vg/vg_map/{identity}.{segment}.map.gam'
    output:
        'output/rpvg/{identity}.{segment}.txt'
    log:
        'output/logs/rpvg/{identity}.{segment}.rpvg.log'
    threads:    
        lambda wildcards, attempt: 10 * attempt
    resources:  
        time = lambda wildcards, attempt: 10 * attempt
    container: 
        rpvg
    shell:
        'rpvg '
        '-g {input.xg} '
        '-p {input.gbwt} '
        '-a {input.gam} '
        '-i transcripts '
        '-o {output} '
        '2> {log}'



# convert gam files to bam files
# rule vg_surject:
    # input:
        # xg = 'output/vg/autoindex/{identity}.{segment}.merged_graph.xg',
        # gam = 'output/vg/vg_map/{identity}.{segment}.map.gam'
    # output:
        # 'output/vg/vg_surject/{identity}.{segment}.map.bam'
    # log:
        # 'output/logs/vg/vg_surject/{identity}.{segment}.vg_surject.log'
    # threads:    
        # lambda wildcards, attempt: 10 * attempt
    # resources:  
        # time = lambda wildcards, attempt: 10 * attempt
    # container: 
        # vg
    # shell:
        # 'vg surject '
        # '-x {input.xg} '
        # '-b {input.gam} '
        # '> {output} '
        # '2> {log}'


# vg view -a 60.100.map.gam


# convert gam files to bam files
rule vg_bam:
	input:
		xg = 'output/vg/autoindex/{identity}.{segment}.merged_graph.xg',
		gam = 'output/vg/vg_map/SRR1645076/{identity}.{segment}.map.gam'
	output:
		'output/vg/vg_bam/SRR1645076/{identity}.{segment}.map.bam'
	log:
		'output/logs/vg/vg_bam/SRR1645076/{identity}.{segment}.log'
	threads:
		lambda wildcards, attempt: 10 * attempt
	resources:
		time = lambda wildcards, attempt: 10 * attempt
	container:
		vg
	shell:
		'vg surject '
		'-x {input.xg} '
		'-b {input.gam} '
		'> {output} '
		'2> {log}'


# get sorted bam files
rule sortedbam:
	input:
		'output/vg/vg_bam/{identity}.{segment}.map.bam'
	output:
		'output/sortedbam/{identity}.{segment}.sorted.bam'
	log:
		'output/logs/sortedbam/{identity}.{segment}.log'
	threads:
		lambda wildcards, attempt: 10 * attempt
	resources:
		time = lambda wildcards, attempt: 10 * attempt
	container:
		samtools
	shell:
		'samtools sort '
		'{input} '
		'-o {output} '
		'2> {log}'


# count the alignments
# rule samtools:
	# input:
		# 'output/vg/vg_sortedbam/{identity}.{segment}.map.bam'
	# output:
		# 'output/samtools/{identity}.{segment}.txt'
	# log:
		# 'output/logs/samtools/{identity}.{segment}.log'
	# resources:
		# time = '0-0:5:00'
	# container: 
		# samtools
	# shell:
		# 'samtools view '
		# '-c '
		# '-F 2 '
    	# '{input} '
    	# '> {output} '
    	# '2> {log}'
        

rule aggregate:
    input:
        aggregate




