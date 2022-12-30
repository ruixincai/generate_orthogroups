#!/usr/bin/env python3

from pathlib import Path

proteomes_directory = 'data/proteomes'

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

# leave this as the first rule
rule target:
    input:
        'output/orthoresult/Results/Orthogroups/Orthogroups.tsv',
        expand('output/separated_transcript_protein_tables/{species}.csv',
               species=list(ref_gff.keys()))
        'output/combined_table.csv',
        'output/ortho_match.csv',
        directory('output/separated_orthogroups')


rule orthofinder:
    input:
        input_files = ref_proteomes.values() # Path()
    output:
        orthoresult = 'output/orthoresult/Results/Orthogroups/Orthogroups.tsv'
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


rule ortho_match:
    input:
        input_orthofinder = 'output/orthoresult/Results/Orthogroups/Orthogroups.tsv',
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


rule separated_orthogroups:
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

