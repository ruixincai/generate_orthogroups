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

rule orthofinder:
    input:
        input_files = ref_proteomes.values()
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
	'docker://ghcr.io/tomharrop/r-containers:latest'
    script:
	'separated_transcript_protein.R'
