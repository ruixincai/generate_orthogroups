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


rule orthofinder:
    input:
        input_files = ref_proteomes.values()
    output:
        output_dir = 'output/orthoresult'
    log:
        'output/logs/orthofinder.log'
    params:
        input_dir = proteomes_directory
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
        '-o {output.output_dir} '
        '-t {threads} '
        '&> {log}'
