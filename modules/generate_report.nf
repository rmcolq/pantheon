process generate_assembly_report {

    label 'process_low'

    publishDir path: "${params.outdir}/", mode: 'copy'

    conda 'bioconda::biopython=1.78 anaconda::Mako=1.2.3'
    container "${params.wf.container}@${params.wf.container_sha}"

    input:
        val unique_id
        path input_summary
        path reference_summary
        path extract_summary
        path assembly_summary
    output:
        path "${unique_id}_assembly_report.html"

    script:
        """
        assembly_report.py --input ${input_summary} \
          --reference_summary ${reference_summary} \
          --extract_summary ${extract_summary} \
          --assembly_summary ${assembly_summary} \
          --min_read_count ${params.reference_min_count} \
          --run ${unique_id} \
          --relative_directory ${launchDir} \
          --version ${workflow.manifest.version}
        """
}

process generate_heatmap_report {

    label 'process_low'

    publishDir path: "${params.outdir}/", mode: 'copy'

    conda 'bioconda::biopython=1.78 anaconda::Mako=1.2.3'
    container "${params.wf.container}@${params.wf.container_sha}"

    input:
        val unique_id
        path index_csv
    output:
        path "${unique_id}_heatmap_report.html"

    script:
        """
        heatmap_report2.py --input ${index_csv} \
          --prefix ${unique_id} \
          --relative_directory ${launchDir} \
          --version ${workflow.manifest.version}
        """
}