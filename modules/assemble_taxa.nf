process filter_references {

    label 'process_low'

    errorStrategy {task.exitStatus == 2 ? 'ignore' : 'terminate'}

    conda 'bioconda::biopython=1.78 bioconda::tabix=1.11 bioconda::mappy=2.26'
    container "${params.wf.container}@${params.wf.container_sha}"

    input:
        tuple val(taxon_id), val(barcode_id), path(reads), path(reference_fa), path(reference_summary)
    output:
        tuple val(taxon_id), val(barcode_id), path(reads),  path("filtered_references.fa")
    script:
        """
        infer_references.py \
          -r ${reference_fa} -s ${reference_summary} \
          -q ${reads} \
          --min_count ${params.reference_min_count} \
          --segment_sep ${params.reference_segment_sep}

        if [ ! -s "filtered_references.fa" ]
        then
            echo "No references for this taxid had sufficient reads to continue"
            exit 2
        fi
        """
}

process medaka_consensus {

    label 'process_low'

    publishDir path: "${params.outdir}/${barcode_id}/assemblies", mode: 'copy', pattern: 'assembly_*.fa'

    conda 'bioconda::medaka=1.8.0'
    container "quay.io/biocontainers/medaka@sha256:777dc707f163a3b56b5e18f066cfac3c3fc96ef029332e7c0bdd1db81ada412f"

    input:
        tuple val(taxon_id), val(barcode_id), path(reads), path(reference)
    output:
        tuple val(taxon_id), val(barcode_id), path(reads),  path("assembly_${taxon_id}.fa")
    script:
        """
        medaka_haploid_variant -i ${reads} \
            -r ${reference} \
            -o "medaka_output" \
            -f -x && \
        medaka stitch "medaka_output/consensus_probs.hdf" ${reference} "medaka_output/consensus.fasta"
        mv medaka_output/consensus.fasta "assembly_${taxon_id}.fa"
        """
}
