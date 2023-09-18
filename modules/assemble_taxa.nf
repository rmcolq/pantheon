process subset_references {

    label 'process_low'

    errorStrategy {task.exitStatus == 2 ? 'ignore' : 'terminate'}

    conda 'bioconda::biopython=1.78 bioconda::tabix=1.11 bioconda::mappy=2.26'
    container "${params.wf.container}@${params.wf.container_sha}"

    input:
        tuple val(taxon_id), val(barcode_id), path(reads), path(reference_fa)
    output:
        tuple val(taxon_id), val(barcode_id), path(reads),  path("${barcode_id}_${taxon_id}_reference.fa"), emit: all
        path "summary.csv", emit: summary
    script:
        """
        subset_references.py \
          -r ${reference_fa} \
          -q ${reads} \
          -o "${barcode_id}_${taxon_id}_reference.fa" \
          --min_count ${params.reference_min_count} \
          --taxid ${taxon_id} --barcode ${barcode_id}
        """
}

process check_subset {

    label 'process_low'

    errorStrategy {task.exitStatus == 2 ? 'ignore' : 'terminate'}

    conda 'bioconda::biopython=1.78 bioconda::tabix=1.11 bioconda::mappy=2.26'
    container "${params.wf.container}@${params.wf.container_sha}"

    input:
        tuple val(taxon_id), val(barcode_id), path(reads), path(reference_fa)
    output:
        tuple val(taxon_id), val(barcode_id), path(reads),  path(reference_fa)
    script:
        """
        if [ ! -s "${reference_fa}" ]
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
