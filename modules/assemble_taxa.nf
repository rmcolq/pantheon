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

process map_to_consensus {

    label 'process_low'

    conda 'bioconda::minimmap2=v2.26 bioconda::samtools=v1.17'
    container "${params.wf.container}@${params.wf.container_sha}"

    input:
        tuple val(taxon_id), val(barcode_id), path(reads), path(assembly)
    output:
        tuple val(taxon_id), val(barcode_id), path(reads), path(assembly), path("out.bam"), path("out.bam.csi")
    script:
        """
        minimap2 -x map-ont -a ${assembly} ${reads} | samtools sort -o out.bam --write-index -
        """
}

process define_mask {

    label 'process_low'

    container "${params.wf.container}@${params.wf.container_sha}"

    input:
        tuple val(taxon_id), val(barcode_id), path(reads), path(assembly), path(bam), path(index)
    output:
        tuple val(taxon_id), val(barcode_id), path(reads), path(assembly), path("mask.bed.tsv")
    script:
        """
        str=\$(head -n1 ${assembly})
        ref="\${str: 1}"
        maskara -d ${params.mask_depth} -r \$ref -o mask.bed ${bam}
        """
}

process apply_mask {

    label 'process_low'

    publishDir path: "${params.outdir}/${barcode_id}/assemblies", mode: 'copy', pattern: 'assembly_*.fa'

    conda 'bioconda::bedtools=v2.31.0'
    container "${params.wf.container}@${params.wf.container_sha}"

    input:
        tuple val(taxon_id), val(barcode_id), path(reads), path(assembly), path(mask)
    output:
        tuple val(taxon_id), val(barcode_id), path(reads), path("assembly_${taxon_id}.masked.fa")
    script:
        """
        bedtools maskfasta -fi ${assembly} -bed ${mask} -fo "assembly_${taxon_id}.masked.fa"
        """
}

workflow assemble {
    take:
        unique_id
        ch_assembly
    main:
        subset_references(ch_assembly)
        check_subset(subset_references.out.all)
        medaka_consensus(check_subset.out)
        map_to_consensus(medaka_consensus.out)
        define_mask(map_to_consensus.out)
        apply_mask(define_mask.out)
        subset_references.out.summary.collectFile(name: "${params.outdir}/${unique_id}/reference_by_barcode.csv", keepHeader: true, skip:1).set{ assembly_summary }
        default_assembly_file = file("resources/assembly_summary.csv")
        out_summary = assembly_summary.ifEmpty(default_assembly_file)
    emit:
        completed = apply_mask.out
        summary = out_summary
    }