process taxid_from_name {
    label 'process_low'

    conda 'bioconda::entrez-direct'
    container "${params.wf.container}@${params.wf.container_sha}"

    input:
        val taxon_name
    output:
        tuple stdout, val(taxon_name)
    """
    export NCBI_API_KEY=${params.ncbi_api_key}
    export EMAIL=${params.email}
    esearch -db taxonomy -query "${taxon_name}" | efetch -format xml  | xtract -pattern Taxon -element TaxId | tr -d '\n'
    """
}

process download_references_by_name {

    label 'process_low'
    maxForks 1

    conda 'bioconda::entrez-direct'
    container "${params.wf.container}@${params.wf.container_sha}"

    input:
        val unique_id
        val taxon_name
    output:
        tuple val(taxon_name), path("taxa_references_${taxon_name}.fa"), path("taxa_refs_${taxon_name}.tsv"), emit: refs
    script:
        """
        export NCBI_API_KEY=${params.ncbi_api_key}
        export EMAIL=${params.email}
        esearch -db genome -query "${taxon_name}"[orgn] | \
            elink -target nuccore | efetch -format docsum | \
            xtract -pattern DocumentSummary -if SourceDb -contains refseq -contains NC_ -element Caption,Title,SourceDb \
            > taxa_refs_${taxon_name}.tsv
        cat taxa_refs_${taxon_name}.tsv | efetch -db nuccore -format fasta > taxa_references_${taxon_name}.fa
        """
}

process download_references_by_taxid {

    label 'process_low'
    maxForks 1
    errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
    maxRetries 2

    conda 'bioconda::entrez-direct'
    container "${params.wf.container}@${params.wf.container_sha}"

    input:
        val taxon_id
    output:
        tuple val(taxon_id), path("taxa_references_${taxon_id}.fa"), path("taxa_refs_${taxon_id}.tsv"), emit: refs
    script:
        """
        export NCBI_API_KEY=${params.ncbi_api_key}
        export EMAIL=${params.email}
        esearch -db genome -query "txid${taxon_id}[Organism:exp]" | \
            elink -target nuccore | efetch -format docsum | \
            xtract -pattern DocumentSummary -if SourceDb -contains refseq -element Caption,Title,SourceDb \
            > taxa_refs_${taxon_id}.tsv
        cat taxa_refs_${taxon_id}.tsv | efetch -db nuccore -format fasta > taxa_references_${taxon_id}.fa
        """
}

process filter_references {

    label 'process_low'

    publishDir path: "${params.outdir}/references", mode: 'copy', pattern: 'references_*.fa'

    errorStrategy {task.exitStatus == 2 ? 'ignore' : 'terminate'}

    conda 'bioconda::biopython=1.78 bioconda::tabix=1.11 bioconda::mappy=2.26'
    container "${params.wf.container}@${params.wf.container_sha}"

    input:
        tuple val(taxon_id), path(reference_fa), path(reference_summary)
    output:
        tuple val(taxon_id), path("references_${taxon_id}.fa"), emit: refs
        tuple val(taxon_id), stdout, emit: counts
    script:
        """
        filter_references.py \
          -r ${reference_fa} -s ${reference_summary} \
          -o "references_${taxon_id}.fa" \
          --segment_sep ${params.reference_segment_sep}

        if [ ! -s "references_${taxon_id}.fa" ]
        then
            echo "No references found for this taxid"
            exit 2
        fi

        grep ">" "references_${taxon_id}.fa" | wc -l | tr -d '\n'
        """
}

workflow download_references {
    main:
        if ( params.taxon_names ) {
            taxa_list = params.taxon_names?.split(',') as List
            ch_taxa = Channel.from(taxa_list)
            taxid_from_name(ch_taxa)
            taxid_from_name.out
                .tap{ ref_pair_ch }
                .map{ it -> it[0] }
                .tap{ taxid_ch }
            download_references_by_taxid(taxid_ch)
            filter_references(download_references_by_taxid.out)

            ref_pair_ch
                .join(filter_references.out.counts)
                .map{ taxid, name, count -> "$name,$taxid,$count" }
                .collectFile(name: "${params.outdir}/references/references.csv", newLine: true)
                .set{ reference_summary }

             default_references_file = file("resources/references_summary.csv")
             out_summary = reference_summary.ifEmpty(default_references_file)

            ch_references = filter_references.out.refs
        } else {
            exit 1, "Must provide a list of taxon names with --taxon_names. These should be a comma separated list, and each taxon item should have quotes to allow effective parsing of spaces"
        }
    emit:
        references = filter_references.out.refs
        summary = out_summary

}


