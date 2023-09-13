process download_references_by_name {

    label 'process_low'

    publishDir path: "${params.outdir}/${unique_id}/references", mode: 'copy'

    conda 'bioconda::biopython=1.78 bioconda::tabix=1.11'
    container "${params.wf.container}@${params.wf.container_sha}"

    input:
        val unique_id
        val taxon_name
    output:
        path "references_${taxon_name}.fa", emit: refs
    script:
        """
        esearch -db genome -query "${taxon_name}"[orgn] | \
            elink -target nuccore | efetch -format docsum | \
            xtract -pattern DocumentSummary -if SourceDb -contains refseq -element Caption,Title,SourceDb | \
            efetch -db nuccore -format fasta > references_${taxon_name}.fa
        """
}

process download_references_by_taxid {

    label 'process_low'

    publishDir path: "${params.outdir}/${unique_id}/references", mode: 'copy'

    conda 'bioconda::biopython=1.78 bioconda::tabix=1.11'
    container "${params.wf.container}@${params.wf.container_sha}"

    input:
        val unique_id
        val taxon_id
    output:
        path "references_${taxon_id}.fa", emit: refs
    script:
        """
        esearch -db genome -query "txid${taxon_id}[Organism:exp]" | \
            elink -target nuccore | efetch -format docsum | \
            xtract -pattern DocumentSummary -if SourceDb -contains refseq -element Caption,Title,SourceDb | \
            efetch -db nuccore -format fasta > references_${taxon_id}.fa
        """
}

