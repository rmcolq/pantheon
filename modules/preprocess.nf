process move_or_compress {

    label "process_low"

    conda "bioconda::tabix==v1.11"
    container "quay.io/biocontainers/tabix:1.11--hdfd78af_0"
    cpus 1

    input:
        // don't stage `input` with a literal because we check the file extension
        tuple val(barcode), path(input)
    output:
        tuple val("${barcode}"), path("${barcode}.all.fastq.gz")
    script:
        """
        for file in $input
        do
            if [[ "\$file" == *.gz ]]
            then
                cat "\$file" | gunzip | bgzip -@ $task.cpus >> "${barcode}.all.fastq.gz"
            else
                cat "\$file" | bgzip -@ $task.cpus >> "${barcode}.all.fastq.gz"
            fi
        done
        """
}

process unpack_taxonomy {
    label "process_single"
    storeDir "${params.store_dir}"
    input:
        path taxonomy
    output:
        path "taxonomy_dir"
    """
    if [[ "${taxonomy}" == *.tar.gz ]]
    then
        mkdir taxonomy_dir
        tar xf "${taxonomy}" -C taxonomy_dir
    elif [ -d "${taxonomy}" ]
    then
        mv "${taxonomy}" taxonomy_dir
    else
        echo "Error: taxonomy is neither .tar.gz nor a dir"
        echo "Exiting".
        exit 1
    fi
    """
}

process extract_reads {

    label 'process_low'

    publishDir path: "${params.outdir}/${unique_id}/reads_by_taxa", mode: 'copy', pattern: "reads_summary.json"

    conda 'bioconda::biopython=1.78 bioconda::tabix=1.11'
    container "${params.wf.container}@${params.wf.container_sha}"

    input:
        tuple val(unique_id), path(fastq), path(kraken_assignments), path(bracken_report)
        path taxonomy_dir
    output:
        tuple val(unique_id), path("reads_*.f*"), path("reads_human.txt"), emit: reads
        path "${unique_id}_summary.json", emit: summary
    script:
        if ( params.taxon_names )
            taxid = "--taxid \"${params.taxon_names}\""
        else
            taxid = ""
        """
        extract_kraken_reads.py \
            -s ${fastq} \
            -k ${kraken_assignments} \
            -r ${bracken_report} \
            -t ${taxonomy_dir} \
            -p reads \
            --include_children \
            --min_count_descendants ${params.extract_min_reads} \
            --rank ${params.extract_rank} \
            --min_percent ${params.extract_min_percent} ${taxid}

        file1=`cat reads_summary.json`
        echo "{"'"${unique_id}"'": "\$file1"}" >> "${unique_id}_summary.json"
        """
}

process check_reads {

    label 'process_low'

    errorStrategy {task.exitStatus == 2 ? 'ignore' : 'terminate'}

    publishDir path: "${params.outdir}/${unique_id}/reads_by_taxa", mode: 'copy'

    conda 'bioconda::biopython=1.78 bioconda::tabix=1.11'
    container "${params.wf.container}@${params.wf.container_sha}"

    input:
        tuple val(unique_id), path(reads), path(human)
    output:
        tuple val(unique_id), path("reads_*.f*.gz")
    script:
        """
        PATTERN=(reads_*.f*)
        if [ -f \${PATTERN[0]} ]; then
            for f in \$(ls reads_*.f*)
              do
                if [ ! -s \$f ]
                  then
                    rm \$f
                  else
                    bgzip --threads $task.cpus \$f
                fi
              done
        fi

        PATTERN=(reads_*.f*.gz)
        if [ ! -f \${PATTERN[0]} ]; then
            echo "Found no output files - maybe there weren't any for this sample"
            exit 2
        fi

        a=\$(< ${human})
        if [ "\$a" -gt ${params.max_human_reads_before_rejection} ]; then
            echo "Too many human reads"
            exit 2
        fi
        """
}

process taxid_from_name {
    label 'process_low'

    conda 'bioconda::entrez-direct'
    container "${params.wf.container}@${params.wf.container_sha}"

    input:
        val taxon_name
    output:
        tuple val(taxon_name), stdout
    """
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
        tuple val(taxon_id), path("references_${taxon_id}.fa")
    script:
        """
        filter_references.py \
          -r ${reference_fa} -s ${reference_summary} \
          -o "references_${taxon_id}.fa" \
          --segment_sep ${params.reference_segment_sep}

        echo

        if [ ! -s "references_${taxon_id}.fa" ]
        then
            echo "No references found for this taxid"
            exit 2
        fi
        """
}


