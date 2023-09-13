process move_or_compress {

    label "process_low"

    conda "bioconda::tabix==v1.11"
    container "biocontainers/tabix:1.11--hdfd78af_0"
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

    errorStrategy {task.exitStatus == 2 ? 'ignore' : 'terminate'}

    publishDir path: "${params.outdir}/${unique_id}/reads_by_taxa", mode: 'copy'

    conda 'bioconda::biopython=1.78 bioconda::tabix=1.11'
    container "${params.wf.container}@${params.wf.container_sha}"

    input:
        tuple val(unique_id), path(fastq), path(kraken_assignments), path(bracken_report)
        path taxonomy_dir
    output:
        path "reads.*.f*.gz", emit: reads
        path "reads_summary.json", emit: summary
    script:
        """
        extract_kraken_reads.py \
            -s ${fastq} \
            -k ${kraken_assignments} \
            -r ${bracken_report} \
            -t ${taxonomy_dir} \
            -p reads \
            --include_children \
            --max_human ${params.max_human_reads_before_rejection} \
            --min_count_descendants ${params.assembly_min_reads} \
            --rank ${params.assembly_rank} \
            --min_percent ${params.assembly_min_percent}

        PATTERN=(reads.*.f*)
        if [ -f \${PATTERN[0]} ]; then
            for f in \$(ls reads.*.f*)
              do
                bgzip --threads $task.cpus \$f
              done
        else
            exit 2
        fi
        """
}
