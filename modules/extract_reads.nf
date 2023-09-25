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

process extract_kraken_reads {

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

workflow extract_reads {
    take:
        unique_id
        ch_input
    main:
        if ( params.wf_dir ) {
            wf_dir = file("${params.wf_dir}", type: "dir", checkIfExists:true)
            ch_kraken = Channel.fromPath("${wf_dir}/output/kraken/*kreport.txt", type: "file", checkIfExists:true).map { [it.simpleName, it]}
            ch_bracken = Channel.fromPath("${wf_dir}/output/bracken/*kreport_bracken_species.txt", type: "file", checkIfExists:true).map { [it.simpleName, it]}
            ch_assignments1 = Channel.fromPath("${params.wf_dir}/work/*/*/*.kraken2.assignments.tsv", type: "file", checkIfExists:true).map { [it.simpleName, it]}
            ch_assignments2 = Channel.fromPath("${params.wf_dir}/output/kraken_reads_assignments/*kraken2.assignments.tsv", type: "file", checkIfExists:true).map { [it.simpleName.replace("_lineages",""), it]}
            ch_assignments = ch_assignments2.ifEmpty( ch_assignments1 )

            input_taxonomy = file("${params.store_dir}/taxonomy_dir", type: "dir")
            if (input_taxonomy.isEmpty()) {
                taxonomy = unpack_taxonomy("${params.default_taxonomy}")
            } else {
                taxonomy = input_taxonomy
            }

            ch_input
                .join(ch_kraken)
                .map{ barcode, reads, kraken -> "$barcode,$reads,$kraken" }.collectFile(name: "${params.outdir}/${unique_id}/input.csv", newLine: true).set{ input_summary }
            default_input_file = file("resources/input_summary.csv")
            out_summary1 = input_summary.ifEmpty(default_input_file)

            ch_input
                .join(ch_assignments)
                .join(ch_bracken)
                .set{ ch_barcode }

            extract_kraken_reads(ch_barcode, taxonomy)

            check_reads(extract_kraken_reads.out.reads)
            check_reads.out.transpose()
                .map{ it -> [it[1].simpleName.replace("reads_",""), it[0], it[1]] }
                .set{ ch_reads }

            extract_kraken_reads.out.summary.collectFile(name: "${params.outdir}/${unique_id}/extract_summary.json").set{ extract_summary }
            default_extract_file = file("resources/extract_summary.json")
            out_summary2 = extract_summary.ifEmpty(default_extract_file)
        } else {
            exit 1, "Must provide wf_dir with --wf_dir."
        }


    emit:
        reads = ch_reads
        input = out_summary1
        summary = out_summary2
}