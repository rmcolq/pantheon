// workflow to extract reads and create assemblies from a wf-metagenomic or scylla run
include { move_or_compress; unpack_taxonomy; extract_reads; check_reads; taxid_from_name; download_references_by_taxid; filter_references } from '../modules/preprocess'
include { subset_references; check_subset; medaka_consensus } from '../modules/assemble_taxa'

EXTENSIONS = ["fastq", "fastq.gz", "fq", "fq.gz"]

ArrayList get_fq_files_in_dir(Path dir) {
    return EXTENSIONS.collect { file(dir.resolve("*.$it"), type: "file") } .flatten()
}



workflow process_barcode {
    take:
        barcode_id
        barcode_fq
        taxonomy
    main:
        if ( params.wf_dir ) {
            wf_dir = file("${params.wf_dir}", type: "dir", checkIfExists:true)
            bracken_report = file("${params.wf_dir}/output/bracken/${barcode_id}.kreport_bracken_species.txt", type: "file", checkIfExists:true)
            kraken_assignments = file("${params.wf_dir}/work/*/*/${barcode_id}.kraken2.assignments.tsv", type: "file", checkIfExists:true)
            wf_report = file("${params.wf_dir}/output/wf-metagenomics-report.html", type: "file", checkIfExists:true)

            extract_reads(barcode_id, barcode_fq, kraken_assignments, bracken_report, taxonomy)
        } else {
            exit 1, "Note wf_dir needs to be provided -- aborting"
        }
    emit:
        barcode_summary = extract_reads.out.summary
        barcode_kreport = bracken_report
        barcode_wfreport = wf_report
        barcode_id = barcode_id
}

workflow process_run {
    take:
        unique_id
    main:
        run_dir = file("${params.run_dir}", type: "dir", checkIfExists:true)
        barcode_input = Channel.fromPath("${run_dir}/*", type: "dir", checkIfExists:true, maxDepth:1).map { [it.baseName, get_fq_files_in_dir(it)]}
        ch_input = move_or_compress(barcode_input)

        if ( params.taxon_names ) {
            taxa_list = params.taxon_names?.split(',') as List
            ch_taxa = Channel.from(taxa_list)
            taxid_from_name(ch_taxa)
            taxid_from_name.out.map{ it -> it[1] }.set{ taxid_ch }
            taxid_from_name.out.map{ name, taxid -> "$name,$taxid" }.collectFile(name: "${params.outdir}/references/references.csv", newLine: true)

        } else {
            exit 1, "Must provide a list of taxon names with --taxon_names. These should be a comma separated list, and each taxon item should have quotes to allow effective parsing of spaces"
        }
        download_references_by_taxid(taxid_ch)
        ch_references = filter_references(download_references_by_taxid.out)

        if ( params.wf_dir ) {
            wf_dir = file("${params.wf_dir}", type: "dir", checkIfExists:true)
            ch_bracken = Channel.fromPath("${wf_dir}/output/bracken/*.kreport_bracken_species.txt", type: "file", checkIfExists:true).map { [it.simpleName, it]}
            ch_assignments = Channel.fromPath("${params.wf_dir}/work/*/*/*.kraken2.assignments.tsv", type: "file", checkIfExists:true).map { [it.simpleName, it]}

            input_taxonomy = file("${params.store_dir}/taxonomy_dir", type: "dir")
            if (input_taxonomy.isEmpty()) {
               taxonomy = unpack_taxonomy("${params.default_taxonomy}")
            } else {
               taxonomy = input_taxonomy
            }

            ch_input
                .join(ch_assignments)
                .join(ch_bracken)
                .set{ ch_barcode }

            extract_reads(ch_barcode, taxonomy)
            check_reads(extract_reads.out.reads)
            check_reads.out.transpose()
              .map{ it -> [it[1].simpleName.replace("reads_",""), it[0], it[1]] }
              .set{ ch_reads }

            ch_assembly = ch_reads.combine(ch_references, by: 0)
            subset_references(ch_assembly)
            check_subset(subset_references.out.all)
            medaka_consensus(check_subset.out)
            subset_references.out.summary.collectFile(name: "${params.outdir}/${unique_id}/reference_by_barcode.csv", keepHeader: true, skip:1)
        }
}