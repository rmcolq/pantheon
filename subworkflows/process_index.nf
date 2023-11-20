// workflow to extract reads and create assemblies from a wf-metagenomic or scylla run
include { download_references } from '../modules/download_references'
include { move_or_compress } from '../modules/preprocess'
include { extract_reads } from '../modules/extract_reads'
include { assemble } from '../modules/assemble_taxa'
include { generate_heatmap_report } from '../modules/generate_report'

EXTENSIONS = ["fastq", "fastq.gz", "fq", "fq.gz"]

ArrayList get_fq_files_in_dir(Path dir) {
    return EXTENSIONS.collect { file(dir.resolve("*.$it"), type: "file") } .flatten()
}


workflow process_index {
    take:
        unique_id
    main:
        index_csv = Channel.fromPath(params.index_csv)
        index = index_csv.splitCsv(header:true)

        index.map{ row-> file(row.kraken_report) }.collect().set{ kreport_ch }
        generate_heatmap_report(unique_id, index_csv, kreport_ch)

        index.map{ row-> tuple(row.sample_id, get_fq_files_in_dir(file(row.read_path))) }
             .transpose()
             .map{ it -> [it[1].simpleName.replace("reads_",""), it[0], it[1]] }
             .view()


}
