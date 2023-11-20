// workflow to extract reads and create assemblies from a wf-metagenomic or scylla run
include { download_references } from '../modules/download_references'
include { move_or_compress } from '../modules/preprocess'
include { extract_reads } from '../modules/extract_reads'
include { assemble } from '../modules/assemble_taxa'
include { generate_heatmap_report;generate_assembly_report } from '../modules/generate_report'

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

        if (params.heatmap){
            generate_heatmap_report(unique_id, index_csv)
        }

        if (params.assemble){
            download_references()
            ch_references = download_references.out.references

            index.map{ row-> tuple(row.sample_id, get_fq_files_in_dir(file(row.read_path))) }
                .transpose()
                .map{ it -> [it[1].simpleName.replace("_1","").replace("_2",""), it[0], it[1]] }
                .combine(ch_references, by: 0)
                .groupTuple(by:[0,1])
                .branch {
                    paired: it[2].size() == 2
                    single: it[2].size() == 1
                }
                .set { read_tuples }

            assemble(unique_id, read_tuples.single)
            default_extract_file = file("resources/extract_summary.json")
            generate_assembly_report(unique_id, index_csv, download_references.out.summary, default_extract_file, assemble.out.summary)

        }



}
