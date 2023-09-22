// workflow to extract reads and create assemblies from a wf-metagenomic or scylla run
include { download_references } from '../modules/download_references'
include { move_or_compress } from '../modules/preprocess'
include { extract_reads } from '../modules/extract_reads'
include { assemble } from '../modules/assemble_taxa'
include { generate_assembly_report } from '../modules/generate_report'

EXTENSIONS = ["fastq", "fastq.gz", "fq", "fq.gz"]

ArrayList get_fq_files_in_dir(Path dir) {
    return EXTENSIONS.collect { file(dir.resolve("*.$it"), type: "file") } .flatten()
}


workflow process_run {
    take:
        unique_id
    main:
        run_dir = file("${params.run_dir}", type: "dir", checkIfExists:true)
        barcode_input = Channel.fromPath("${run_dir}/*", type: "dir", checkIfExists:true, maxDepth:1).map { [it.baseName, get_fq_files_in_dir(it)]}
        ch_input = move_or_compress(barcode_input)

        download_references()
        ch_references = download_references.out.references

        extract_reads(unique_id, ch_input)
        ch_reads = extract_reads.out.reads

        ch_read_ref_pairs = ch_reads.combine(ch_references, by: 0)
        assemble(unique_id, ch_read_ref_pairs)

        generate_assembly_report(unique_id, extract_reads.out.input, download_references.out.summary, extract_reads.out.summary, assemble.out.summary)

}