include { process_run } from './subworkflows/process_run'
include { process_index } from './subworkflows/process_index'

workflow {
    unique_id = "${params.unique_id}"

    if (unique_id == "null") {
        if (params.run_dir) {
            run_dir = file(params.run_dir, type: "dir", checkIfExists:true)
            unique_id = "${run_dir.simpleName}"
        } else if (params.index_csv) {
            index = file(params.index_csv, type: "file", checkIfExists:true)
            unique_id = "${index.simpleName}"
        } else {
            exit 1, "Note run_dir AND wf_dir need to be provided -- aborting"
        }
    }

    if (params.run_dir)
        process_run(unique_id)
    else if (params.index_csv)
        process_index(unique_id)
}