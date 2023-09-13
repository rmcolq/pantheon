include { process_run } from './subworkflows/process_run'

workflow {
    unique_id = "${params.unique_id}"

    if (unique_id == "null") {
        if (params.run_dir) {
            run_dir = file(params.run_dir, type: "dir", checkIfExists:true)
            unique_id = "${run_dir.simpleName}"
        } else {
            exit 1, "Note run_dir AND wf_dir need to be provided -- aborting"
        }
    }

    process_run(unique_id)
}