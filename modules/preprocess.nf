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