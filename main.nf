nextflow.enable.dsl = 2

process CREATE_HELLO_FILE {

    tag 'hello-genomeops'

    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    output:
    path 'hello.txt', emit: hello_file

    script:
    """
    echo "GenomeOps pipeline works" > hello.txt
    """
}

workflow {
    CREATE_HELLO_FILE()
}
