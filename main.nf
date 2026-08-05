nextflow.enable.dsl = 2


process VALIDATE_SAMPLESHEET {

    tag "${input_samplesheet.simpleName}"

    publishDir "${params.outdir}/pipeline_info", mode: 'copy', overwrite: true

    input:
    path input_samplesheet

    output:
    path 'validated_samplesheet.csv', emit: validated_samplesheet

    script:
    def checkFilesArgument = params.check_files
        ? '--check-files'
        : ''

    """
    python "${projectDir}/bin/validate_samplesheet.py" \
        --input "${input_samplesheet}" \
        --output validated_samplesheet.csv \
        --base-dir "${projectDir}" \
        ${checkFilesArgument}
    """
}


process SUMMARIZE_SAMPLESHEET {

    tag 'validated-input'

    publishDir "${params.outdir}/pipeline_info", mode: 'copy', overwrite: true

    input:
    path validated_samplesheet

    output:
    path 'samplesheet_summary.txt', emit: summary

    script:
    """
    awk 'END { print "Validated rows: " (NR - 1) }' \
        "${validated_samplesheet}" \
        > samplesheet_summary.txt
    """
}


process FASTQC {

    tag "${meta.sample}/${meta.lane}"

    label 'process_low'

    container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'

    publishDir {
        "${params.outdir}/fastqc/${meta.sample}/${meta.lane}"
    }, mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*_fastqc.html'), emit: html

    tuple val(meta), path('*_fastqc.zip'), emit: zip

    path 'fastqc.version.txt', emit: version

    script:
    """
    fastqc \
        --threads ${task.cpus} \
        --outdir . \
        ${reads.join(' ')}

    fastqc --version > fastqc.version.txt
    """
}


process MULTIQC {

    tag 'fastqc-summary'

    label 'process_low'

    container 'multiqc/multiqc:v1.35'

    publishDir "${params.outdir}/multiqc",
        mode: 'copy',
        overwrite: true

    input:
    path 'fastqc/*'

    output:
    path 'multiqc_report.html',
        emit: report

    path 'multiqc_data',
        emit: data

    path 'multiqc.version.txt',
        emit: version

    script:
    """
    multiqc \
        --force \
        --module fastqc \
        --require-logs \
        --outdir . \
        fastqc

    multiqc --version > multiqc.version.txt
    """
}


workflow {

    if (!params.input) {
        error(
            "Missing required parameter: --input\n" + "Example:\n" + "nextflow run main.nf " + "--input assets/samplesheet.example.csv " + "-profile docker"
        )
    }

    input_samplesheet_ch = Channel.fromPath(
        params.input,
        checkIfExists: true
    )

    VALIDATE_SAMPLESHEET(
        input_samplesheet_ch
    )

    validated_samplesheet_ch = VALIDATE_SAMPLESHEET.out.validated_samplesheet

    SUMMARIZE_SAMPLESHEET(
        validated_samplesheet_ch
    )

    fastq_pairs_ch = validated_samplesheet_ch
        .splitCsv(header: true)
        .map { row ->

            def meta = [id: "${row.sample}.${row.lane}", sample: row.sample, lane: row.lane, library: row.library, platform: row.platform, center: row.center]

            def reads = [file(
                row.fastq_1,
                checkIfExists: true
            ), file(
                row.fastq_2,
                checkIfExists: true
            )]

            tuple(meta, reads)
        }

    FASTQC(
        fastq_pairs_ch
    )

    fastqc_archives_ch = FASTQC.out.zip
        .map { meta, zip_files ->
            zip_files
        }
        .flatten()
        .collect()

    MULTIQC(
        fastqc_archives_ch
    )
}
