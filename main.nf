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
        "${params.outdir}/fastqc/" + "${meta.sample}/${meta.lane}"
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

    publishDir "${params.outdir}/multiqc", mode: 'copy', overwrite: true

    input:
    path 'fastqc/*'

    output:
    path 'multiqc_report.html', emit: report

    path 'multiqc_data', emit: data

    path 'multiqc.version.txt', emit: version

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

process PREPARE_REFERENCE_METADATA {

    tag "${reference.simpleName}"

    label 'process_medium'

    container 'quay.io/biocontainers/samtools:1.24--h9dcdb79_1'

    publishDir "${params.outdir}/reference/core",
        mode: 'copy',
        overwrite: true

    input:
    path reference

    output:
    tuple path(reference),
        path("${reference.name}.fai"),
        path("${reference.simpleName}.dict"),
        path('reference.contigs.tsv'),
        path('reference.dictionary.contigs.tsv'),
        path('reference.validation.tsv'),
        path('reference.sha256'),
        emit: bundle

    path 'samtools.reference.version.txt',
        emit: version

    script:
    def dictName = "${reference.simpleName}.dict"

    """
    set -euo pipefail

    rm -f \
        "${reference}.fai" \
        "${dictName}"

    samtools faidx \
        "${reference}"

    samtools dict \
        -a "${params.reference_assembly}" \
        -s "${params.reference_species}" \
        -u "${reference.name}" \
        -o "${dictName}" \
        "${reference}"

    cut \
        -f 1,2 \
        "${reference}.fai" \
        > reference.contigs.tsv

    awk -F '\\t' '
        \$1 == "@SQ" {
            name = ""
            seq_length = ""

            for (i = 2; i <= NF; i++) {
                if (\$i ~ /^SN:/) {
                    name = substr(\$i, 4)
                }

                if (\$i ~ /^LN:/) {
                    seq_length = substr(\$i, 4)
                }
            }

            print name "\\t" seq_length
        }
    ' "${dictName}" \
        > reference.dictionary.contigs.tsv

    diff \
        -u \
        reference.contigs.tsv \
        reference.dictionary.contigs.tsv \
        > reference.dictionary.diff \
    || {
        echo "[ERROR] FASTA index and dictionary differ." >&2
        cat reference.dictionary.diff >&2
        exit 1
    }

    sequence_count=\$(
        wc \
            -l \
            < reference.contigs.tsv
    )

    {
        printf 'check\\tstatus\\n'
        printf 'fai_vs_dictionary\\tPASS\\n'
        printf 'sequence_count\\t%s\\n' \
            "\${sequence_count}"
        printf 'assembly\\t%s\\n' \
            "${params.reference_assembly}"
        printf 'species\\t%s\\n' \
            "${params.reference_species}"
    } > reference.validation.tsv

    sha256sum \
        "${reference}" \
        "${reference}.fai" \
        "${dictName}" \
        reference.contigs.tsv \
        reference.dictionary.contigs.tsv \
        reference.validation.tsv \
        > reference.sha256

    samtools --version \
        > samtools.reference.version.txt
    """
}

process BWA_MEM2_INDEX {

    tag "${reference.simpleName}"

    label 'process_medium'

    container 'quay.io/biocontainers/' + 'bwa-mem2:2.2.1--hd03093a_2'

    publishDir "${params.outdir}/reference/bwa-mem2", mode: 'copy', overwrite: true

    input:
    path reference

    output:
    tuple path(reference), path("${reference.name}.*"), emit: indexed_reference

    path 'bwa-mem2.index.version.txt', emit: version

    script:
    """
    bwa-mem2 index "${reference}"

    bwa-mem2 version \
        > bwa-mem2.index.version.txt \
        2>&1
    """
}


process BWA_MEM2_ALIGN {

    tag "${meta.sample}/${meta.lane}"

    label 'process_medium'

    container 'quay.io/biocontainers/' + 'bwa-mem2:2.2.1--hd03093a_2'

    publishDir {
        "${params.outdir}/alignment/" + "${meta.sample}/${meta.lane}"
    }, mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(reads), path(reference), path(index_files)

    output:
    tuple val(meta), path("${meta.id}.sam"), emit: sam

    tuple val(meta), path("${meta.id}.bwa-mem2.log"), emit: log

    tuple val(meta), path("${meta.id}.bwa-mem2.version.txt"), emit: version

    script:
    def readGroup = ['@RG', "ID:${meta.id}", "SM:${meta.sample}", "LB:${meta.library}", "PL:${meta.platform}", "PU:${meta.lane}", "CN:${meta.center}"].join('\\t')

    """
    bwa-mem2 mem \
        -M \
        -t ${task.cpus} \
        -R '${readGroup}' \
        "${reference}" \
        "${reads[0]}" \
        "${reads[1]}" \
        > "${meta.id}.sam" \
        2> "${meta.id}.bwa-mem2.log"

    bwa-mem2 version \
        > "${meta.id}.bwa-mem2.version.txt" \
        2>&1
    """
}

process SAMTOOLS_SORT_INDEX_QC {

    tag "${meta.sample}/${meta.lane}"

    label 'process_medium'

    container 'quay.io/biocontainers/samtools:1.24--h9dcdb79_1'

    publishDir {
        "${params.outdir}/alignment/" +
        "${meta.sample}/${meta.lane}"
    },
        mode: 'copy',
        overwrite: true

    input:
    tuple val(meta), path(sam)

    output:
    tuple val(meta), path("${meta.id}.sorted.bam"),
        path("${meta.id}.sorted.bam.bai"),
        emit: bam

    tuple val(meta), path("${meta.id}.flagstat.txt"),
        emit: flagstat

    tuple val(meta), path("${meta.id}.idxstats.txt"),
        emit: idxstats

    tuple val(meta), path("${meta.id}.stats.txt"),
        emit: stats

    tuple val(meta), path("${meta.id}.quickcheck.txt"),
        emit: quickcheck

    path 'samtools.version.txt',
        emit: version

    script:
    """
    samtools sort \
        -@ ${task.cpus} \
        -O BAM \
        -o "${meta.id}.sorted.bam" \
        "${sam}"

    samtools index \
        -@ ${task.cpus} \
        "${meta.id}.sorted.bam"

    samtools quickcheck \
        -v \
        "${meta.id}.sorted.bam"

    printf 'OK\\n' \
        > "${meta.id}.quickcheck.txt"

    samtools flagstat \
        -@ ${task.cpus} \
        "${meta.id}.sorted.bam" \
        > "${meta.id}.flagstat.txt"

    samtools idxstats \
        "${meta.id}.sorted.bam" \
        > "${meta.id}.idxstats.txt"

    samtools stats \
        -@ ${task.cpus} \
        "${meta.id}.sorted.bam" \
        > "${meta.id}.stats.txt"

    samtools --version \
        > samtools.version.txt
    """
}

process SAMTOOLS_MARKDUP_LIBRARY {

    tag "${meta.sample}/${meta.library}"

    label 'process_medium'

    container 'quay.io/biocontainers/samtools:1.24--h9dcdb79_1'

    publishDir {
        "${params.outdir}/alignment/" +
        "${meta.sample}/libraries/${meta.library}"
    },
        mode: 'copy',
        overwrite: true,
        pattern: '*.txt'

    input:
    tuple val(meta), path(lane_bams)

    output:
    tuple val(meta),
        path("${meta.sample}.${meta.library}.markdup.bam"),
        emit: bam

    tuple val(meta),
        path("${meta.sample}.${meta.library}.markdup.metrics.txt"),
        emit: metrics

    path 'samtools.markdup.version.txt',
        emit: version

    script:
    def bamList = lane_bams instanceof List
        ? lane_bams
        : [lane_bams]

    def bamInputs = bamList
        .collect { "\"${it}\"" }
        .join(' ')

    def prefix = "${meta.sample}.${meta.library}"

    """
    set -o pipefail

    samtools merge \
        -@ ${task.cpus} \
        -f \
        "${prefix}.merged.bam" \
        ${bamInputs}

    samtools collate \
        -@ ${task.cpus} \
        -O \
        -u \
        "${prefix}.merged.bam" \
    | samtools fixmate \
        -@ ${task.cpus} \
        -m \
        -u \
        - \
        - \
    | samtools sort \
        -@ ${task.cpus} \
        -u \
        - \
    | samtools markdup \
        -@ ${task.cpus} \
        -s \
        -f "${prefix}.markdup.metrics.txt" \
        - \
        "${prefix}.markdup.bam"

    samtools quickcheck \
        -v \
        "${prefix}.markdup.bam"

    samtools --version \
        > samtools.markdup.version.txt
    """
}

process SAMTOOLS_FINALIZE_SAMPLE {

    tag "${sample}"

    label 'process_medium'

    container 'quay.io/biocontainers/samtools:1.24--h9dcdb79_1'

    publishDir {
        "${params.outdir}/alignment/${sample}/final"
    },
        mode: 'copy',
        overwrite: true

    input:
    tuple val(sample),
        path(library_bams),
        path(reference)

    output:
    tuple val(sample),
        path("${sample}.markdup.bam"),
        path("${sample}.markdup.bam.bai"),
        emit: bam

    tuple val(sample),
        path("${sample}.markdup.cram"),
        path("${sample}.markdup.cram.crai"),
        emit: cram

    tuple val(sample),
        path("${sample}.markdup.flagstat.txt"),
        emit: flagstat

    tuple val(sample),
        path("${sample}.markdup.stats.txt"),
        emit: stats

    tuple val(sample),
        path("${sample}.format-check.txt"),
        emit: check

    path 'samtools.final.version.txt',
        emit: version

    script:
    def bamList = library_bams instanceof List
        ? library_bams
        : [library_bams]

    def bamInputs = bamList
        .collect { "\"${it}\"" }
        .join(' ')

    """
    set -o pipefail

    samtools merge \
        -@ ${task.cpus} \
        -f \
        "${sample}.markdup.bam" \
        ${bamInputs}

    samtools index \
        -@ ${task.cpus} \
        "${sample}.markdup.bam"

    samtools view \
        -@ ${task.cpus} \
        -C \
        -T "${reference}" \
        -o "${sample}.markdup.cram" \
        "${sample}.markdup.bam"

    samtools index \
        -@ ${task.cpus} \
        "${sample}.markdup.cram"

    samtools quickcheck \
        -v \
        "${sample}.markdup.bam" \
        "${sample}.markdup.cram"

    bam_records=\$(
        samtools view \
            -c \
            "${sample}.markdup.bam"
    )

    cram_records=\$(
        samtools view \
            -T "${reference}" \
            -c \
            "${sample}.markdup.cram"
    )

    test "\${bam_records}" -eq "\${cram_records}"

    {
        echo "BAM records: \${bam_records}"
        echo "CRAM records: \${cram_records}"
        echo "Status: OK"
    } > "${sample}.format-check.txt"

    samtools flagstat \
        -@ ${task.cpus} \
        "${sample}.markdup.bam" \
        > "${sample}.markdup.flagstat.txt"

    samtools stats \
        -@ ${task.cpus} \
        "${sample}.markdup.bam" \
        > "${sample}.markdup.stats.txt"

    samtools --version \
        > samtools.final.version.txt
    """
}

process CHECK_REFERENCE_COMPATIBILITY {

    tag "${sample}"

    label 'process_low'

    container 'quay.io/biocontainers/samtools:1.24--h9dcdb79_1'

    publishDir {
        "${params.outdir}/reference/validation/${sample}"
    },
        mode: 'copy',
        overwrite: true

    input:
    tuple val(sample),
        path(bam),
        path(bai),
        path(reference),
        path(fai),
        path(dict),
        path(reference_contigs),
        path(dictionary_contigs),
        path(metadata_validation),
        path(checksums)

    output:
    tuple val(sample),
        path("${sample}.reference-compatibility.tsv"),
        emit: report

    tuple val(sample),
        path("${sample}.bam.contigs.tsv"),
        emit: bam_contigs

    tuple val(sample),
        path("${sample}.reference-checksums.txt"),
        emit: checksum_report

    path 'samtools.compatibility.version.txt',
        emit: version

    script:
    """
    set -euo pipefail

    samtools quickcheck \
        -v \
        "${bam}"

    samtools idxstats \
        "${bam}" \
        > /dev/null

    samtools view \
        -H \
        "${bam}" \
    | awk -F '\\t' '
        \$1 == "@SQ" {
            name = ""
            seq_length = ""

            for (i = 2; i <= NF; i++) {
                if (\$i ~ /^SN:/) {
                    name = substr(\$i, 4)
                }

                if (\$i ~ /^LN:/) {
                    seq_length = substr(\$i, 4)
                }
            }

            print name "\\t" seq_length
        }
    ' > "${sample}.bam.contigs.tsv"

    diff \
        -u \
        "${reference_contigs}" \
        "${sample}.bam.contigs.tsv" \
        > "${sample}.reference-contigs.diff" \
    || {
        echo "[ERROR] BAM and reference contigs differ." >&2
        cat "${sample}.reference-contigs.diff" >&2
        exit 1
    }

    sha256sum \
        -c \
        "${checksums}" \
        > "${sample}.reference-checksums.txt"

    reference_sequences=\$(
        wc \
            -l \
            < "${reference_contigs}"
    )

    bam_sequences=\$(
        wc \
            -l \
            < "${sample}.bam.contigs.tsv"
    )

    test \
        "\${reference_sequences}" \
        -eq \
        "\${bam_sequences}"

    {
        printf 'check\\tstatus\\n'
        printf 'fai_vs_dictionary\\tPASS\\n'
        printf 'fai_vs_bam_header\\tPASS\\n'
        printf 'reference_checksums\\tPASS\\n'
        printf 'bam_quickcheck\\tPASS\\n'
        printf 'bam_index\\tPASS\\n'
        printf 'sequence_count\\t%s\\n' \
            "\${reference_sequences}"
    } > "${sample}.reference-compatibility.tsv"

    samtools --version \
        > samtools.compatibility.version.txt
    """
}

workflow {

    if (!params.input) {
        error(
            "Missing required parameter: --input\n" + "Example:\n" + "nextflow run main.nf " + "--input assets/samplesheet.example.csv " + "--reference tests/data/reference/" + "test_reference.fa " + "-profile docker"
        )
    }

    if (!params.reference) {
        error(
            "Missing required parameter: --reference\n" + "Example:\n" + "nextflow run main.nf " + "--input assets/samplesheet.example.csv " + "--reference tests/data/reference/" + "test_reference.fa " + "-profile docker"
        )
    }

    input_samplesheet_ch = Channel.fromPath(
        params.input,
        checkIfExists: true
    )

    reference_ch = Channel.fromPath(
        params.reference,
        checkIfExists: true
    )

    reference_for_cram_ch = Channel.fromPath(
        params.reference,
        checkIfExists: true
    )

    reference_metadata_ch = Channel.fromPath(
        params.reference,
        checkIfExists: true
    )

    VALIDATE_SAMPLESHEET(
        input_samplesheet_ch
    )

    PREPARE_REFERENCE_METADATA(
        reference_metadata_ch
    )

    BWA_MEM2_INDEX(
        reference_ch
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
        .map { meta, zipFiles ->
            zipFiles
        }
        .flatten()
        .collect()

    MULTIQC(
        fastqc_archives_ch
    )

    alignment_inputs_ch = fastq_pairs_ch.combine(
        BWA_MEM2_INDEX.out.indexed_reference
    )

        BWA_MEM2_ALIGN(
        alignment_inputs_ch
    )

        SAMTOOLS_SORT_INDEX_QC(
        BWA_MEM2_ALIGN.out.sam
    )

    library_bams_ch = SAMTOOLS_SORT_INDEX_QC.out.bam
        .map { laneMeta, bam, bai ->

            def libraryMeta = [
                sample: laneMeta.sample,
                library: laneMeta.library
            ]

            tuple(
                libraryMeta,
                bam
            )
        }
        .groupTuple()

    SAMTOOLS_MARKDUP_LIBRARY(
        library_bams_ch
    )

    sample_bams_ch = SAMTOOLS_MARKDUP_LIBRARY.out.bam
        .map { libraryMeta, bam ->
            tuple(
                libraryMeta.sample,
                bam
            )
        }
        .groupTuple()

    sample_finalize_inputs_ch = sample_bams_ch
        .combine(reference_for_cram_ch)

        SAMTOOLS_FINALIZE_SAMPLE(
        sample_finalize_inputs_ch
    )

    reference_compatibility_inputs_ch =
        SAMTOOLS_FINALIZE_SAMPLE.out.bam
            .combine(
                PREPARE_REFERENCE_METADATA.out.bundle
            )

    CHECK_REFERENCE_COMPATIBILITY(
        reference_compatibility_inputs_ch
    )
}
