#!/usr/bin/env python3

from __future__ import annotations

from pathlib import Path
from random import Random


PROJECT_ROOT = Path(__file__).resolve().parents[1]

FASTQ_DIR = PROJECT_ROOT / "tests/data/fastq"
REFERENCE_DIR = PROJECT_ROOT / "tests/data/reference"
VARIANT_DIR = PROJECT_ROOT / "tests/data/variants"

REFERENCE_PATH = REFERENCE_DIR / "test_reference.fa"
TRUTH_VCF_PATH = VARIANT_DIR / "synthetic_truth.vcf"
TRUTH_TSV_PATH = VARIANT_DIR / "synthetic_truth.tsv"

SEED = 20260806
REFERENCE_LENGTH = 10_000
VARIANT_POSITION = 5_000  # One-based VCF position.
READ_LENGTH = 150
FRAGMENT_LENGTH = 350
PAIRS_PER_LANE = 30

SAMPLE = "HG002"
LANES = ("L001", "L002")
DNA = "ACGT"


def reverse_complement(sequence: str) -> str:
    table = str.maketrans(
        "ACGT",
        "TGCA",
    )
    return sequence.translate(table)[::-1]


def write_fasta(
    path: Path,
    sequence: str,
) -> None:
    wrapped = "\n".join(
        sequence[index : index + 80]
        for index in range(
            0,
            len(sequence),
            80,
        )
    )

    path.write_text(
        f">chrTest\n{wrapped}\n",
        encoding="utf-8",
    )


def write_fastq_record(
    handle: object,
    *,
    name: str,
    sequence: str,
    mate: int,
) -> None:
    quality = "I" * len(sequence)

    handle.write(
        f"@{name}/{mate}\n"
        f"{sequence}\n"
        "+\n"
        f"{quality}\n"
    )


def main() -> None:
    rng = Random(SEED)

    FASTQ_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )
    REFERENCE_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )
    VARIANT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    reference = "".join(
        rng.choice(DNA)
        for _ in range(REFERENCE_LENGTH)
    )

    variant_index = VARIANT_POSITION - 1
    reference_allele = reference[variant_index]

    alternate_allele = {
        "A": "C",
        "C": "G",
        "G": "T",
        "T": "A",
    }[reference_allele]

    write_fasta(
        REFERENCE_PATH,
        reference,
    )

    total_pairs = 0
    alternate_pairs = 0

    for lane_index, lane in enumerate(LANES):
        r1_path = (
            FASTQ_DIR
            / f"{SAMPLE}_{lane}_R1.fastq"
        )
        r2_path = (
            FASTQ_DIR
            / f"{SAMPLE}_{lane}_R2.fastq"
        )

        with (
            r1_path.open(
                "w",
                encoding="utf-8",
            ) as r1_handle,
            r2_path.open(
                "w",
                encoding="utf-8",
            ) as r2_handle,
        ):
            for pair_index in range(
                PAIRS_PER_LANE
            ):
                global_pair_index = (
                    lane_index * PAIRS_PER_LANE
                    + pair_index
                )

                fragment_start = (
                    4_850
                    + global_pair_index
                )

                fragment_end = (
                    fragment_start
                    + FRAGMENT_LENGTH
                )

                r1_sequence = reference[
                    fragment_start :
                    fragment_start + READ_LENGTH
                ]

                r2_template = reference[
                    fragment_end - READ_LENGTH :
                    fragment_end
                ]

                r2_sequence = reverse_complement(
                    r2_template
                )

                carries_alternate = (
                    global_pair_index % 2 == 0
                )

                local_variant_index = (
                    variant_index
                    - fragment_start
                )

                if not (
                    0
                    <= local_variant_index
                    < READ_LENGTH
                ):
                    raise RuntimeError(
                        "Synthetic R1 does not cover "
                        "the configured variant."
                    )

                if carries_alternate:
                    r1_bases = list(r1_sequence)

                    r1_bases[
                        local_variant_index
                    ] = alternate_allele

                    r1_sequence = "".join(
                        r1_bases
                    )

                    alternate_pairs += 1

                read_name = (
                    f"{SAMPLE}_{lane}_"
                    f"{pair_index + 1:04d}"
                )

                write_fastq_record(
                    r1_handle,
                    name=read_name,
                    sequence=r1_sequence,
                    mate=1,
                )

                write_fastq_record(
                    r2_handle,
                    name=read_name,
                    sequence=r2_sequence,
                    mate=2,
                )

                total_pairs += 1

    TRUTH_VCF_PATH.write_text(
        "##fileformat=VCFv4.2\n"
        "##source=GenomeOpsSyntheticData\n"
        f"##contig=<ID=chrTest,length="
        f"{REFERENCE_LENGTH}>\n"
        "##INFO=<ID=SYNTHETIC,Number=0,"
        "Type=Flag,Description=\"Synthetic truth "
        "variant\">\n"
        "##FORMAT=<ID=GT,Number=1,Type=String,"
        "Description=\"Genotype\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\t"
        "FILTER\tINFO\tFORMAT\tHG002\n"
        f"chrTest\t{VARIANT_POSITION}\t"
        f"synthetic_snp\t{reference_allele}\t"
        f"{alternate_allele}\t60\tPASS\t"
        "SYNTHETIC\tGT\t0/1\n",
        encoding="utf-8",
    )

    TRUTH_TSV_PATH.write_text(
        "chrom\tposition\tref\talt\tgenotype\n"
        f"chrTest\t{VARIANT_POSITION}\t"
        f"{reference_allele}\t"
        f"{alternate_allele}\t0/1\n",
        encoding="utf-8",
    )

    print(
        f"[OK] Reference: {REFERENCE_PATH}"
    )
    print(
        f"[OK] Reference length: "
        f"{REFERENCE_LENGTH}"
    )
    print(
        f"[OK] Truth SNP: "
        f"chrTest:{VARIANT_POSITION} "
        f"{reference_allele}>"
        f"{alternate_allele}"
    )
    print(
        f"[OK] Read pairs: {total_pairs}"
    )
    print(
        f"[OK] Alternate-supporting pairs: "
        f"{alternate_pairs}"
    )


if __name__ == "__main__":
    main()