# GenomeOps synthetic test data

The files in `fastq/` contain small synthetic paired-end sequencing
reads created exclusively for automated pipeline tests.

They do not originate from a real sequencing experiment and must not
be used for biological or clinical interpretation.

Test layout:

- sample: `HG002`
- lanes: `L001` and `L002`
- sequencing type: paired-end
- reads per FASTQ file: 2
- read length: 40 bases

## Synthetic variant dataset

The synthetic test dataset is generated with the following command:

    python scripts/generate_synthetic_variant_data.py

The dataset contains:

- one synthetic contig named `chrTest`;
- a reference length of 10,000 bp;
- one sample named `HG002`;
- two sequencing lanes: `L001` and `L002`;
- 30 paired-end read pairs per lane;
- reads with a length of 150 bp;
- one synthetic heterozygous SNP at `chrTest:5000`;
- truth files stored in `tests/data/variants/`.

All files in this dataset are artificial and must not be used for biological or clinical interpretation.