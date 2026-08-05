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