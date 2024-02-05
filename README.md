# nf_insilico_dilution

A nextflow pipeline to perform in silico dilution between two fastq samples.
Example:
nextflow run insilico_dilution.nf --fqA_R1 Sample1_1.fastq --fqA_R2 Sample1_2.fastq --fqB_R1 Sample2_50bp_R1.fastq --fqB_R2 Sample2_50bp_R2.fastq --totalbases 10000 --fqA_fraction 0.15