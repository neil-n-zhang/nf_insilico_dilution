/*
 * Downsample sample A and sample B and combine the reads
 */
params.fqA_R1 = "./Sample1_1.fastq"
params.fqA_R2 = "./Sample1_2.fastq"

params.fqB_R1 = "./Sample2_50bp_R1.fastq"
params.fqB_R2 = "./Sample2_50bp_R2.fastq"

// The number of bases in the final fastq files
params.totalbases = 20000
// The fraction of bases in the final fastq that come from sample A
params.fqA_fraction = 0.1

params.outdir = "./output"

include { Downsample_R1R2 as Downsample_fqA_R1R2 } from './fastq_downsample.nf'
include { Downsample_R1R2 as Downsample_fqB_R1R2 } from './fastq_downsample.nf'

process Merge_fastqs_R1 {
    publishDir params.outdir, mode:'copy',overwrite:true

    input:
    path(fqA_R1)
    path(fqB_R1)

    output:
    path 'diluted_R1.fq'

    script:
    """
    cat $fqA_R1 $fqB_R1 > diluted_R1.fq
    """
}

process Merge_fastqs_R2 {
    publishDir params.outdir, mode:'copy',overwrite:true

    input:
    path(fqA_R2)
    path(fqB_R2)

    output:
    path 'diluted_R2.fq'

    script:
    """
    cat $fqA_R2 $fqB_R2 > diluted_R2.fq
    """
}

workflow {

    // Need to use () since we have two channels here, otherwise, we will get error information: Unexpected input: '{'
    (ds_fqA_R1, ds_fqA_R2)=Downsample_fqA_R1R2(params.fqA_R1, params.fqA_R2, params.totalbases*params.fqA_fraction)
    (ds_fqB_R1, ds_fqB_R2)=Downsample_fqB_R1R2(params.fqB_R1, params.fqB_R2, params.totalbases*(1-params.fqA_fraction))
    Merge_fastqs_R1(ds_fqA_R1,ds_fqB_R1)
    Merge_fastqs_R2(ds_fqA_R2,ds_fqB_R2)


}