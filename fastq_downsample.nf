/*
 * Downsample the fastq to contain a certain number of bases
 * Need to perform adapter trimming after downsampling to make sure all reads still have the same length
 */
params.R1 = "./Sample1_1.fastq"
params.R2 = "./Sample1_2.fastq"
params.dsbases = 10000

params.outdir = "./output"

seqtk_path="/Users/nzhang/bioinfotools/seqtk/seqtk"

include { Readlength_Statistics as Readlength_Statistics_R1 } from './readlen_stat.nf'
include { Readlength_Statistics as Readlength_Statistics_R2 } from './readlen_stat.nf'

process Readlength_Check_R1R2 {

    input:
    path(readlength_stat_R1)
    path(readlength_stat_R2)

    output:
    true

    script:
    """
    if cmp -s $readlength_stat_R1 $readlength_stat_R2; then
        echo "R1 and R2 files match."
    else
        echo "Error: Read1 and read2 files have different read length."
        exit 1
    fi

    if [ \$(wc -l $readlength_stat_R1 | awk '{print \$1}') != 1 ]; then 
        echo "Error: Reads don't have the same length."
        exit 1
    fi

    if [ $params.dsbases -gt \$(cat $readlength_stat_R1 | awk '{print \$1*\$2}') ]; then 
        echo "Error: $params.R1 doesn't have enough reads."
        exit 1
    fi
    """
}

process Seqtk_downsample_R1 {
    publishDir params.outdir, mode:'copy',overwrite:true

    input:
    path(fastq_address)
    path(readlength_stat)
    val(readlength_check)
    val(dsbases)

    output:
    path '*_ds_*fq'

    // Shell block: it uses the exclamation mark ! character, instead of the usual dollar $ character, to denote Nextflow variables.
    shell:
    '''
    dsbases_int=$(echo "!{dsbases}/1" | bc) 
    reads_length=$(cat !{readlength_stat} | awk '{print $2}')
    output_fq_name=$(basename !{fastq_address}_ds_$dsbases_int.fq)
    # Celling of read number
    reads_number=$(awk "BEGIN { printf \\"%.2f\\", ($dsbases_int / 2 / $reads_length) }" | awk '{printf("%d",$0+=$0<0?0:0.999)}')
    "!{seqtk_path}" sample -s100 !{fastq_address} $reads_number > $output_fq_name

    '''
}

process Seqtk_downsample_R2 {
    publishDir params.outdir, mode:'copy',overwrite:true

    input:
    path(fastq_address)
    path(readlength_stat)
    val(readlength_check)
    val(dsbases)

    output:
    path '*_ds_*fq'

    shell:
    '''
    dsbases_int=$(echo "!{dsbases}/1" | bc) 
    reads_length=$(cat !{readlength_stat} | awk '{print $2}')
    output_fq_name=$(basename !{fastq_address}_ds_$dsbases_int.fq)
    # Celling of read number
    reads_number=$(awk "BEGIN { printf \\"%.2f\\", ($dsbases_int / 2 / $reads_length) }" | awk '{printf("%d",$0+=$0<0?0:0.999)}')
    "!{seqtk_path}" sample -s100 !{fastq_address} $reads_number > $output_fq_name

    '''
}

workflow Downsample_R1R2{
    take:
    R1
    R2
    dsbases

    main:
    readlength_stat_ch_R1 = Readlength_Statistics_R1(Channel.fromPath(R1))
    readlength_stat_ch_R2 = Readlength_Statistics_R2(Channel.fromPath(R2))
    Readlength_Check_R1R2(readlength_stat_ch_R1,readlength_stat_ch_R2)
    Seqtk_downsample_R1(Channel.fromPath(R1),readlength_stat_ch_R1,Readlength_Check_R1R2,dsbases)
    Seqtk_downsample_R2(Channel.fromPath(R2),readlength_stat_ch_R1,Readlength_Check_R1R2,dsbases)

    emit:
    Seqtk_downsample_R1.out
    Seqtk_downsample_R2.out
}

workflow {

    Downsample_R1R2(params.R1,params.R2,params.dsbases)

}