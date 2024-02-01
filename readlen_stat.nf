/*
 * pipeline input parameters
 */
params.path = "/Users/nzhang/Desktop/downsample/Sample1_1.fastq"
params.outdir = "/Users/nzhang/Desktop/downsample"


process Readlength_Statistics {
    publishDir params.outdir, mode:'copy',overwrite:true

    input:
    path(fastq_address)

    output:
    path '*.stat.txt'

    script:
    """
    if [[ $fastq_address == *gz ]];
    then
        gzcat $fastq_address | awk '{if(NR%4==2) print length(\$1)}' | sort -n | uniq -c > \$(basename $fastq_address).stat.txt
    else
        cat $fastq_address | awk '{if(NR%4==2) print length(\$1)}' | sort -n | uniq -c > \$(basename $fastq_address).stat.txt
    fi
    """
}

workflow {
    Channel
        .fromPath( params.path )
        .set { fastq_ch }

    readlength_stat_ch = Readlength_Statistics(fastq_ch)
}
