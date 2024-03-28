#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process sort_and_tabix {

    publishDir "${params.results}/tabixed"
    container 'docker.io/porchard/general:20220406125608'
    memory '20 GB'
    //clusterOptions '--partition=topmed-working --exclude=topmed,topmed[2-10]'
    executor 'local'

    input:
    path(direct) 
    path(hg38_fasta)
    path(chain)

    output:
    path("direct.txt.gz")
    path("direct.txt.gz.tbi")

    """
    lift-trans-10e4.py $direct $hg38_fasta $chain > direct.unsorted.txt
    head -n 1 direct.unsorted.txt > direct.txt
    cat direct.unsorted.txt | grep -v "^#" | sort -k1,1 -k2n,2 -S 10G >> direct.txt
    bgzip direct.txt
    tabix -p bed direct.txt.gz
    rm direct.unsorted.txt
    """
}




workflow {
    sort_and_tabix(Channel.fromPath(params.direct), Channel.fromPath(params.hg38_fasta), Channel.fromPath(params.chain))
}