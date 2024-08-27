#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { GET_HLA_TYPE } from './subworkflows/hla_typing.nf'
include { REFORMAT_DATA } from './subworkflows/reformat_data.nf'
include { RUN_NEOANTIMON } from './modules/modules.nf'

workflow {
    
    /*************************/
    /**** Load input data ****/
    /*************************/

    bams = Channel.fromPath(params.bam_files, checkIfExists: true) // BAM files.
        .take(1)
        .map { file -> tuple(file.simpleName, file)}
    vcfs = Channel.fromPath(params.vcf_files, checkIfExists: true) // VCF files.
        .take(1)
        .map { file -> tuple(file.simpleName, file)}
    bed = file(params.bed_file, checkIfExists: true)               // BED file.

    /**********************/
    /**** Run analysis ****/
    /**********************/

    // Run HLA typing (MHC class I only).
    GET_HLA_TYPE(bams, bed)

    // Reformat data for Neoantimon.
    REFORMAT_DATA(GET_HLA_TYPE.out.hla_type, vcfs)

    // Run Neoantimon analysis.
    RUN_NEOANTIMON(
        REFORMAT_DATA.out.reformatted_hla_table,
        REFORMAT_DATA.out.reformatted_vcf
    )

}
