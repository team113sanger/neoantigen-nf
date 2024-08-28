#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { DOWNLOAD_FILES } from './subworkflows/download_files.nf'
include { GET_HLA_TYPE } from './subworkflows/hla_typing.nf'
include { REFORMAT_DATA } from './subworkflows/reformat_data.nf'
include { RUN_NEOANTIMON } from './modules/modules.nf'

workflow {
    
    /*************************/
    /**** Load input data ****/
    /*************************/

    bam_and_index = Channel.fromPath(params.bam_and_index)             // BAM+BAI files.
        .splitCsv(header: true, sep: ',')
        .map { row -> tuple(row.sample, row.bam, row.bai) }
        .view()

    vcfs = Channel.fromPath(params.vcf_files, checkIfExists: true)     // VCF files.
        .map { file -> tuple(file.simpleName, file)}

    bed = file(params.bed_file, checkIfExists: true)                   // BED file.

    /**********************/
    /**** Run analysis ****/
    /**********************/

    // Download files required by Neoantimon.
    DOWNLOAD_FILES()

    // Run HLA typing (MHC class I only).
    GET_HLA_TYPE(bam_and_index, bed)

    // Reformat data for Neoantimon.
    REFORMAT_DATA(GET_HLA_TYPE.out.hla_type, vcfs)

    // Run Neoantimon analysis.
    RUN_NEOANTIMON(
        REFORMAT_DATA.out.reformatted_hla_table,
        REFORMAT_DATA.out.reformatted_vcf,
        DOWNLOAD_FILES.out.refflat,
        DOWNLOAD_FILES.out.refmrna,
        params.outdir
    )

}
