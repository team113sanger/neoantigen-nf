#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { DOWNLOAD_FILES } from "./subworkflows/download_files.nf"
include { GET_HLA_TYPE } from "./subworkflows/hla_typing.nf"
include { REFORMAT_DATA } from "./subworkflows/reformat_data.nf"
include { PREPROCESS_RNA_EXPRESSION; RUN_NEOANTIMON } from "./modules/modules.nf"

workflow {
    
    /*************************/
    /**** Load input data ****/
    /*************************/

    // VCF, BAM, and BAI files.
    data = Channel.fromPath(params.data_files)       
        .splitCsv(header: true, sep: ",")
        .map { row -> tuple(row.sample, row.vcf, row.bam, row.bai) }

    // BED file.
    bed = file(params.bed_file, checkIfExists: true)          

    // RNA expression file.
    rna_expression = file(params.rna_expression, checkIfExists: true)       

    // NetMHCpan program.
    netMHCpan = file(params.net_mhc_pan, checkIfExists: true)   

    /**********************/
    /**** Run analysis ****/
    /**********************/

    // Download files required by Neoantimon.
    DOWNLOAD_FILES()

    // Run HLA typing (MHC class I only).
    GET_HLA_TYPE(data, bed)

    // Pre-process RNA expression table.
    PREPROCESS_RNA_EXPRESSION(rna_expression)

    // Reformat data for Neoantimon.
    REFORMAT_DATA(
        GET_HLA_TYPE.out.data_for_hla_reformatting,
        PREPROCESS_RNA_EXPRESSION.out.preprocessed_rna_expression_table
    )

    // Create output directory for Neoantimon.
    neoantimon_dir = file("${params.outdir}/neoantimon")
    neoantimon_dir.mkdirs()

    // Run Neoantimon analysis.
    RUN_NEOANTIMON(
        REFORMAT_DATA.out.data_for_neoantimon,
        DOWNLOAD_FILES.out.refflat,
        DOWNLOAD_FILES.out.refmrna,
        DOWNLOAD_FILES.out.dna_refseq,
        netMHCpan,
        neoantimon_dir
    )

}
