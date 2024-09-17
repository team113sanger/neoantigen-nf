include { REFORMAT_HLA; REFORMAT_VCF; EXTRACT_RNA_EXPRESSION } from '../modules/modules.nf'

workflow REFORMAT_DATA {

    take:
        data // VCF file and HLA table.
        rna_expression // Count matrix.

    main:

        // Reformat the HLA table for Neoantimon.
        REFORMAT_HLA(data)

        // Reformat the VEP-annotated VCF for Neoantimon.
        REFORMAT_VCF(REFORMAT_HLA.out.data_for_vcf_reformatting)

        // Extract
        EXTRACT_RNA_EXPRESSION(
            REFORMAT_VCF.out.data_for_gene_expression_extraction,
            rna_expression
        )

    emit:
        data_for_neoantimon = EXTRACT_RNA_EXPRESSION.out.data_for_neoantimon

}

