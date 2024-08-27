include { REFORMAT_HLA; REFORMAT_VCF } from './modules/modules.nf'

workflow REFORMAT_DATA {

    take:
        hla_table      // HLA table.
        annotated_vcf  // VEP-annotated VCF.
    
    main:

        // Reformat the HLA table for Neoantimon.
        REFORMAT_HLA(hla_table)

        // Reformat the VEP-annotated VCF for Neoantimon.
        REFORMAT_VCF(annotated_vcf)

    emit:
        // The output of this subworkflow will be the sample's reformatted HLA table and
        // VEP-annotated VCF.
        reformatted_hla_table = REFORMAT_HLA.out.reformatted_hla_table
        reformatted_vcf = REFORMAT_HLA.out.reformatted_vcf

}

