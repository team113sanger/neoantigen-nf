include { REFORMAT_HLA; REFORMAT_VCF } from '../modules/modules.nf'

workflow REFORMAT_DATA {

    take:
        data // VCF file and HLA table.
    
    main:

        // Reformat the HLA table for Neoantimon.
        REFORMAT_HLA(data)

        // Reformat the VEP-annotated VCF for Neoantimon.
        REFORMAT_VCF(REFORMAT_HLA.out.data_for_vcf_reformatting)

    emit:
        data_for_neoantimon = REFORMAT_VCF.out.data_for_neoantimon

}

