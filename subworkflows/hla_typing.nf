include { FILTER_HLA; HLA_TYPING } from '../modules/modules.nf'

workflow GET_HLA_TYPE {

    take:
        data // VCF, BAM, and BAI files.
        bed  // BED file with HLA data.
    
    main:

        // Filter BAM file, keeping only HLA-related sequences.
        // Additionally, convert the BAM file to FASTQ format.
        FILTER_HLA(data, bed)

        // Carry out HLA typing (MHC class I only) with OptiType.
        HLA_TYPING(FILTER_HLA.out.data_for_hla_typing)

    emit:
        data_for_hla_reformatting = HLA_TYPING.out.data_for_hla_reformatting

}

