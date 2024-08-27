include { FILTER_HLA; HLA_TYPING } from './modules/modules.nf'

workflow GET_HLA_TYPE {

    take:
        bam // BAM file.
        bed // BED file with HLA data.
    
    main:

        // Filter BAM file, keeping only HLA-related sequences.
        // Additionally, convert the BAM file to FASTQ format.
        FILTER_HLA(bams, bed)

        // Carry out HLA typing (MHC class I only) with OptiType.
        HLA_TYPING(FILTER_HLA.out.fastq_1, FILTER_HLA.out.fastq_2)

    emit:
        // The output of this subworkflow will be the sample's HLA table.
        hla_type = HLA_TYPING.out.hla_table

}

