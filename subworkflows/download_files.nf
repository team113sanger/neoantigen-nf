include { DOWNLOAD_REFFLAT; DOWNLOAD_REFMRNA; DOWNLOAD_DNA_REFSEQ} from '../modules/modules.nf'

workflow DOWNLOAD_FILES {

    main:

        DOWNLOAD_REFFLAT() // Download refFlat file for Neoantimon.
        DOWNLOAD_REFMRNA() // Download refMrna file for Neoantimon.
        DOWNLOAD_DNA_REFSEQ() // Download DNA RefSeq file for Neoantimon.

    emit:
        refflat = DOWNLOAD_REFFLAT.out.refflat
        refmrna = DOWNLOAD_REFMRNA.out.refmrna
        dna_refseq = DOWNLOAD_DNA_REFSEQ.out.dna_refseq
}

