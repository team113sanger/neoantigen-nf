include { DOWNLOAD_REFFLAT; DOWNLOAD_REFMRNA } from '../modules/modules.nf'

workflow DOWNLOAD_FILES {

    main:

        DOWNLOAD_REFFLAT() //Download refFlat file for Neoantimon.
        DOWNLOAD_REFMRNA() //Download refMrna file for Neoantimon.

    emit:
        refflat = DOWNLOAD_REFFLAT.out.refflat
        refmrna = DOWNLOAD_REFMRNA.out.refmrna
}

