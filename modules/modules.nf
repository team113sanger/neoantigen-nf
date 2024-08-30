process DOWNLOAD_REFFLAT {
    publishDir "${params.outdir}"

    output:
        path("*.grch38.txt"), emit: refflat

    script:
        """
        wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refFlat.txt.gz
        gunzip refFlat.txt.gz
        mv refFlat.txt refFlat.grch38.txt
        """

    stub:
        """
        touch refFlat.grch38.txt
        """
}

process DOWNLOAD_REFMRNA {
    publishDir "${params.outdir}"

    output:
        path("*.grch38.fa"), emit: refmrna

    script:
        """
        wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/refMrna.fa.gz
        gunzip refMrna.fa.gz
        mv refMrna.fa refMrna.grch38.fa
        """

    stub:
        """
        touch refMrna.grch38.fa
        """
}

process FILTER_HLA {
    publishDir "${params.outdir}/fastq", pattern: "*.fq"

    input:
        tuple val(SAMPLE_ID), path(VCF), path(BAM), path(BAI)
        path(BED)

    output:
        tuple val(SAMPLE_ID), path(VCF), path("*_1.fq"), path("*_2.fq"), emit: data_for_hla_typing

    script:
        """
        samtools view -F 4 -b -h ${BAM} --regions-file ${BED} | \
        samtools collate - -u -O | \
        samtools fastq \
            -c 6 \
            -@ 8 \
            -1 HLA_${SAMPLE_ID}_1.fq \
            -2 HLA_${SAMPLE_ID}_2.fq \
            -0 /dev/null \
            -s /dev/null \
            -n
        """

    stub:
        """
        touch HLA_${SAMPLE_ID}_1.fq
        touch HLA_${SAMPLE_ID}_2.fq
        """
}

process HLA_TYPING {
    publishDir "${params.outdir}/hla_typing", pattern: "*/*.tsv"
    errorStrategy "ignore"

    input:
        tuple val(SAMPLE_ID), path(VCF), path(FASTQ1), path(FASTQ2)

    output:
        tuple val(SAMPLE_ID), path(VCF), path("*/*.tsv"), emit: data_for_hla_reformatting
        tuple val(SAMPLE_ID), path(VCF), path("*/*.pdf")

    script:
        """
        OptiTypePipeline.py \
            -i ${FASTQ1} ${FASTQ2} \
            --dna \
            --verbose \
            --outdir \$(pwd)
        """

    stub:
        """
        pwd
        mkdir hla_typing
        touch hla_typing/hla.tsv
        touch hla_typing/hla.pdf
        """
}

process REFORMAT_HLA {
    publishDir "${params.outdir}/hla_typing", pattern: "*_HLA_reformatted.txt"

    input:
        tuple val(SAMPLE_ID), path(VCF), path(HLA)

    output:
        tuple val(SAMPLE_ID), path(VCF), path("*_HLA_reformatted.txt"), emit: data_for_vcf_reformatting

    script:
        """
        #!/usr/bin/env Rscript

        # Read the HLA table.
        readr::read_tsv("${HLA}", col_names = TRUE) |>

            # Create a column called "Name", under which there will be the sample ID.
            dplyr::mutate(Name = "${SAMPLE_ID}") |>

            # Select the columns expected by Neoantimon.
            dplyr::select(Name, A1, A2, B1, B2, C1, C2) |>

            # Save the reformatted HLA table to an output file.
            readr::write_delim(
                file = paste0("${SAMPLE_ID}", "_HLA_reformatted.txt"),
                delim = "\t",
                quote = "none",
                append = FALSE,
                col_names = TRUE
            )
        """

    stub:
        """
        touch ${SAMPLE_ID}_HLA_reformatted.txt
        """
}

process REFORMAT_VCF {
    publishDir "${params.outdir}/vcf", pattern: "*_VCF_reformatted.txt"

    input:
        tuple val(SAMPLE_ID), path(VCF), path(HLA)

    output:
        tuple val(SAMPLE_ID), path("*_VCF_reformatted.txt"), path(HLA), emit: data_for_neoantimon

    script:
        """
        # Create a placeholder for the reformatted VCF data.
        touch ${SAMPLE_ID}_VCF_reformatted.txt

        # Add a header.
        echo -e "Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra" \
            >> ${SAMPLE_ID}_VCF_reformatted.txt

        # Extract information from the VCF.
        bcftools +split-vep ${VCF} \
            -f "${SAMPLE_ID}_VCF_reformatted.txt\t%CHROM:%POS\t%Allele\t%Gene\t%Feature\t%Feature_type\t%Consequence\t%cDNA_position\t%CDS_position\t%Protein_position\t%Amino_acids\t%Codons\t%Existing_variation\tIMPACT=%IMPACT;DISTANCE=%DISTANCE;STRAND=%STRAND\n" \
            --duplicate >> \
            ${SAMPLE_ID}_VCF_reformatted.txt

        # Remove the string "chr" ("Location" column) so that the data looks exactly like what Neoantimon is expecting.
        sed -i "s/${SAMPLE_ID}.smartphase.vep.filt.snpflagged.vcf.gz\tchr/${SAMPLE_ID}.smartphase.vep.filt.snpflagged.vcf.gz\t/g" ${SAMPLE_ID}_VCF_reformatted.txt

        # Remove dinucleotide variation results as these are not compatible with Neoantimon.
        sed -i '/CC/d' ${SAMPLE_ID}_VCF_reformatted.txt
        sed -i '/AA/d' ${SAMPLE_ID}_VCF_reformatted.txt
        sed -i '/TT/d' ${SAMPLE_ID}_VCF_reformatted.txt
        sed -i '/GG/d' ${SAMPLE_ID}_VCF_reformatted.txt
        """

    stub:
        """
        touch ${SAMPLE_ID}_VCF_reformatted.txt
        """
}

process RUN_NEOANTIMON {

    input:
        tuple val(SAMPLE_ID), path(REFORMATTED_VCF), path(REFORMATTED_HLA)
        path(REFFLAT)
        path(REFMRNA)
        path(NET_MHC_PAN)
        path(OUTDIR)

    script:
         """
        #!/usr/bin/env Rscript
        library(Neoantimon)
        library(biomaRt)

        dir.create(${OUTDIR}/neoantimon/${SAMPLE_ID})

        Neoantimon::MainSNVClass1(
            input_vep_format_file = "${REFORMATTED_VCF}",
            hla_file = "${REFORMATTED_HLA}",
            file_name_in_hla_table = "${SAMPLE_ID}",
            export_dir = "${OUTDIR}/neoantimon/${SAMPLE_ID}",
            job_id = "${SAMPLE_ID}",
            refflat_file = "${REFFLAT}",
            refmrna_file = "${REFMRNA}",
            netMHCpan_dir = "${NET_MHC_PAN}"
        )
        """

}