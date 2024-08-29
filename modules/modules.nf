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
    publishDir "${params.outdir}/fastq"

    input:
        tuple val(SAMPLE_ID), path(BAM), path(BAI)
        path(BED)

    output:
        tuple val(SAMPLE_ID), path("*_1.fq"), emit: fastq_1
        tuple val(SAMPLE_ID), path("*_2.fq"), emit: fastq_2

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
    publishDir "${params.outdir}/hla_typing"

    input:
        tuple val(SAMPLE_ID), path(FASTQ1)
        tuple val(SAMPLE_ID), path(FASTQ2)

    output:
        tuple val(SAMPLE_ID), path("\$(pwd)/*/*.tsv"), emit: hla_table
        tuple val(SAMPLE_ID), path("\$(pwd)/*/*.pdf"), emit: figure

    script:
        """
        OptiTypePipeline.py \
            -i ${FASTQ1} ${FASTQ2} \
            --dna \
            --verbose \
            --outdir \$(pwd)
        """

}

process REFORMAT_HLA {
    publishDir "${params.outdir}/hla_typing"

    input:
        tuple val(SAMPLE_ID), path(HLA)

    output:
        tuple val(SAMPLE_ID), path("*.txt"), emit: reformatted_hla_table

    script:
        """
        #!/usr/bin/env Rscript

        readr::read_tsv("${HLA}", col_names = TRUE) |>
        dplyr::mutate(Name = "${SAMPLE_ID}") |>
        dplyr::select(Name, A1, A2, B1, B2, C1, C2) |>
        readr::write_delim(
            file = paste0("${SAMPLE_ID}", "_HLA_reformatted.txt"),
            delim = "\t",
            quote = "none",
            append = FALSE,
            col_names = TRUE
        )
        """

}

process REFORMAT_VCF {
    publishDir "${params.outdir}/vcf"

    input:
        tuple val(SAMPLE_ID), path(VCF)

    output:
        tuple val(SAMPLE_ID), path("*.txt"), emit: reformatted_vcf

    script:
        """
        touch ${SAMPLE_ID}_VCF_reformatted.txt
        echo -e "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra" \
            >> ${SAMPLE_ID}_VCF_reformatted.txt
        bcftools +split-vep ${VCF} \
            -f "\$(basename ${VCF})\t%CHROM:%POS\t%Allele\t%Gene\t%Feature\t%Feature_type\t%Consequence\t%cDNA_position\t%CDS_position\t%Protein_position\t%Amino_acids\t%Codons\t%Existing_variation\tIMPACT=%IMPACT;DISTANCE=%DISTANCE;STRAND=%STRAND\n" \
            --duplicate >> \
            ${SAMPLE_ID}_VCF_reformatted.txt
        sed -i "s/${SAMPLE_ID}.smartphase.vep.vcf.gz\tchr/${SAMPLE_ID}.smartphase.vep.vcf.gz\t/g" ${SAMPLE_ID}_VCF_reformatted.txt
        sed -i '/CC/d' ${SAMPLE_ID}_VCF_reformatted.txt
        sed -i '/AA/d' ${SAMPLE_ID}_VCF_reformatted.txt
        sed -i '/TT/d' ${SAMPLE_ID}_VCF_reformatted.txt
        sed -i '/GG/d' ${SAMPLE_ID}_VCF_reformatted.txt
        """
}

process RUN_NEOANTIMON {

    input:
        tuple val(SAMPLE_ID), path(REFORMATTED_HLA)
        tuple val(SAMPLE_ID), path(REFORMATTED_VCF)
        path(REFFLAT)
        path(REFMRNA)
        path(OUTDIR)

    script:
        def NET_MHC_PAN = "/lustre/scratch126/casm/team113da/users/jb62/projects/PDX_neoantigen_analysis/software/NetMHCpan/netMHCpan-4.1/netMHCpan"
        """
        #!/usr/bin/env Rscript

        MainSNVClass1(
            input_vep_format_file = "${REFORMATTED_VCF}",
            hla_file = "${REFORMATTED_HLA}",
            file_name_in_hla_table = "${SAMPLE_ID}",
            export_dir = "${OUTDIR}",
            job_id = "${SAMPLE_ID}",
            refflat_file = "${REFFLAT}",
            refmrna_file = "${REFMRNA}",
            netMHCpan_dir = "${NET_MHC_PAN}"
        )
        """

}