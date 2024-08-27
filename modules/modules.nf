process FILTER_HLA {
    publishDir "${params.outdir}/fastq"
    conda "/lustre/scratch124/casm/team113/users/jb62/projects/neoantigen-nf/conda_env"

    input:
        tuple val(SAMPLE_ID), path(BAM)
        path(BED)

    output:
        tuple val(SAMPLE_ID), path("*_1.fq.gz"), emit: fastq_1
        tuple val(SAMPLE_ID), path("*_2.fq.gz"), emit: fastq_2

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
        touch ${SAMPLE_ID}_1.fq
        touch ${SAMPLE_ID}_2.fq
        """
}

process HLA_TYPING {
    publishDir "${params.outdir}/hla_typing"
    conda "/lustre/scratch124/casm/team113/users/jb62/projects/neoantigen-nf/conda_env"

    input:
        tuple val(SAMPLE_ID), path(FASTQ1)
        tuple val(SAMPLE_ID), path(FASTQ2)
        path(OUTDIR)

    output:
        tuple val(SAMPLE_ID), path("${OUTDIR}/*/*.tsv"), emit: hla_table

    script:
        """
        OptiTypePipeline.py \
            -i ${FASTQ1} ${FASTQ2} \
            --dna \
            --verbose \
            --outdir $(pwd)
        """

}

process REFORMAT_HLA {
    publishDir "${params.outdir}/hla_typing"
    conda "/lustre/scratch124/casm/team113/users/jb62/projects/neoantigen-nf/conda_env"

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
    conda "/lustre/scratch124/casm/team113/users/jb62/projects/neoantigen-nf/conda_env"

    input:
        tuple val(SAMPLE_ID), path(VCF)

    output:
        tuple val(SAMPLE_ID), path("*.txt"), emit: reformatted_vcf

    script:
        """
        touch ${SAMPLE_ID}_VCF_reformatted.txt
        echo -e "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra" >> ${SAMPLE_ID}_VCF_reformatted.txt
        bcftools +split-vep ${VCF} \
            -f "$(basename ${VCF})\t%CHROM:%POS\t%Allele\t%Gene\t%Feature\t%Feature_type\t%Consequence\t%cDNA_position\t%CDS_position\t%Protein_position\t%Amino_acids\t%Codons\t%Existing_variation\tIMPACT=%IMPACT;DISTANCE=%DISTANCE;STRAND=%STRAND\n" \
            --duplicate >> \
            ${SAMPLE_ID}_VCF_reformatted.txt
        sed -i "s/HUMAN_GRCh38_full_analysis_set_plus_decoy_hla_pulldown_${SAMPLE_ID}.smartphase.vep.vcf.gz\tchr/HUMAN_GRCh38_full_analysis_set_plus_decoy_hla_pulldown_${SAMPLE_ID}.smartphase.vep.vcf.gz\t/g" ${SAMPLE_ID}_VCF_reformatted.txt
        sed -i '/CC/d' ${SAMPLE_ID}_VCF_reformatted.txt
        sed -i '/AA/d' ${SAMPLE_ID}_VCF_reformatted.txt
        sed -i '/TT/d' ${SAMPLE_ID}_VCF_reformatted.txt
        sed -i '/GG/d' ${SAMPLE_ID}_VCF_reformatted.txt
        """
}

process RUN_NEOANTIMON {
    conda "/lustre/scratch124/casm/team113/users/jb62/projects/neoantigen-nf/conda_env"

    input:
        tuple val(SAMPLE_ID), path(REFORMATTED_HLA)
        tuple val(SAMPLE_ID), path(REFORMATTED_VCF)
        path(OUTDIR)

    script:
        def REFFLAT_FILE = "/lustre/scratch126/casm/team113da/users/jb62/projects/PDX_neoantigen_analysis/data/refFlat.grch38.txt"
        def REFMRNA_FILE = "/lustre/scratch126/casm/team113da/users/jb62/projects/PDX_neoantigen_analysis/data/refMrna.grch38.fa"
        def NET_MHC_PAN = "/lustre/scratch126/casm/team113da/users/jb62/projects/PDX_neoantigen_analysis/software/NetMHCpan/netMHCpan-4.1/netMHCpan"
        """
        #!/usr/bin/env Rscript

        MainSNVClass1(
            input_vep_format_file = "${REFORMATTED_VCF}",
            hla_file = "${REFORMATTED_HLA}",
            file_name_in_hla_table = "${SAMPLE_ID}",
            export_dir = "${OUTDIR}",
            job_id = "${SAMPLE_ID}",
            refflat_file = "${REFFLAT_FILE}",
            refmrna_file = "${REF_MRNA_FILE}",
            netMHCpan_dir = "${NET_MHC_PAN}"
        )
        """

}