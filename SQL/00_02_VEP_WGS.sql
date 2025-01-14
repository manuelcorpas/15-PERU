USE `peru`;

DROP TABLE IF EXISTS `peru`.`00_02_VEP_WGS`;

CREATE TABLE 00_02_VEP_WGS (
    ID INT AUTO_INCREMENT PRIMARY KEY,
    Uploaded_variation VARCHAR(255),
    Location VARCHAR(255),
    Allele CHAR(1),
    Consequence VARCHAR(255),
    IMPACT VARCHAR(255),
    SYMBOL VARCHAR(255),
    Gene VARCHAR(255),
    Feature_type VARCHAR(255),
    Feature VARCHAR(255),
    BIOTYPE VARCHAR(255),
    EXON VARCHAR(255),
    INTRON VARCHAR(255),
    HGVSc VARCHAR(255),
    HGVSp VARCHAR(255),
    cDNA_position VARCHAR(255),
    CDS_position VARCHAR(255),
    Protein_position VARCHAR(255),
    Amino_acids VARCHAR(255),
    Codons VARCHAR(255),
    Existing_variation VARCHAR(2000),
    REF_ALLELE VARCHAR(255),            -- New column
    UPLOADED_ALLELE VARCHAR(255),       -- New column
    DISTANCE VARCHAR(255),
    STRAND VARCHAR(255),
    FLAGS VARCHAR(255),
    SYMBOL_SOURCE VARCHAR(255),
    HGNC_ID VARCHAR(255),
    SIFT VARCHAR(255),
    PolyPhen VARCHAR(255),
    AF VARCHAR(255),
    gnomADg_AF VARCHAR(255),
    gnomADg_AFR_AF VARCHAR(255),
    gnomADg_AMI_AF VARCHAR(255),
    gnomADg_AMR_AF VARCHAR(255),
    gnomADg_ASJ_AF VARCHAR(255),
    gnomADg_EAS_AF VARCHAR(255),
    gnomADg_FIN_AF VARCHAR(255),
    gnomADg_MID_AF VARCHAR(255),
    gnomADg_NFE_AF VARCHAR(255),
    gnomADg_OTH_AF VARCHAR(255),
    gnomADg_SAS_AF VARCHAR(255),
    CLIN_SIG VARCHAR(255),
    SOMATIC VARCHAR(255),
    PHENO VARCHAR(255),
    MOTIF_NAME VARCHAR(255),
    MOTIF_POS VARCHAR(255),
    HIGH_INF_POS VARCHAR(255),
    MOTIF_SCORE_CHANGE VARCHAR(255),
    TRANSCRIPTION_FACTORS VARCHAR(255),
    LOEUF VARCHAR(255),                 -- New column
    CADD_PHRED VARCHAR(255),            -- New column
    CADD_RAW VARCHAR(255),              -- New column
    am_class VARCHAR(255),              -- New column
    am_pathogenicity VARCHAR(255),      -- New column
    INDEX(Location),
    INDEX(IMPACT),
    INDEX(SYMBOL)
)
ENGINE = InnoDB;

