--        Genotype        Phenotype       Haplotype1      Haplotype2      AlternativePhase        VariantData     CNV
--9512704011_R01C01_9512704011_R01C01     *2/*2   Normal Metabolizer      *2;     *2;     ;       *2:default;
--

USE `peru`;

DROP TABLE IF EXISTS `peru`.`04_RUN_PYGPX` ;
-- 13:115064342-115064342 {'CHOPCCAS': 0, 'CUSCO': 0, 'IQUITOS': 2, 'MATZES': 0, 'MOCHES': 0, 'TRUJILLO': 0, 'UROS': 0}
CREATE TABLE IF NOT EXISTS `peru`.`04_RUN_PGX` (
	    `GENE`        VARCHAR(100) NOT NULL,
            `SAMPLE`      VARCHAR(100) NOT NULL,
            `GENOTYPE`    VARCHAR(100) NOT NULL,
            `PHENOTYPE`   VARCHAR(100) NOT NULL,
            `HAPLOTYPE1`  VARCHAR(100) NOT NULL,
            `HAPLOTYPE2`  VARCHAR(100) NOT NULL,
            `ALTERNATIVEPHASE` VARCHAR(100) NOT NULL, 
            `VARIANTDATA` VARCHAR(100) NOT NULL, 
	    `CLIN_SIG`      VARCHAR(100) NOT NULL,
            `CHOPCCAS`      INT NOT NULL,
	    `CUSCO`         INT NOT NULL,
	    `IQUITOS`       INT NOT NULL,
	    `MATZES`        INT NOT NULL,
            `MOCHES`        INT NOT NULL,
            `TRUJILLO`      INT NOT NULL,
            `UROS`          INT NOT NULL,
	    PRIMARY KEY (`LOCATION`)
	)
ENGINE = InnoDB;
