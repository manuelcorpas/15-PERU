
USE `peru`;

DROP TABLE IF EXISTS `peru`.`00_05_ALT_ALLELE_POP_COUNT` ;
-- 13:115064342-115064342 {'CHOPCCAS': 0, 'CUSCO': 0, 'IQUITOS': 2, 'MATZES': 0, 'MOCHES': 0, 'TRUJILLO': 0, 'UROS': 0}
CREATE TABLE IF NOT EXISTS `peru`.`00_05_ALT_ALLELE_POP_COUNT` (
	    `LOCATION`      VARCHAR(100) NOT NULL,
            `REF`           VARCHAR(100) NOT NULL,
            `ALT`           VARCHAR(100) NOT NULL,
            `SYMBOL`        VARCHAR(100) NOT NULL,
            `CADD_PHRED`    VARCHAR(100) NOT NULL,
            `Consequence`   VARCHAR(500) NOT NULL,
	    `CLIN_SIG`      VARCHAR(100) NOT NULL,
            `CHOPCCAS`      INT NOT NULL,
	    `CUSCO`         INT NOT NULL,
	    `IQUITOS`       INT NOT NULL,
	    `MATZES`        INT NOT NULL,
            `MOCHES`        INT NOT NULL,
            `TRUJILLO`      INT NOT NULL,
            `UROS`          INT NOT NULL,
            `GNOMAD_AF`     FLOAT NULL,
            `GNOMAD_AF_AMR` FLOAT NULL,
            `GNOMAD_AN`     INT NULL,
            `GNOMAD_AN_AMR` INT NULL,
	    PRIMARY KEY (`LOCATION`)
	)
ENGINE = InnoDB;
