-- POPULATION,ALLELE,ALLELES_OBSERVED,ALLELES_TOTAL,FREQUENCY

DROP TABLE IF EXISTS `peru`.`07_UKBB_POP_AF_INSERT` ;

CREATE TABLE IF NOT EXISTS `peru`.`07_UKBB_POP_AF_INSERT` (
            `ID`               INT NOT NULL AUTO_INCREMENT,
            `GENE`             VARCHAR(30) NOT NULL,
            `POPULATION`       VARCHAR(200) DEFAULT NULL,
            `ALLELE`           VARCHAR(200) NOT NULL,
            `ALLELES_OBSERVED` INT NOT NULL,
            `ALLELES_TOTAL`    INT NOT NULL,
            `FREQUENCY`        FLOAT(20) NOT NULL,
            PRIMARY KEY (`ID`)
        )
ENGINE = InnoDB;
