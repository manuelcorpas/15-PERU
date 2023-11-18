-- ABCB1   Chopccas        *2      38

DROP TABLE IF EXISTS `peru`.`06_SUMMARISE_AF_BY_POP` ;

CREATE TABLE IF NOT EXISTS `peru`.`06_SUMMARISE_AF_BY_POP` (
            `ID`         INT NOT NULL AUTO_INCREMENT,
            `GENE`       VARCHAR(20) NOT NULL,
            `POPULATION` VARCHAR(100) DEFAULT NULL,
            `ALLELE`     VARCHAR(50) NOT NULL,
            `COUNTER`    VARCHAR(500) NOT NULL,
            PRIMARY KEY (`ID`)
        )
ENGINE = InnoDB;
