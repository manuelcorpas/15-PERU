DROP TABLE IF EXISTS `peru`.`00_04_VEP_WGS`;

CREATE TABLE IF NOT EXISTS `peru`.`00_04_VEP_WGS` (
    `id`                INT NOT NULL AUTO_INCREMENT,
    `Chromosome`        VARCHAR(40) NOT NULL,
    `Chr_position`      INT NOT NULL,
    `REF`               VARCHAR(255) NOT NULL,      
    `ALT`               VARCHAR(255) NOT NULL,
    `Consequence`       VARCHAR(255) NOT NULL,
    `IMPACT`            VARCHAR(255) NOT NULL,
    `SYMBOL`            VARCHAR(255) NOT NULL,
    `BIOTYPE`           VARCHAR(255) NOT NULL,
    `Existing_variation` VARCHAR(500) NOT NULL,
    `SIFT`              VARCHAR(255) NOT NULL,
    `PolyPhen`          VARCHAR(255) NOT NULL,
    `AF`                VARCHAR(255) NOT NULL,
    `CLIN_SIG`          VARCHAR(255) NOT NULL,
    `LOEUF`             VARCHAR(255) NOT NULL,      -- New column
    `CADD_PHRED`        VARCHAR(255) NOT NULL,      -- New column
    `am_class`          VARCHAR(255) NOT NULL,      -- New column
    `ZYG`               TEXT NOT NULL,
    PRIMARY KEY (`id`),
    UNIQUE INDEX `unique_ref_alt_position` (`Chromosome`, `Chr_position`, `REF`, `ALT`)
)
ENGINE = InnoDB;



