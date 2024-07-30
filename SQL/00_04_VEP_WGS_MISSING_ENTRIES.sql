DROP TABLE IF EXISTS `peru`.`00_04_VEP_WGS_MISSING_ENTRIES`;

-- Create the new table to store missing entries
CREATE TABLE IF NOT EXISTS `peru`.`00_04_VEP_WGS_MISSING_ENTRIES` (
    `Chromosome`    VARCHAR(40) NOT NULL,
    `Chr_position`  INT NOT NULL,
    `REF`           VARCHAR(255) NOT NULL,
    `ALT`           VARCHAR(255) NOT NULL,
    `ZYG`           TEXT NOT NULL,
    PRIMARY KEY (`Chromosome`, `Chr_position`, `REF`, `ALT`)
)
ENGINE = InnoDB;
