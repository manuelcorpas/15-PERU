-- drug    fda     ema     swissmedic      hcsc    pmda
USE `peru`;

DROP TABLE IF EXISTS `peru`.`00_PHARMGKB` ;

CREATE TABLE IF NOT EXISTS `peru`.`00_PHARMGKB` (
            `ID_00_PHARMGKB` INT AUTO_INCREMENT,
            `drug`        VARCHAR(500) NOT NULL,
            `fda`         VARCHAR(500) NOT NULL,
            `ema`         VARCHAR(500) NOT NULL,
            `swissmedic`  VARCHAR(500) NOT NULL,
            `hcsc`        VARCHAR(500) NOT NULL,
            `pmda`        VARCHAR(500) NOT NULL,
            PRIMARY KEY (`ID_00_PHARMGKB`)
        )
ENGINE = InnoDB;



