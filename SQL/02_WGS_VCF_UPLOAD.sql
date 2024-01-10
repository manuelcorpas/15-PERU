DROP TABLE IF EXISTS `peru`.`02_WGS_VCF_UPLOAD` ;

CREATE TABLE IF NOT EXISTS `peru`.`02_WGS_VCF_UPLOAD` (
	`Chromosome`    VARCHAR(40) NOT NULL,
	`Chr_position`  INT NOT NULL,
	`REF`           VARCHAR(255) NOT NULL,
	`ALT`           VARCHAR(255) NOT NULL,
	`ZYG`           TEXT NOT NULL,
	PRIMARY KEY (`Chromosome`,`Chr_position`)
)
ENGINE = InnoDb;
