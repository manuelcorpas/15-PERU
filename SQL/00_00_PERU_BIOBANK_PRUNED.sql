USE `peru2`;

DROP TABLE IF EXISTS `peru2`.`00_00_PERU_BIOBANK_PRUNED` ;

CREATE TABLE 00_00_PERU_BIOBANK_PRUNED (
    codigo VARCHAR(255),
    cod_microarray VARCHAR(255),
    cod_genoma VARCHAR(255),
    107_etnia VARCHAR(255),
    110_cuando_nacio VARCHAR(255),
    112_genero VARCHAR(255),
    115_talla VARCHAR(255),
    116_peso VARCHAR(255),
    123a_fuma_dia VARCHAR(255),
    malaria VARCHAR(255),
    leishmaniosis VARCHAR(255),
    micosis VARCHAR(255),
    tuberculosis VARCHAR(255),
    vih VARCHAR(255),
    artrosis VARCHAR(255),
    diabetes VARCHAR(255),
    hipertension VARCHAR(255),
    acv VARCHAR(255),
    perfil_lipidico VARCHAR(255),
    asma_epoc VARCHAR(255),
    cancer VARCHAR(255),
    enf_mentales VARCHAR(255),
    convulsiones VARCHAR(255),
    a15_observaciones VARCHAR(255),
    b1_glucosa VARCHAR(255),
    b2_colest_total VARCHAR(255),
    b3_triglic VARCHAR(255),
    b4_dhl_col VARCHAR(255),
    b5_ldl_col VARCHAR(255),
    b6_vldl_col VARCHAR(255),
    b7_tgo VARCHAR(255),
    b8_tgp VARCHAR(255),
    b9_bilirr_total VARCHAR(255),
    b10_bilirr_directa VARCHAR(255),
    b11_bilirr_indirecta VARCHAR(255),

PRIMARY KEY (codigo)
)
ENGINE = InnoDB;


