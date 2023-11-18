import csv
import mariadb
import os
import configparser
import pypgx

config   = configparser.ConfigParser()
config.read_file(open(r'CONF/mariadb.conf'))

try:
    db = mariadb.connect(
        host   = config.get('peru','host'),
        user   = config.get('peru','user'     ),
        passwd = config.get('peru','password' ),
        db     = config.get('peru','database'))

except mariadb.Error as e:
    print(f"Error connecting to MariaDB Platform: {e}")
    sys.exit(1)

cursor = db.cursor(dictionary=True)

# SELECT GENE, EFFECT, COUNT(EFFECT) AS EffectCount FROM 05_PGX_INSERT GROUP BY GENE, EFFECT;

sql1 = "SELECT SAMPLE,GENE,ALLELE,EFFECT FROM 05_PGX_INSERT"
sql2 = "SELECT POPULATION FROM 00_SAMPLE WHERE ARRAY_CODE = '{0}'"

cursor.execute(sql1)
pgx_effects = cursor.fetchall()

gene ={}

for effect in pgx_effects:
    Smpl    = effect['SAMPLE'].split('_')
    SAMPLE  = Smpl[0] + '_' + Smpl[1]
    GENE    = effect['GENE']
    EFFECT  = effect['EFFECT']
    ALLELES = effect['ALLELE'].split('/')
    
    cursor.execute(sql2.format(SAMPLE))
    POPULATION = ''
    populations = cursor.fetchall()
    for pop in populations:
        POPULATION = pop['POPULATION']
    #print(GENE,SAMPLE,POPULATION,*ALLELES,EFFECT,sep='\t')
    if not GENE in gene:
        gene[GENE] = {}
    if not POPULATION in gene[GENE]:
        gene[GENE][POPULATION] = {}
    if not ALLELES[0] in gene[GENE][POPULATION]:
        gene[GENE][POPULATION][ALLELES[0]] = 1
    else:
        gene[GENE][POPULATION][ALLELES[0]] += 1
    if not ALLELES[1] in gene[GENE][POPULATION]:
        gene[GENE][POPULATION][ALLELES[1]] = 1
    else:
        gene[GENE][POPULATION][ALLELES[1]] += 1 

for g in gene:
    sorted_pop = sorted(gene[g].keys())
    for p in sorted_pop:
        for a in gene[g][p]:
            print(g,p,a,gene[g][p][a],sep='\t')





