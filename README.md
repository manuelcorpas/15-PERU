# Peruvian Genome Project: Scripts Documentation

**Manuscript**: *Genomic Insights into 30 Peruvian Indigenous and Mestizo Populations*

This repository contains Python scripts for analysing whole-genome sequencing (WGS) and SNP array data from 30 Peruvian populations (n = 736 unrelated individuals after quality control). The analyses support population genetics, pharmacogenomics, and variant discovery findings reported in the manuscript.

---

## Table of Contents

1. [Project Overview](#project-overview)
2. [Directory Structure](#directory-structure)
3. [Script Documentation](#script-documentation)
   - [00-01: Genetic Diversity and SNP Analysis](#00-01-genetic-diversity-and-snp-analysis)
   - [00-02: VCF Annotation](#00-02-vcf-annotation)
   - [00-03: Population Structure](#00-03-population-structure)
   - [00-04: Pharmacogenomics](#00-04-pharmacogenomics)
   - [00-05: gnomAD Validation](#00-05-gnomad-validation)
   - [00-06: Variant Annotation and Statistics](#00-06-variant-annotation-and-statistics)
   - [00-07: Figure Generation](#00-07-figure-generation)
4. [Dependencies](#dependencies)
5. [Reproducibility](#reproducibility)

---

## Project Overview

This project analyses genomic data from 30 Peruvian populations representing Peru's three ecological regions:
- **Coastal**: Tumbes, Trujillo, Lima, Arequipa, Moquegua, Tacna, Lambayeque
- **Highland (Andean)**: Cusco, Ayacucho, Huaraz, Puno, Chopccas, Chachapoyas, Jaqarus, Queros, Uros
- **Amazonian**: Matzes, Machiguenga, Ashaninka, Awajun, Shipibo, Candoshi, Nahua, Lamas, Iquitos

Plus one Afro-Peruvian population (AFRODESCENDENTS) sampled from Chincha province.

**Sample sizes**:
- Initial recruitment: 1,149 individuals
- After quality control and kinship filtering: 736 unrelated individuals
- WGS samples: 150 (109 after kinship filtering)
- Array-genotyped samples: 873 (627 after kinship filtering)

---

## Directory Structure

```
HEINNER-GUIO/
├── PYTHON/
│   ├── ARCHIVE/                    # Deprecated scripts
│   ├── 00-01-*.py                  # Genetic diversity analyses
│   ├── 00-02-*.py                  # VCF annotation
│   ├── 00-03-*.py                  # Population structure (IBD, ADMIXTURE, PCA, Fst, f3)
│   ├── 00-04-*.py                  # Pharmacogenomics (PyPGx)
│   ├── 00-05-*.py                  # gnomAD validation
│   ├── 00-06-*.py                  # Variant annotation and statistics
│   └── 00-07-*.py                  # Figure generation
│
├── ANALYSIS/
│   ├── 00-00-00-DEMOGRAPHICS/      # Sample metadata
│   ├── 00-00-WGS-STATS/            # WGS quality metrics
│   ├── 00-01-GEN-DIV/              # Genetic diversity outputs
│   ├── 00-05-PyPGx/                # Pharmacogenomics results
│   ├── 00-06-IBD/                  # Identity-by-descent analysis
│   ├── 00-07-ADMIXTURE/            # Ancestry estimation
│   ├── 00-08-PCA/                  # Principal component analysis
│   ├── 00-09-FST/                  # Fst differentiation
│   ├── 00-10-GNOMAD-CHECK/         # gnomAD validation
│   ├── 00-11-GNOMAD-NOVEL/         # Novel variant verification
│   ├── 00-12-PGX-STATS/            # Pharmacogenomics statistics
│   ├── 00-13-DSTATS/               # Outgroup f3 statistics
│   └── 00-17-FIGURE1/              # Main figure outputs
│
└── DATA/                           # Input data (not included)
```

---

## Script Documentation

### 00-01: Genetic Diversity and SNP Analysis

Scripts for computing variant counts, heterozygosity, and population-specific metrics.

| Script | Purpose | Output | Figure/Table |
|--------|---------|--------|--------------|
| `00-01-00-gen-div.py` | Computes genetic diversity metrics (heterozygosity, allele richness) per population | Summary statistics | Table 1, Figure 2 |
| `00-01-01-snp-counts.py` | Tallies SNPs per sample and population | Count matrices | Figure 2A |
| `00-01-02-plot-snp-boxplots.py` | Generates boxplots of SNP distribution across populations | PNG/PDF figures | Figure 2A |
| `00-01-03-pairwise-snp-tests.py` | Statistical comparisons (t-tests, ANOVA) of SNP counts between populations | p-values, effect sizes | Results text |
| `00-01-04-private-variant-count.py` | Identifies population-specific (private) variants | Variant lists | Figure 2B |
| `00-01-05-visualise-private-var-count.py` | Visualises private variant counts as bar charts | PNG/PDF figures | Figure 2B |
| `00-01-06-wgs-qc-metrics.py` | Calculates WGS quality control metrics (depth, Ti/Tv ratio, call rates) | QC summary table | Suppl. Table S1, Methods |

**Key findings**: Uros and Matzes show reduced heterozygosity consistent with demographic isolation; Ti/Tv ratio of 2.09 confirms high-quality variant calls.

---

### 00-02: VCF Annotation

| Script | Purpose | Output | Figure/Table |
|--------|---------|--------|--------------|
| `00-02-annotate-vcf-zygosity.py` | Annotates variants with zygosity status (homozygous/heterozygous) | Annotated VCF | Methods |

---

### 00-03: Population Structure

Scripts for IBD filtering, ancestry estimation, PCA, and population differentiation.

| Script | Purpose | Output | Figure/Table |
|--------|---------|--------|--------------|
| `00-03-00-common-variants-array-wgs.py` | Identifies overlapping variants between WGS and array platforms | Merged dataset (936,301 SNPs) | Methods |
| `00-03-01-remove-duplicates-IBD.py` | Removes duplicate/related samples using IBD (PI_HAT > 0.95) | 736 unrelated individuals | Table 2, Methods |
| `00-03-02-00-admixture-analysis.py` | Runs unsupervised ADMIXTURE for K=2-14 | Q matrices, CV errors | Figure 4, Suppl. Fig. S4 |
| `00-03-02-01-admixture-figures.py` | Generates ADMIXTURE bar plots and cross-validation plots | PNG/PDF figures | Figure 4, Suppl. Fig. S4 |
| `00-03-03-common-vars-peru-sdgp.py` | Merges Peruvian data with SGDP reference panel | Joint dataset for global PCA | Figure 3D |
| `00-03-04-pca-composite-analysis.py` | Performs PCA with regional/linguistic annotations | Multi-panel PCA figure | Figure 3A-D |
| `00-03-05-fst-pairwise-analysis.py` | Computes pairwise Weir & Cockerham Fst for all 28 populations | Fst matrix, heatmap, dendrogram | Suppl. Fig. S5 |
| `00-03-06-dstats-analysis.py` | Calculates outgroup f3 statistics using AFRODESCENDENTS as outgroup | f3 matrix, heatmap, dendrogram | Suppl. Fig. S6 |

**Key findings**:
- ADMIXTURE K=5 resolves Amazonian, Andean-Altiplano, other Andean, coastal/European-admixed, and African-enriched components
- PCA explains 51.5% (PC1), 9.4% (PC2), 6.7% (PC3) of variance
- Mean Fst = 0.014 (range 0.001-0.089); AFRODESCENDENTS most differentiated
- f3 statistics reveal tri-regional clustering: Amazonian (highest shared drift) → Highland → Coastal (lowest)

---

### 00-04: Pharmacogenomics

Scripts for star-allele calling and clinical interpretation using PyPGx.

| Script | Purpose | Output | Figure/Table |
|--------|---------|--------|--------------|
| `00-04-00-run-pypgx.py` | Runs PyPGx pipeline on all 736 samples across 56 pharmacogenes | Diplotypes, phenotypes | Figure 5, Methods |
| `00-04-01-integrated-pgx-figures.py` | Generates integrated pharmacogenomics figures (phenotype frequencies, FDA triggers, heatmaps) | Multi-panel Figure 5 | Figure 5A-C |
| `00-04-02-pgx-clinical-interpretation.py` | Maps genotypes to FDA recommendations and clinical guidelines | Drug-gene associations | Figure 5B, Results |

**Key findings**:
- CYP2C19 poor metaboliser frequency elevated in Indigenous populations
- CYP3A5 expresser status higher in Amazonian groups
- FDA guidelines may have limited applicability to Indigenous Peruvian populations

---

### 00-05: gnomAD Validation

Scripts for verifying novel variant claims against gnomAD.

| Script | Purpose | Output | Figure/Table |
|--------|---------|--------|--------------|
| `00-05-00-gnomad-filtered-check.py` | Verifies high-impact variants against gnomAD v4.1 including filtered variants | Validation report | Results text, R1.3 response |

**Key findings**: 94% of novel variants (1,536,198 of 1,638,862) confirmed absent from gnomAD v2.1.1.

---

### 00-06: Variant Annotation and Statistics

Scripts for gnomAD annotation and statistical testing.

| Script | Purpose | Output | Figure/Table |
|--------|---------|--------|--------------|
| `00-06-00-gnomad-download.py` | Downloads gnomAD reference files | Reference data | Methods |
| `00-06-01-gnomad-annotate.py` | Annotates variants with gnomAD allele frequencies | Annotated VCF | Table 1 |
| `00-06-00-02-gnomad-novel-crossref.py` | Cross-references novel variants with gnomAD and dbSNP | Novel variant counts | Table 1, Results |
| `00-06-01-pgx-statistical-tests.py` | Statistical comparisons of pharmacogenomic phenotypes between populations | Chi-square tests, p-values | Results text |

---

### 00-07: Figure Generation

Scripts for generating publication-ready figures.

| Script | Purpose | Output | Figure/Table |
|--------|---------|--------|--------------|
| `00-07-00-figure1-map-w-admixture-k5.py` | Creates Figure 1 map with sampling locations, regional shading, and ADMIXTURE K=5 pie charts | Figure 1B | Figure 1 |

**Features**:
- Geographic map of Peru with ecological region shading (Coastal, Highland, Amazonian)
- ADMIXTURE K=5 ancestry pie charts at each sampling location
- Circle size proportional to sample size
- South American inset map for context

---

## Dependencies

### Python (3.11+)

```
numpy
pandas
matplotlib
seaborn
scipy
geopandas
tqdm
```

### External Tools

| Tool | Version | Purpose |
|------|---------|---------|
| PLINK | 1.9 / 2.0 | Genotype QC, IBD, Fst |
| ADMIXTURE | 1.3.0 | Ancestry estimation |
| PyPGx | latest | Pharmacogenomics |
| bcftools | 1.17+ | VCF manipulation |
| VEP | 110 | Variant effect prediction |

### Reference Data

- gnomAD v2.1.1 (GRCh37)
- dbSNP (build 155)
- SGDP reference panel
- PharmVar database

---

## Reproducibility

### Hardware

All analyses were performed on:
- **Mac Studio M3 Ultra** (32 cores, 256 GB RAM)
- macOS Sonoma

Scripts are optimised for parallel processing (typically 28 workers).

### Reference Genome

GRCh37 (hg19) was retained for:
1. Compatibility with PyPGx and PharmVar database coordinates
2. Consistency with 1000 Genomes Phase 3 reference panels
3. Direct comparison with previous Peruvian Genome Project publications
4. Avoiding liftover errors in complex pharmacogene regions (e.g., CYP2D6)

### Running the Pipeline

```bash
# 1. Genetic diversity
python3.11 00-01-00-gen-div.py
python3.11 00-01-01-snp-counts.py
python3.11 00-01-02-plot-snp-boxplots.py

# 2. Population structure
python3.11 00-03-01-remove-duplicates-IBD.py
python3.11 00-03-02-00-admixture-analysis.py
python3.11 00-03-04-pca-composite-analysis.py
python3.11 00-03-05-fst-pairwise-analysis.py
python3.11 00-03-06-dstats-analysis.py

# 3. Pharmacogenomics
python3.11 00-04-00-run-pypgx.py
python3.11 00-04-01-integrated-pgx-figures.py

# 4. Figure generation
python3.11 00-07-00-figure1-map-w-admixture-k5.py
```

---

## Citation

If you use these scripts or data, please cite:

> Corpas M, Guio H, et al. (2026). Genomic Insights into 30 Peruvian Indigenous and Mestizo Populations. *Nature Human Behaviour* (in revision).

---

## Contact

- **Manuel Corpas** - Corresponding author
- **Heinner Guio** - Corresponding author

---

## License

This project is licensed under [LICENSE] - see the LICENSE file for details.
