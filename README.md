

# Scripts Documentation for the Paper:
**Genome Interpretation of Peruvian Inca Empire Descendants Reveals Actionable Insights**

In our manuscript, we describe multiple analyses on whole-genome sequencing (WGS) and SNP array data from 30 Peruvian populations. The Python scripts in `PYTHON/` were pivotal for generating the results, figures, and tables appearing in the paper. Below, we outline each script’s role and where its outputs appear in the manuscript.

## 1. Genetic Diversity and SNP Analysis

### 00-01-00-gen-div.py
- **Purpose in the Paper**: Computes general genetic diversity metrics (e.g., heterozygosity, allele richness) used in the *“Whole-Genome Mutational Landscape”* and *“Genetic Diversity and Structure”* sections.  
- **Relevant Figures/Results**: Contributed to summary statistics of biallelic SNPs across populations, described in the text near *Figure 2* (SNP distributions) and *Table 1* (SNV counts).

### 00-01-01-snp-counts.py
- **Purpose in the Paper**: Tallies SNPs per sample or population, providing the raw counts featured in the *“Whole-Genome Mutational Landscape”* subsection.  
- **Relevant Figures/Results**: Underlies *Figure 2A*, which presents the boxplot of SNP counts per individual.

### 00-01-02-plot-snp-boxplots.py
- **Purpose in the Paper**: Generates boxplots of SNP distribution, helping visualize how many variants each individual harbors.  
- **Relevant Figures/Results**: Directly used to create the boxplot shown in *Figure 2A*, illustrating how certain populations (e.g., Uros, Matzes) have fewer SNPs compared to more admixed groups.

### 00-01-03-pairwise-snp-tests.py
- **Purpose in the Paper**: Performs statistical comparisons (e.g., t-tests, ANOVA) across populations to detect significant differences in SNP counts.  
- **Relevant Figures/Results**: Though not highlighted in a standalone figure, its p-values and effect sizes support statements in the *“Whole-Genome Mutational Landscape”* discussing variation among populations.

### 00-01-04-private-variant-count.py
- **Purpose in the Paper**: Identifies private (population-specific) variants in each group. This analysis informed the text about unique bottlenecks and founder effects.  
- **Relevant Figures/Results**: Ties into *Figure 2B* (bar chart of private variants per population). Also referenced in the *“High-Impact Variant Discovery”* section when discussing population-specific alleles.

### 00-01-05-visualise-private-var-count.py
- **Purpose in the Paper**: Produces plots (barplots/boxplots) of private variant counts.  
- **Relevant Figures/Results**: Contributed to the final version of *Figure 2B*, which highlights how certain groups (e.g., Trujillo, Cusco) contain higher unique variant counts than more isolated populations (e.g., Matzes).

---

## 2. VCF Annotation and Zygosity

### 00-03-annotate-vcf-zygosity.py
- **Purpose in the Paper**: Annotates each variant with zygosity (homozygous/heterozygous), supporting the *“Whole-Genome Mutational Landscape”* methods.  
- **Relevant Figures/Results**: Used internally during data preparation (not visualized directly as a figure), but essential to the final curated dataset described in *“Online Methods – Variant Calling and Quality Control.”*

---

## 3. PharmCAT Pipeline

### 00-04-00-pharmcat-pipeline.py
- **Purpose in the Paper**: Runs PharmCAT on multiple WGS or VCF files, generating `.report.html` outputs used for pharmacogenomic profiling.  
- **Relevant Figures/Results**: Feeds into the *“Pharmacogenomic Analysis”* section. Data here were ultimately referenced in *Figure 3* (FDA recommendations per population).

### 00-04-01-parse-pharmcat-results.py
- **Purpose in the Paper**: Reads the PharmCAT `.report.html` files, extracting gene-allele-phenotype information.  
- **Relevant Figures/Results**: The summarized gene-phenotype tables appear in the *“Pharmacogenomic Analysis”* subsection, helping to quantify how frequently certain metabolizer statuses occur (e.g., poor metabolizers in CYP2C19).

---

## 4. PyPGx and FDA PGx Integration

### 00-05-00-run-pypgx.py
- **Purpose in the Paper**: Launches PyPGx to compute star-allele haplotypes, diplotypes, and predicted metabolizer phenotypes.  
- **Relevant Figures/Results**: Key for *“Pharmacogenomic Analysis.”* The outputs assisted in identifying clinically actionable variants (e.g., in CYP2C19, CYP2D6) discussed in *Figure 3* and text.

### 00-05-01-compile-pypgx-results.py
- **Purpose in the Paper**: Aggregates multiple PyPGx result files into a single summary table of genotype-phenotype calls.  
- **Relevant Figures/Results**: Facilitated cross-population comparisons in the *“Pharmacogenomic Analysis”* section, where differences in metabolizer status among the seven main populations were noted.

### 00-05-02-00-match-fda-recommendation.py
- **Purpose in the Paper**: Matches local genotypes/diplotypes to official FDA-labeled pharmacogenetic recommendations, highlighting potential dosing or warning flags.  
- **Relevant Figures/Results**: Populates the data seen in *Figure 3*, where each population’s gene–drug match count is displayed (e.g., “CYP2C19–clopidogrel poor metabolizer”).

### 00-05-02-01-parse-fda-pgx-assoc.py
- **Purpose in the Paper**: Parses FDA’s official pharmacogenomic tables or guidelines to create an internal reference for drug-gene associations.  
- **Relevant Figures/Results**: Used in combination with the above script to generate the *FDA Recommendation Mapping* dataset behind *Figure 3*.

### 00-05-02-02-visualise-fda-pgx-assoc.py
- **Purpose in the Paper**: Produces charts illustrating how many relevant gene–drug associations are triggered in each population.  
- **Relevant Figures/Results**: Directly creates the bar chart of FDA guidelines (shown in *Figure 3*), illustrating that certain subgroups have disproportionately high risk of drug–gene interactions.

---

## 5. Genotype Frequencies by Gene

### 00-05-03-genotype-freq-by-gene.py
- **Purpose in the Paper**: Calculates genotype frequencies for each gene across the study cohort.  
- **Relevant Figures/Results**: Supports the text discussing population-level differences in allele and genotype frequency (e.g., *“Genetic Diversity and Structure”* and part of *“Pharmacogenomic Analysis”*).

### 00-05-04-plot-predicted-pheno-by-gene.py
- **Purpose in the Paper**: Visualizes predicted phenotypes (e.g., poor/intermediate metabolizer) by gene across individuals.  
- **Relevant Figures/Results**: Could be adapted for supplementary figures illustrating distribution of metabolizer statuses or used to cross-check main text statements about prevalence.

### 00-05-05-genotype-freq-by-gene-and-population.py
- **Purpose in the Paper**: Extends the above script to show genotype frequencies *per population,* capturing how certain groups have unique allele distributions.  
- **Relevant Figures/Results**: Underpins statements about how “coastal mestizos” vs. “Andean highlanders” differ in genotype frequencies, further supporting *Figure 4A–C* on genetic structure and *Figure 3* on PGx differences.

---

## How These Scripts Integrate with the Paper

1. **Whole-Genome and SNP Analysis**:  
   Scripts `00-01-xx` create the variant count data and plots featured in *Figure 2* (SNP distribution and private variants).  
2. **PharmCAT & PyPGx**:  
   Scripts `00-04-xx` and `00-05-xx` drive the pharmacogenomic pipeline, parsing `.report.html` files (PharmCAT) and star-allele calls (PyPGx). Results feed into *Figure 3* (FDA recommendations) and text.  
3. **Population Genetics**:  
   Zygosity annotation (`00-03-annotate-vcf-zygosity.py`) and genotype frequencies (`00-05-03` through `00-05-05`) help characterize each group’s genomic profile, linking to PCA-based discussions (*Figure 4*) and high-impact variants (*Table 1*).

These scripts collectively underpin the manuscript’s core findings: the discovery of novel/rare variants, the identification of actionable pharmacogenomic markers, and the characterization of genetic structure in Peruvian populations. By integrating them with complementary tools (GATK, PLINK, VEP, etc.), we provide a reproducible framework for analyzing large genomic datasets and highlighting clinical implications for underserved communities.



