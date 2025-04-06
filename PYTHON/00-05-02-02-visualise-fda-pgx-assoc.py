#!/usr/bin/env python3

"""
patient_and_population_recommendations.py

Generates:
  1) patient_fda_matches.csv (detailed table)
  2) patient_recommendations_bar.png (one bar per patient, sorted desc)
  3) population_recommendations_bar.png (one bar per population, sorted desc)

All data shown (no "Others" grouping).

REQUIRES:
  1) Genotype TSV (GENOTYPE_TSV) with columns:
       - "SampleID" (or "Unnamed: 0")
       - "Gene"
       - "Phenotype"
       - optionally "Population"
         If missing, we parse from SampleID using POP_MAP (uppercase codes).
  2) FDA CSV (FDA_CSV) with columns:
       - "Gene"
       - "Drug"
       - "Affected Subgroups+"
       - "Description of Gene-Drug Interaction"

USAGE:
  python3 patient_and_population_recommendations.py
"""

import os
import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# File paths (adjust if needed)
GENOTYPE_TSV = "ANALYSIS/00-05-PyPGx/all_genes_results.tsv"
FDA_CSV      = "ANALYSIS/00-05-PyPGx/fda_pharmacogenomics_associations.csv"
OUT_DIR      = "ANALYSIS/00-05-PyPGx"

# If your genotype data doesn't have a "Population" column,
# define uppercase codes found in SampleID => population name:
POP_MAP = {
    "CHOPCCAS": "Chopccas",
    "MOCHES":   "Moches",
    "TRUJILLO": "Trujillo",
    "IQUITOS":  "Iquitos",
    "UROS":     "Uros",
    "CUSCO":    "Cusco",
    "MATZES":   "Matzes"
}

def parse_population(sample_id):
    """Derive population from uppercase code in sample_id."""
    sid_up = sample_id.upper()
    for code, pop_name in POP_MAP.items():
        if code in sid_up:
            return pop_name
    return "Unknown"

def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    # ========== 1) Load genotype data ==========
    try:
        df_geno = pd.read_csv(GENOTYPE_TSV, sep="\t")  # If CSV, remove sep="\t"
    except Exception as e:
        print(f"Error reading genotype file {GENOTYPE_TSV}:\n{e}")
        sys.exit(1)

    if "Unnamed: 0" in df_geno.columns and "SampleID" not in df_geno.columns:
        df_geno.rename(columns={"Unnamed: 0": "SampleID"}, inplace=True)

    needed_geno = {"SampleID", "Gene", "Phenotype"}
    missing_g = needed_geno - set(df_geno.columns)
    if missing_g:
        print(f"ERROR: Genotype file missing columns: {missing_g}")
        sys.exit(1)

    # If no population column, parse from sample IDs
    if "Population" not in df_geno.columns:
        df_geno["Population"] = df_geno["SampleID"].apply(parse_population)

    # ========== 2) Load FDA data ==========
    try:
        df_fda = pd.read_csv(FDA_CSV)
    except Exception as e:
        print(f"Error reading FDA CSV {FDA_CSV}:\n{e}")
        sys.exit(1)

    needed_fda = {"Gene", "Drug", "Affected Subgroups+", "Description of Gene-Drug Interaction"}
    missing_fda = needed_fda - set(df_fda.columns)
    if missing_fda:
        print(f"ERROR: FDA file missing columns: {missing_fda}")
        sys.exit(1)

    # ========== 3) Merge + Filter for substring match ==========
    df_geno["phenotype_lower"] = df_geno["Phenotype"].str.lower()
    df_fda["subgroup_lower"]   = df_fda["Affected Subgroups+"].str.lower()

    df_merged = pd.merge(df_geno, df_fda, on="Gene", how="inner")

    def phenotype_matches_subgroup(phen, subgroup):
        return phen in subgroup  # substring check

    df_merged["MatchFlag"] = df_merged.apply(
        lambda row: phenotype_matches_subgroup(row["phenotype_lower"], row["subgroup_lower"]),
        axis=1
    )
    df_matched = df_merged[df_merged["MatchFlag"]].copy()

    if df_matched.empty:
        print("No FDA recommendations found. Possibly 'Indeterminate' phenotypes only.")
        sys.exit(0)

    # Make a final table of interest
    df_matched.rename(columns={"Description of Gene-Drug Interaction": "Recommendation"}, inplace=True)
    final_cols = ["SampleID", "Population", "Gene", "Phenotype", "Drug", "Recommendation"]
    df_final = df_matched[final_cols].drop_duplicates()

    # Save it
    out_csv = os.path.join(OUT_DIR, "patient_fda_matches.csv")
    df_final.to_csv(out_csv, index=False)
    print(f"Saved detailed matches to {out_csv}")

    # ========== 4) Build a palette for populations ==========
    unique_pops = df_final["Population"].unique().tolist()
    color_list = sns.color_palette("Set2", n_colors=len(unique_pops))
    pop_palette = dict(zip(unique_pops, color_list))

    # ========== 5) Patient-level bar chart (all patients) ==========
    rec_counts = df_final.groupby(["SampleID","Population"]).size().reset_index(name="NumRecommendations")
    # Sort descending
    rec_counts = rec_counts.sort_values("NumRecommendations", ascending=False)

    # This might get very long if you have many patients
    plt.figure(figsize=(10, max(6, 0.4*len(rec_counts))))
    sns.barplot(
        data=rec_counts,
        x="NumRecommendations",
        y="SampleID",
        hue="Population",
        palette=pop_palette,
        dodge=False
    )
    plt.title("FDA Recommendations by Patient (All Patients)")
    plt.xlabel("Number of Recommendations")
    plt.ylabel("Patient (SampleID)")
    plt.legend(title="Population", bbox_to_anchor=(1.02, 1.0), loc="upper left")
    plt.tight_layout()

    patient_bar = os.path.join(OUT_DIR, "patient_recommendations_bar.png")
    plt.savefig(patient_bar, dpi=300)
    plt.close()
    print(f"Saved patient-level bar chart to {patient_bar}")

    # ========== 6) Population-level bar chart ==========
    pop_counts = df_final.groupby("Population").size().reset_index(name="NumRecommendations")
    pop_counts = pop_counts.sort_values("NumRecommendations", ascending=False)

    plt.figure(figsize=(8, max(4, 0.4*len(pop_counts))))
    sns.barplot(
        data=pop_counts,
        x="NumRecommendations",
        y="Population",
        hue="Population",
        palette=pop_palette,
        dodge=False
    )
    plt.title("FDA Recommendations by Population")
    plt.xlabel("Number of Recommendations")
    plt.ylabel("Population")
    # Hide legend if desired (each row is a single pop)
    plt.legend(title="Population", bbox_to_anchor=(1.02, 1.0), loc="upper left")
    plt.tight_layout()

    pop_bar = os.path.join(OUT_DIR, "population_recommendations_bar.png")
    plt.savefig(pop_bar, dpi=300)
    plt.close()
    print(f"Saved population-level bar chart to {pop_bar}")

    print("All done! You have a full patient chart and a population chart in:", OUT_DIR)

if __name__ == "__main__":
    main()

