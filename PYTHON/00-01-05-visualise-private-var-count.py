#!/usr/bin/env python3

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def plot_snp_counts_and_private_variants(
    snp_counts_csv="ANALYSIS/00-01-GEN-DIV/individual_snp_counts.csv",
    private_variants_csv="ANALYSIS/00-01-GEN-DIV/private_variant_counts.csv",
    snp_counts_output="ANALYSIS/00-01-GEN-DIV/snp_counts_boxplot.png",
    private_variants_output="ANALYSIS/00-01-GEN-DIV/private_variant_barplot.png"
):
    """
    Creates two plots using a consistent color palette for each group:
      1) A boxplot of individual SNP counts (from 'SNP Count').
      2) A bar plot of total private variants (from 'PrivateVariantCount').

    CSV File Requirements:
      - snp_counts_csv: Must have columns: ['Group', 'Individual', 'SNP Count']
      - private_variants_csv: Must have columns: ['Group', 'PrivateVariantCount']

    Parameters:
        snp_counts_csv (str): Path to CSV with columns ['Group', 'Individual', 'SNP Count'].
        private_variants_csv (str): Path to CSV with columns ['Group', 'PrivateVariantCount'].
        snp_counts_output (str): Output path to save the SNP count boxplot (PNG).
        private_variants_output (str): Output path to save the private variant barplot (PNG).
    """

    # 1. Load the data
    df_snp = pd.read_csv(snp_counts_csv)       # columns: ['Group', 'Individual', 'SNP Count']
    df_priv = pd.read_csv(private_variants_csv)  # columns: ['Group', 'PrivateVariantCount']

    # 2. Combine all group names into one set, then sort
    all_groups = sorted(set(df_snp["Group"]).union(df_priv["Group"]))

    # 3. Create a consistent color palette
    palette = sns.color_palette("Set2", n_colors=len(all_groups))
    group_to_color = dict(zip(all_groups, palette))

    ################################################################
    # Plot 1: Boxplot of Individual SNP Counts
    ################################################################
    # Order groups by median 'SNP Count'
    group_order_snp = (
        df_snp.groupby("Group")["SNP Count"]
              .median()
              .sort_values()
              .index
              .tolist()
    )

    plt.figure(figsize=(8, 6))
    sns.set_style("whitegrid")

    # Create the boxplot
    sns.boxplot(
        x="Group",
        y="SNP Count",
        data=df_snp,
        order=group_order_snp,
        palette=[group_to_color[g] for g in group_order_snp]
    )

    plt.title("Distribution of SNP Counts per Individual")
    plt.xlabel("Group")
    plt.ylabel("SNP Count")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()

    if snp_counts_output:
        plt.savefig(snp_counts_output, dpi=300)
        print(f"Saved SNP count boxplot to: {snp_counts_output}")
    else:
        plt.show()

    ################################################################
    # Plot 2: Bar Chart of Total Private Variants
    ################################################################
    # Order groups by total 'PrivateVariantCount'
    group_order_priv = (
        df_priv.groupby("Group")["PrivateVariantCount"]
               .sum()
               .sort_values()
               .index
               .tolist()
    )

    plt.figure(figsize=(8, 6))
    sns.set_style("whitegrid")

    # Sort df_priv to match group_order_priv
    df_priv_sorted = df_priv.set_index("Group").loc[group_order_priv].reset_index()

    # Create the bar plot
    sns.barplot(
        x="Group",
        y="PrivateVariantCount",
        data=df_priv_sorted,
        palette=[group_to_color[g] for g in group_order_priv]
    )

    plt.title("Total Private Variants per Group")
    plt.xlabel("Group")
    plt.ylabel("Private Variant Count")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()

    if private_variants_output:
        plt.savefig(private_variants_output, dpi=300)
        print(f"Saved private variant bar chart to: {private_variants_output}")
    else:
        plt.show()

if __name__ == "__main__":
    # Update file paths if needed, then run:
    #   python3 00-01-05-visualise-private-var-count.py
    plot_snp_counts_and_private_variants()

