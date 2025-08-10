import pandas as pd
from scipy.stats import mannwhitneyu
from itertools import combinations
import statsmodels.stats.multitest as multitest
import seaborn as sns
import matplotlib.pyplot as plt

# File path to the input CSV file
input_file = "ANALYSIS/00-01-GEN-DIV/individual_snp_counts.csv"

def pairwise_tests(input_file, alpha=0.05):
    """
    Perform pairwise statistical tests to compare SNP counts between groups.
    Uses the Mann-Whitney U test for pairwise comparisons.
    """
    # Read the SNP counts data
    df = pd.read_csv(input_file)

    # Check for necessary columns
    if not {"Group", "SNP Count"}.issubset(df.columns):
        raise ValueError("The input file must contain 'Group' and 'SNP Count' columns.")

    # Ensure SNP Count is numeric
    df["SNP Count"] = pd.to_numeric(df["SNP Count"], errors="coerce")
    df = df.dropna(subset=["SNP Count"])  # Remove rows with invalid/missing SNP counts

    # Get unique groups
    groups = df["Group"].unique()

    # Store pairwise results
    results = []
    p_value_matrix = pd.DataFrame(index=groups, columns=groups, dtype=float)

    # Perform pairwise comparisons
    for group1, group2 in combinations(groups, 2):
        # Extract SNP counts for the two groups
        data1 = df[df["Group"] == group1]["SNP Count"]
        data2 = df[df["Group"] == group2]["SNP Count"]

        # Perform Mann-Whitney U test
        stat, p_value = mannwhitneyu(data1, data2, alternative="two-sided")
        results.append((group1, group2, stat, p_value))

        # Fill matrix with raw p-values
        p_value_matrix.loc[group1, group2] = p_value
        p_value_matrix.loc[group2, group1] = p_value

    # Adjust p-values for multiple testing using the Bonferroni method
    p_values = [result[3] for result in results]
    adjusted_p_values = multitest.multipletests(p_values, method="bonferroni")[1]

    # Add adjusted p-values to results
    results = [
        (group1, group2, stat, p_value, adj_p)
        for (group1, group2, stat, p_value), adj_p in zip(results, adjusted_p_values)
    ]

    # Update matrix with adjusted p-values
    for (group1, group2, _, _, adj_p) in results:
        p_value_matrix.loc[group1, group2] = adj_p
        p_value_matrix.loc[group2, group1] = adj_p

    # Convert results to a DataFrame
    results_df = pd.DataFrame(
        results, columns=["Group 1", "Group 2", "Statistic", "P-Value", "Adjusted P-Value"]
    )

    # Mark significant results
    results_df["Significant"] = results_df["Adjusted P-Value"] < alpha

    return results_df, p_value_matrix

from scipy.cluster.hierarchy import linkage, leaves_list
import numpy as np
from matplotlib.colors import Normalize, LinearSegmentedColormap
def plot_clustered_heatmap(p_value_matrix, alpha=0.05):

    """
    Plot a heatmap of adjusted p-values where smaller (more significant) p-values are represented by more intense red shades.
    """
    import numpy as np
    from matplotlib.lines import Line2D

    # Replace NaN or infinite values with a large placeholder (1.0 for no significance)
    filled_matrix = p_value_matrix.fillna(1.0).replace([np.inf, -np.inf], 1.0)

    # Perform hierarchical clustering
    linkage_matrix = linkage(filled_matrix, method="average", optimal_ordering=True)
    ordered_indices = leaves_list(linkage_matrix)

    # Reorder matrix
    ordered_matrix = filled_matrix.iloc[ordered_indices, ordered_indices]

    # Create a mask for non-significant cells
    significance_mask = ordered_matrix >= alpha

    # Annotate with adjusted p-values (only for significant cells)
    annotation_matrix = ordered_matrix.applymap(
        lambda x: f"{x:.2e}*" if x < alpha else ""
    )

    # Set vmin as the smallest value and vmax as alpha
    vmin = ordered_matrix.min().min()
    vmax = alpha

    # Define a custom colormap from light to intense red
    reds_cmap = LinearSegmentedColormap.from_list("CustomReds", ["#ffe6e6", "#ff0000"])  # Light pink to dark red

    # Create the heatmap
    plt.figure(figsize=(12, 10))
    sns.heatmap(
        ordered_matrix,
        annot=annotation_matrix,
        fmt="",  # Display p-values only for significant cells
        cmap=reds_cmap,  # Custom reds colormap
        mask=significance_mask,  # Mask non-significant cells
        norm=Normalize(vmin=vmin, vmax=vmax),  # Correct normalization
        cbar=False  # Remove the color bar
    )

    # Set x-axis and y-axis ticks with proper alignment
    plt.xticks(ticks=np.arange(len(ordered_matrix.columns)) + 0.5, labels=ordered_matrix.columns, rotation=45, ha="right")
    plt.yticks(ticks=np.arange(len(ordered_matrix.index)) + 0.5, labels=ordered_matrix.index, rotation=0)

    # Add a detailed legend
    legend_elements = [
        Line2D([0], [0], marker="s", color="w", label="Smaller P-Values = Lighter Red", markerfacecolor="#ff0000", markersize=10),
        Line2D([0], [0], marker="s", color="w", label="Not Significant (p >= 0.05)", markerfacecolor="white", markersize=10),
        Line2D([0], [0], linestyle="None", marker="*", color="black", label="Star (*) = Adjusted P-Value"),
    ]
    legend = plt.legend(
        handles=legend_elements,
        loc="upper right",
        title="Legend",
        frameon=True,
    )

    # Customize legend background and border
    legend.get_frame().set_facecolor("white")
    legend.get_frame().set_edgecolor("black")
    legend.get_frame().set_alpha(1.0)

    # Add a title
    plt.title("Pairwise Differences in SNP Counts Between Populations", fontsize=16)
    plt.xlabel("Population", fontsize=12)
    plt.ylabel("Population", fontsize=12)
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    # Perform pairwise tests
    results_df, p_value_matrix = pairwise_tests(input_file)

    # Save the results to a CSV file
    output_file = "pairwise_snp_tests_results.csv"
    results_df.to_csv(output_file, index=False)

    # Print summary
    print("Pairwise comparison results:")
    print(results_df)
    print(f"\nResults saved to {output_file}")

    # Plot heatmap of adjusted p-values
    plot_clustered_heatmap(p_value_matrix)

