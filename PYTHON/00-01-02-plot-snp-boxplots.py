import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# File path to the CSV file
input_file = "ANALYSIS/00-01-GEN-DIV/individual_snp_counts.csv"
output_png = "ANALYSIS/00-01-GEN-DIV/Figure2A_SNP_boxplots.png"
output_pdf = "ANALYSIS/00-01-GEN-DIV/Figure2A_SNP_boxplots.pdf"

def plot_ordered_boxplots(input_file, save_png=None, save_pdf=None):
    """
    Reads the SNP counts from the CSV file and plots ordered boxplots for each group.
    
    Updated to address reviewer comment: "Please show all data points in the graphs"
    - Added jittered stripplot overlay showing individual data points
    - Added sample size annotations
    """
    # Read the CSV file into a pandas DataFrame
    df = pd.read_csv(input_file)

    # Check the structure of the data
    if not {"Group", "Individual", "SNP Count"}.issubset(df.columns):
        raise ValueError("The input file must contain 'Group', 'Individual', and 'SNP Count' columns.")

    # Convert SNP Count to numeric if not already
    df["SNP Count"] = pd.to_numeric(df["SNP Count"], errors="coerce")

    # Drop rows with invalid or missing SNP counts
    df = df.dropna(subset=["SNP Count"])

    # Calculate the median SNP count for each group
    group_order = df.groupby("Group")["SNP Count"].median().sort_values().index

    # Set the style of the plot
    sns.set_style("whitegrid")

    # Create figure
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Create a boxplot with ordered groups
    sns.boxplot(
        data=df,
        x="Group",
        y="SNP Count",
        order=group_order,
        palette="Set2",
        showfliers=False,  # Hide default outlier markers since we show all points
        width=0.6,
        linewidth=1.5,
        ax=ax
    )

    # Overlay individual data points using stripplot (addresses reviewer request)
    sns.stripplot(
        data=df,
        x="Group",
        y="SNP Count",
        order=group_order,
        color="black",
        alpha=0.5,
        size=5,
        jitter=0.15,
        ax=ax
    )

    # Add sample size annotations above each boxplot
    for i, group in enumerate(group_order):
        n = len(df[df["Group"] == group])
        max_val = df[df["Group"] == group]["SNP Count"].max()
        ax.annotate(f'n={n}', xy=(i, max_val), xytext=(0, 5),
                   textcoords='offset points', ha='center', va='bottom',
                   fontsize=9, color='dimgray')

    # Add a title and labels
    ax.set_title("Distribution of SNP Counts by Population (Ordered by Median)", fontsize=14)
    ax.set_xlabel("Population", fontsize=12)
    ax.set_ylabel("SNP Count", fontsize=12)
    plt.xticks(rotation=45, ha='right')

    # Adjust layout
    plt.tight_layout()

    # Save figures
    if save_png:
        plt.savefig(save_png, dpi=300, bbox_inches='tight')
        print(f"Saved: {save_png}")
    
    if save_pdf:
        plt.savefig(save_pdf, format='pdf', bbox_inches='tight')
        print(f"Saved: {save_pdf}")

    # Show the plot
    plt.show()
    
    return fig, ax

if __name__ == "__main__":
    # Plot the boxplots with individual data points
    plot_ordered_boxplots(input_file, save_png=output_png, save_pdf=output_pdf)