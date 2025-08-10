import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# File path to the CSV file
input_file = "ANALYSIS/00-01-GEN-DIV/individual_snp_counts.csv"

def plot_ordered_boxplots(input_file):
    """
    Reads the SNP counts from the CSV file and plots ordered boxplots for each group.
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
    sns.set(style="whitegrid")

    # Create a boxplot with ordered groups
    plt.figure(figsize=(10, 6))
    sns.boxplot(
        data=df,
        x="Group",
        y="SNP Count",
        order=group_order,  # Order groups by ascending median SNP count
        palette="Set2",
        showfliers=False
    )

    # Add a title and labels
    plt.title("Distribution of SNP Counts by Population (Ordered by Median)", fontsize=16)
    plt.xlabel("Population", fontsize=12)
    plt.ylabel("SNP Count", fontsize=12)
    plt.xticks(rotation=45)

    # Show the plot
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    # Plot the boxplots
    plot_ordered_boxplots(input_file)

