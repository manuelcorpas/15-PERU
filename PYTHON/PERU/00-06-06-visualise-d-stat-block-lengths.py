
import pyreadr
import matplotlib.pyplot as plt

# Function to read RDS file and return as Pandas DataFrame
def read_rds(file_path):
    result = pyreadr.read_r(file_path)
    return result[None]  # Extract the DataFrame from the result

# Function to plot histogram using matplotlib
def plot_histogram(data, title, xlab):
    plt.hist(data, bins='auto')
    plt.title(title)
    plt.xlabel(xlab)
    plt.ylabel('Frequency')

# Load RDS files
path = 'ANALYSIS/15-DI-STAT/D-Stat-out'
block_lengths_f2 = read_rds(f'{path}/block_lengths_f2.rds')
block_lengths_ap = read_rds(f'{path}/block_lengths_ap.rds')
block_lengths_fst = read_rds(f'{path}/block_lengths_fst.rds')

# Plotting
#plot_histogram(block_lengths_f2, 'Histogram of Block Lengths for F2', 'Block Lengths')
#plot_histogram(block_lengths_ap, 'Histogram of Block Lengths for AP', 'Block Lengths')
plot_histogram(block_lengths_fst, 'Histogram of Block Lengths for FST', 'Block Lengths')

# Show the plots
plt.show()

