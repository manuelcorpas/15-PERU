import pandas as pd

# Function to find repeating patterns (from the first file)
def find_repeating_pattern_advanced(s, separator='_'):
    """
    Function to find the repeating pattern in a string, considering both continuous and non-continuous repetitions.
    It looks for the longest pattern that repeats and prints it once.
    """
    parts = s.split(separator)

    # Iterate over possible pattern lengths
    for length in range(1, len(parts)):
        for start in range(len(parts) - length + 1):
            pattern = separator.join(parts[start:start + length])
            # Check if the pattern is repeating
            if s.count(pattern) > 1:
                # To ensure that the pattern is the repeating unit and not a subset of a larger repeating pattern
                if s.replace(pattern, '') == separator * (s.count(pattern) - 1):
                    return pattern

    return None

# Path to your first file (with repeating patterns)
file_path_patterns = 'ANALYSIS/10-PCA/common_variants_ibd_samples.txt'  # Replace with the actual path

# Reading and storing patterns from the first file
patterns = set()
with open(file_path_patterns, 'r') as file:
    for line in file:
        line = line.strip()
        pattern = find_repeating_pattern_advanced(line)
        if pattern:
            patterns.add(pattern)

# Path to your Excel file (second file)
file_path_excel = 'ANALYSIS/03-SAMPLE-ANALYSIS/DemographicData_INS.xlsx'  # Replace with the actual path

# Reading and processing the Excel file
matched_patterns = set()
try:
    df_excel = pd.read_excel(file_path_excel)
    for index, row in df_excel.iterrows():
        array_code = row['ARRAY_CODE']  # Replace with your actual column name for ARRAY_CODE
        if array_code in patterns:
            matched_patterns.add(array_code)
            # Print pattern and '107 etnia' column value if there's a match
            try:
                etnia_value = row['POPULATION_CODE'].upper()  # Replace with your actual column name for '107 etnia'
                if etnia_value == 'MATSES':
                    etnia_value = 'MATZES'
                if etnia_value == 'MOCHE':
                    etnia_value = 'MOCHES'
                if '_' in array_code:
                      column1 = array_code + '_' + array_code + '_' + array_code + '_' + array_code
                print(column1, etnia_value,sep='\t')
            except KeyError:
                # If '107 etnia' column is not found, print just the pattern
                #print()
                x=0
        else:
            # If no pattern is found, print just the array code
            #print()
            x=0
    
    # Print patterns from the first file that were not matched
    unmatched_patterns = patterns - matched_patterns
    if unmatched_patterns:
        for pattern in unmatched_patterns:
            if '-' in pattern:
                etnia_value = pattern.split('-')[0]
                column1     = pattern + '_' + pattern
                print(column1,etnia_value,sep='\t')
    
except Exception as e:
    print("Error reading Excel file:", e)




'''
def find_repeating_pattern_advanced(s, separator='_'):
    """
    Function to find the repeating pattern in a string, considering both continuous and non-continuous repetitions.
    It looks for the longest pattern that repeats and prints it once.
    """
    parts = s.split(separator)

    # Iterate over possible pattern lengths
    for length in range(1, len(parts)):
        for start in range(len(parts) - length + 1):
            pattern = separator.join(parts[start:start + length])
            # Check if the pattern is repeating
            if s.count(pattern) > 1:
                # To ensure that the pattern is the repeating unit and not a subset of a larger repeating pattern
                if s.replace(pattern, '') == separator * (s.count(pattern) - 1):
                    return pattern

    return None

# Path to your file
file_path = 'ANALYSIS/10-PCA/common_variants_ibd_samples.txt'

# Reading and processing the file
with open(file_path, 'r') as file:
    for line in file:
        line = line.strip()
        pattern = find_repeating_pattern_advanced(line)
        if pattern:
            print(pattern)


'''
