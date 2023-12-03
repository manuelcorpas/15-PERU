import pandas as pd

# This is to generate the ii_population for the PCA analysis, maximising the total number of genetic samples, 747 and their common variants.

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
file_path_excel = 'ANALYSIS/00-BIOBANK/peru-biobank-complete.xlsx'  # Replace with the actual path

file_to_output = 'ANALYSIS/10-PCA/ii_languages.txt'

with open(file_to_output, 'w') as file:
    print("IID\tFID\tPopulation", file=file)

# Reading and processing the Excel file
df_excel = pd.read_excel(file_path_excel)
for index, row in df_excel.iterrows():
    array_code = row['cod microarray']  # Replace with your actual column name for cod microarray
    genome_code = row['cod genoma']        
    if array_code in patterns:
        # Print pattern and '114 cual es la lengua' column value if there's a match
        lengua = str(row['114 cual es la lengua']).strip()  # Replace with your actual column name for '114 cual es la lengua'
        lengua.rstrip('.')
        if lengua == 'CASTELLANO' or lengua == '-9':
            continue
        if 'CASTELLANO' in lengua:
            lengua = lengua.replace('CASTELLANO', '')
        cleaned_elements = sorted([element.strip() for element in lengua.split(',') if element.strip()])
        lengua = ','.join(cleaned_elements).strip('.')
        column1 = array_code + '_' + array_code + '_' + array_code + '_' + array_code
        #print(column1,column1,lengua,sep='\t')
        
        with open(file_to_output, 'a') as file:
            print(column1, column1, lengua, sep='\t', file=file)

        #print(column1,column1,lengua,sep='\t',file=open(file_to_output,'a'))
    if genome_code in patterns:
        # Print pattern and '114 cual es la lengua' column value if there's a match
        lengua = str(row['114 cual es la lengua']).strip()  # Replace with your actual column name for '114 cual es la lengua'
        lengua.rstrip('.')
        if lengua == 'CASTELLANO' or lengua == '-9':
            continue
        if 'CASTELLANO' in lengua:
            lengua = lengua.replace('CASTELLANO', '')
        cleaned_elements = sorted([element.strip('.').strip() for element in lengua.split(',') if element.strip()])
        lengua = ','.join(cleaned_elements)
        column1 = genome_code + '_' + genome_code
        with open(file_to_output, 'a') as file:
            print(column1, column1, lengua, sep='\t', file=file)

        #print(column1,column1,lengua,sep='\t')
        #print(column1,column1,lengua,sep='\t',file=open(file_to_output,'a'))
    # Print patterns from the first file that were not matched
    #unmatched_patterns = patterns - matched_patterns
    #if unmatched_patterns:
        #for pattern in unmatched_patterns:
        #for index, row in df_excel.iterrows():
            #genome_code = row['cod genoma']
            #if genome_code in unmatched_patterns:
                #lengua  = row['114 cual es la lengua'].upper()
                #column1 = genome_code + '_' + genome_code
                #print(column1,column1,lengua,sep='\t',file=open(file_to_output,'a'))
    

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
