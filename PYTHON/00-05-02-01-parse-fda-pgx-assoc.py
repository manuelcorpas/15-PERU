import requests
from bs4 import BeautifulSoup
import pandas as pd

# URL of the FDA pharmacogenetic associations table
url = "https://www.fda.gov/medical-devices/precision-medicine/table-pharmacogenetic-associations"

# Fetch the webpage content
response = requests.get(url)
if response.status_code != 200:
    raise Exception(f"Failed to fetch URL: {url}")

# Parse the HTML content
soup = BeautifulSoup(response.text, 'html.parser')

# Locate the table (based on observation; adjust selector if needed)
table = soup.find('table')

# Parse table headers
headers = [header.text.strip() for header in table.find_all('th')]

# Parse table rows
rows = []
for row in table.find_all('tr')[1:]:  # Skip the header row
    cells = [cell.text.strip() for cell in row.find_all('td')]
    if cells:
        rows.append(cells)

# Convert to a pandas DataFrame
df = pd.DataFrame(rows, columns=headers)

# Save to CSV
output_file = "ANALYSIS/00-05-PyPGx/fda_pharmacogenetic_associations.csv"
df.to_csv(output_file, index=False)
print(f"Table saved to {output_file}")

