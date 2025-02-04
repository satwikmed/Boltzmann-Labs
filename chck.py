import pandas as pd

# Read the CSV file
data = pd.read_csv("/Users/satwikmedipalli/single_cell/combined_data.csv")
print(data.dtypes)

# Loop through each column (except the first one) and attempt to convert to float64
for column in data.columns[1:]:
    data[column] = pd.to_numeric(data[column], errors='coerce')

# Save the updated DataFrame to a CSV file
data.to_csv("/Users/satwikmedipalli/single_cell/updated.csv", index=False)
print(data.dtypes)

print("File 'updated.csv' has been saved successfully.")


