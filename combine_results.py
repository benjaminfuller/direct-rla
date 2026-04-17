import pandas as pd
import os
import glob
import sys

mode = sys.argv[1]
# Get all CSV files from the current directory
csv_files = [f for f in glob.glob('*.csv') if any(num in f for num in ['456000', '1000000', '1679423', '5007954', '7963659', '16140044']) and str(mode) in f]

# Read and concatenate all CSV files
dataframes = [pd.read_csv(file) for file in csv_files]
combined_df = pd.concat(dataframes, ignore_index=True)



# Write to output file
combined_df.to_csv('results_combined'+str(mode)+'.csv', index=False)

for csv_file in csv_files:
    os.remove(csv_file)

print(f"Combined {len(csv_files)} CSV files into results_combined{mode}.csv")