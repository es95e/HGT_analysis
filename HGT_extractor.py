import os
import pandas as pd
import time

input_folder = 'HGT_simple_classification'
output_folder = 'filtered_results'
log_file = os.path.join(output_folder, "filtered_results_log.txt")

os.makedirs(output_folder, exist_ok=True)

start_time = time.time()

with open(log_file, "w") as log:
    for filename in os.listdir(input_folder):
        if filename.endswith('.csv'):
            file_path = os.path.join(input_folder, filename)

            df = pd.read_csv(file_path)

            if 'HGT_predicted' in df.columns:
                df_filtered = df[df['HGT_predicted'].astype(str).str.lower() == 'true']
                output_filename = filename.replace('.csv', '_filtered.csv')
                output_path = os.path.join(output_folder, output_filename)
                df_filtered.to_csv(output_path, index=False)

                line = f"Processed: {filename} -> {len(df_filtered)} row extracted"
            else:
                line = f"Coloum 'HGT_predicted' not found {filename}"

            print(line)
            log.write(line + "\n")

combined_df = pd.DataFrame()

for filtered_file in os.listdir(output_folder):
    if filtered_file.endswith('_filtered.csv'):
        filtered_path = os.path.join(output_folder, filtered_file)
        df_filtered = pd.read_csv(filtered_path)
        df_filtered['source_file'] = filtered_file 
        combined_df = pd.concat([combined_df, df_filtered], ignore_index=True)

combined_output_path = os.path.join(output_folder, 'all_filtered_combined.csv')
combined_df.to_csv(combined_output_path, index=False)

line = f"\n all_filtered_combined.csv with {len(combined_df)} total rows"
print(line)
with open(log_file, "a") as log:
    log.write(line + "\n")

total_time = time.time() - start_time
summary = f"\n Filtering completed in {total_time:.2f} seconds ({total_time / 60:.2f} minutes)"
print(summary)
with open(log_file, "a") as log:
    log.write(summary + "\n")

