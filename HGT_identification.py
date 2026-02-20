import pandas as pd
import os
import numpy as np
import time

input_dir = "genomes_csv"
output_dir = "HGT_simple_classification"
log_file = os.path.join(output_dir, "HGT_analysis_log.txt")

os.makedirs(output_dir, exist_ok=True)

start_time = time.time()

with open(log_file, "w") as log:
    for file in os.listdir(input_dir):
        if not file.endswith("_CUB.csv"):
            continue

        df = pd.read_csv(os.path.join(input_dir, file))
        df = df[df["Quality_note"] == "OK"].copy() 

        thresholds = {
            "CAI": df["CAI"].quantile(0.10),
            "ENC": df["ENC"].quantile(0.90),
            "GC_zscore_self": df["GC_zscore_self"].abs().quantile(0.90),
            "Codon_Avg_Dist_self": df["Codon_Avg_Dist_self"].quantile(0.90),
        }

        df["suspect_CAI"] = df["CAI"] < thresholds["CAI"]
        df["suspect_ENC"] = df["ENC"] > thresholds["ENC"]
        df["suspect_GC"] = df["GC_zscore_self"].abs() > thresholds["GC_zscore_self"]
        df["suspect_RSCU"] = df["Codon_Avg_Dist_self"] > thresholds["Codon_Avg_Dist_self"]

        df["HGT_suspect_score"] = (
            df["suspect_CAI"].astype(int) +
            df["suspect_ENC"].astype(int) +
            df["suspect_GC"].astype(int) +
            df["suspect_RSCU"].astype(int)
        )

        df["HGT_predicted"] = df["HGT_suspect_score"] >= 2

        output_path = os.path.join(output_dir, file.replace("_CUB.csv", "_HGT_simple.csv"))
        df.to_csv(output_path, index=False)

        line = f"Processed: {file} -> {df['HGT_predicted'].sum()} potential HGT genes"
        print(line)
        log.write(line + "\n")

    total_time = time.time() - start_time
    summary = f"\nAnalysis complete in {total_time:.2f} seconds ({total_time / 60:.2f} minutes)"
    print(summary)
    log.write(summary + "\n")
