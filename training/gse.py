import pandas as pd

df = pd.read_csv(".data/processed/mp_logs_logp_1212_with_info.csv")

def calculate_gse_logS(melting_point, logP):
    return 0.5 - 0.01 * (melting_point - 25) - logP

df["LogS-GSE"] = df.apply(lambda row: calculate_gse_logS(row["MP-Measured"], row["LogP-Measured"]), axis=1)

# 必要な列だけでもOK
df[["LogS-GSE"]].to_csv(".models/logs_gse.csv", index=False)