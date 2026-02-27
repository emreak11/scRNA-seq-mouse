from pathlib import Path
import numpy
import pandas as pd
import re

# Resolve data path relative to this script so it works when run from any cwd
data_path = Path(__file__).resolve().parent / 'geneid.txt'

# Read file (CSV with header). Keep empty strings as-is.
df = pd.read_csv(data_path, dtype=str, keep_default_na=False)

# Locate a column that looks like the gene name column (fallback to last column)
col = next((c for c in df.columns if 'name' in c.lower()), df.columns[-1])

# Locate a column that looks like gene stable ID
id_col = next((c for c in df.columns if 'stable id' in c.lower()), None)

# Regex: gene name starts with mt- or mt. (case-insensitive)
pat = re.compile(r'(?i)^mt[.-]')

# Select matching rows and print the gene name token(s)
mask = df[col].str.match(pat, na=False)
matches = df[mask]
if matches.empty:
    print('No matches found')
else:
    for _, row in matches.iterrows():
        if id_col:
            print(row[id_col])
        else:
            print(row[col])

           



if id_col:
    mt_id = {id_col: matches[id_col].tolist(), col: matches[col].tolist()}
else:
    mt_id = {col: matches[col].tolist()}

# Check number of matched rows (not number of dict keys)
#if len(matches) > 35:
    #print("works")

data=pd.read_fwf('mt_id_list.txt')
#if data is not None and not data.empty:
    #print("works")
    #print(len(data))

counts_path = Path("counts_raw.txt")  
counts = pd.read_csv(counts_path, sep="\t", index_col=0)
counts = counts.apply(pd.to_numeric, errors="coerce").fillna(0).astype("int64")



print("shape:", counts.shape)
print("index head:", counts.index[:3].tolist())
print("columns head:", counts.columns[:3].tolist())

mt_ids = (
    pd.read_csv("mt_id_list.txt", header=None, dtype=str)
      .iloc[:, 0]
      .str.strip()
)
mt_ids = mt_ids[mt_ids != ""].unique().tolist()

print("mt_ids:", len(mt_ids), "first 5:", mt_ids[:5])

mt_set = set(mt_ids)
idx_set = set(counts.index.astype(str))

overlap = mt_set & idx_set
print("overlap:", len(overlap), "/", len(mt_set))

mt_counts = counts.loc[counts.index.intersection(mt_ids)]
print("mt_counts shape:", mt_counts.shape)

mt_sum_per_cell = mt_counts.sum(axis=0)
total_sum_per_cell = counts.sum(axis=0)

mt_frac = mt_sum_per_cell / total_sum_per_cell

print(mt_sum_per_cell.describe())
print(mt_frac.describe())

import matplotlib.pyplot as plt

#plt.hist(mt_frac, bins=50)
#plt.xlabel("Fraction of mitochondrial genes")
#plt.ylabel("Number of cells")
#plt.title("Distribution of mitochondrial gene fractions")
#plt.show()

# Şüpheli hücreler
#suspect = mt_frac > 0.065

#print("# of suspects:", suspect.sum())
#print("Summary of these cells:")
#print(total_sum_per_cell[suspect].describe())


thereshold = 0.065
good_cells = mt_frac <= thereshold
counts_qc = counts.loc[:, good_cells]

import numpy as np

normalized= counts_qc.div(counts_qc.sum(axis=0), axis=1) * 10**4
log_transform =np.log(1 + normalized)

def visualize():
    import matplotlib.pyplot as plt
    plt.subplot(1, 3, 1)
    plt.hist(mt_frac, bins=50)
    plt.axhline(thereshold, color='red', linestyle='--')
    plt.xlabel("Fraction of mitochondrial genes")
    plt.ylabel("Number of cells")
    plt.title("Distribution of mitochondrial gene fractions")
    plt.grid()
    
    plt.subplot(1, 3, 2)
    plt.scatter(total_sum_per_cell, mt_frac)
    plt.axvline(thereshold, color='red', linestyle='--')
    plt.xlabel("Mitochondrial fraction")
    plt.ylabel("Number of cells")
    plt.title("Mitochondrial fraction vs total counts per cell")
    plt.grid()

visualize()

def plot_log_vs_normalized(normalized, log_transform, n_samples=100000):
    # 1) Matrisi 1D hale getir
    norm_flat = normalized.to_numpy().ravel()
    log_flat  = log_transform.to_numpy().ravel()

    # 2) Rastgele aynı indeksleri seç
    idx = np.random.choice(len(norm_flat), size=n_samples, replace=False)

    norm_sample = norm_flat[idx]
    log_sample  = log_flat[idx]

    # 3) Histogramları üst üste çiz
    plt.figure(figsize=(6,4))
    plt.hist(norm_sample, bins=50, alpha=0.5, label="Normalized")
    plt.hist(log_sample,  bins=50, alpha=0.5, label="Log-normalized")

    plt.xlabel("Expression value")
    plt.ylabel("Frequency")
    plt.title("Normalized vs Log-normalized expression distribution")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.show()
plot_log_vs_normalized(normalized, log_transform)


