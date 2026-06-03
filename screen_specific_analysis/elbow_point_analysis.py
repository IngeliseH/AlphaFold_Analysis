import pandas as pd
import numpy as np
from pathlib import Path
from kneed import KneeLocator
import matplotlib.pyplot as plt

csv_path = Path.home() / "Desktop" / "all_interface_analysis_2025.11.26_absolute.csv"
df = pd.read_csv(csv_path)

# IDENTIFYING AND FILTERING OUT SMALL INTERFACES
# drop all rows where size is na (those where a prediction had no interfaces)
df = df.dropna(subset=['size'])
size_counts = df['size'].value_counts().sort_index()
# find elbow point in size data density to determine cutoff for filtering out small interfaces
x = size_counts.index.values
y = size_counts.values

# initial plot to check curve shape
fig, ax = plt.subplots(figsize=(6, 6))
ax.plot(size_counts.index, size_counts.values, marker='o')
ax.set_xlabel('Interface size (residue pairs)')
ax.set_ylabel('Count')
ax.set_title('Interface Size Data Density')
plt.tight_layout()
plt.show()

kn = KneeLocator(x, y, curve='convex', direction='decreasing')
print(kn.knee)
#elbow point for size = 15

fig, ax = plt.subplots(figsize=(6, 6))
ax.plot(size_counts.index, size_counts.values, marker='o')
ax.axvline(x=kn.knee, color='red', linestyle='--', label=f'Cutoff (size = {kn.knee})')
ax.set_xlabel('Interface size (residue pairs)')
ax.set_ylabel('Count')
ax.set_title('Interface Size Data Density')
ax.legend()
plt.tight_layout()
plt.show()

# FILTER OUT SMALL INTERFACES USING ELBOW POINT
df = df[df['size'] >= kn.knee]

# FINDING ELBOW POINTS OF MIN AND AVG PAE
# Min and avg pae are not integers. Rounding these up, so a whole number threshold can be set for filtering
min_pae_binned = np.ceil(df['min_pae'].round(0).astype(int))
min_pae_value_counts = min_pae_binned.value_counts().sort_index()

avg_pae_binned = np.ceil(df['avg_pae'].round(0).astype(int))
avg_pae_value_counts = avg_pae_binned.value_counts().sort_index()

# initial plot to check curve shape for min and avg pae
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))

ax1.plot(min_pae_value_counts.index, min_pae_value_counts.values, marker='o')
ax1.set_xlabel('min_pae (binned)')
ax1.set_ylabel('Count')
ax1.set_title('Data Density')

ax2.plot(avg_pae_value_counts.index, avg_pae_value_counts.values, marker='o')
ax2.set_xlabel('avg_pae (binned)')
ax2.set_ylabel('Count')
ax2.set_title('Data Density')

plt.tight_layout()
plt.show()

# find elbows
# for min pae, only consider x values 1-14 (portion of curve with exponentially increasing shape - )
x_min_pae = min_pae_value_counts[(min_pae_value_counts.index >= 1) & (min_pae_value_counts.index <= 14)].index.values
y_min_pae = min_pae_value_counts[(min_pae_value_counts.index >= 1) & (min_pae_value_counts.index <= 14)].values
kn_min_pae = KneeLocator(x_min_pae, y_min_pae, curve='convex', direction='increasing')
print(f"Elbow point for min_pae: {kn_min_pae.knee}")
#Elbow point for min_pae: 11.0

# for avg pae, only consider x 1-21 (portion of curve with exponentially increasing shape)
x_avg_pae = avg_pae_value_counts[(avg_pae_value_counts.index >= 1) & (avg_pae_value_counts.index <= 21)].index.values
y_avg_pae = avg_pae_value_counts[(avg_pae_value_counts.index >= 1) & (avg_pae_value_counts.index <= 21)].values
kn_avg_pae = KneeLocator(x_avg_pae, y_avg_pae, curve='convex', direction='increasing')
print(f"Elbow point for avg_pae: {kn_avg_pae.knee}")
#Elbow point for avg_pae: 14.0

# Plotting with elbow points
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
ax1.plot(min_pae_value_counts.index, min_pae_value_counts.values, marker='o')
ax1.axvline(x=kn_min_pae.knee, color='red', linestyle='--', label=f'Elbow (min_pae = {kn_min_pae.knee})')
ax1.axvspan(14, min_pae_value_counts.index.max(), alpha=0.2, color='grey', label='Region not considered')
ax1.set_xlabel('min_pae (binned)')
ax1.set_ylabel('Count')
ax1.set_title('Data Density')
ax1.legend()

ax2.plot(avg_pae_value_counts.index, avg_pae_value_counts.values, marker='o')
ax2.axvline(x=kn_avg_pae.knee, color='red', linestyle='--', label=f'Elbow (avg_pae = {kn_avg_pae.knee})')
ax2.axvspan(21, avg_pae_value_counts.index.max(), alpha=0.2, color='grey', label='Region not considered')
ax2.set_xlabel('avg_pae (binned)')
ax2.set_ylabel('Count')
ax2.set_title('Data Density')
ax2.legend()

plt.tight_layout()
plt.show()

# FILTERING DATA USING ELBOW POINTS FOR MIN AND AVG PAE
df = df[(df['min_pae'] <= kn_min_pae.knee) & (df['avg_pae'] <= kn_avg_pae.knee)]

# FINDING ELBOW POINT OF ROP 
rop_value_counts = df['rop'].value_counts().sort_index()

# initial plot to check curve shape
fig, ax = plt.subplots(figsize=(6, 6))
ax.plot(rop_value_counts.index, rop_value_counts.values, marker='o')
ax.set_xlabel('rop')
ax.set_ylabel('Count')
ax.set_title('Data Density')
plt.tight_layout()
plt.show()

# not subsetting data for rop, as whole curve follows a standard exponentially decreasing pattern
x_rop = rop_value_counts.index.values
y_rop = rop_value_counts.values
kn_rop = KneeLocator(x_rop, y_rop, curve='convex', direction='decreasing')
print(f"Elbow point for rop: {kn_rop.knee}")
#Elbow point for rop: 1

# Plotting with elbow point
fig, ax = plt.subplots(figsize=(6, 6))
ax.plot(rop_value_counts.index, rop_value_counts.values, marker='o')
ax.axvline(x=kn_rop.knee, color='red', linestyle='--', label=f'Elbow (rop = {kn_rop.knee})')
ax.set_xlabel('rop')
ax.set_ylabel('Count')
ax.set_title('Data Density')
ax.legend()
plt.tight_layout()
plt.show()
