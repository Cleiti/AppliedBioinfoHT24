import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

#inputs
input_csv = "/proj/uppstore2019013/nobackup/private/1000spider_master_project/busco/busco_comparison_results.csv"
output_dir = "/proj/uppstore2019013/nobackup/private/1000spider_master_project/busco/plots"




df = pd.read_csv(input_csv)

#scatterplot: Completeness (Unfiltered vs Filtered)
plt.figure(figsize=(10, 6))
sns.scatterplot(x=df["Completeness_unfiltered"], y=df["Completeness_filtered"])
plt.title("Completeness: Unfiltered vs Filtered")
plt.xlabel("Completeness (Unfiltered)")
plt.ylabel("Completeness (Filtered)")
plt.savefig(os.path.join(output_dir, "scatter_completeness.png"))
plt.close()

# Scatterplot: Single Copy (Unfiltered vs Filtered)
plt.figure(figsize=(10, 6))
sns.scatterplot(x=df["Single_copy_unfiltered"], y=df["Single_copy_filtered"])
plt.title("Single-Copy BUSCOs: Unfiltered vs Filtered")
plt.xlabel("Single-Copy (Unfiltered)")
plt.ylabel("Single-Copy (Filtered)")
plt.savefig(os.path.join(output_dir, "scatter_single_copy.png"))
plt.close()

# Boxplot: Completeness Comparison
plt.figure(figsize=(10, 6))
df_melted_completeness = df.melt(
    id_vars=["Species"],
    value_vars=["Completeness_unfiltered", "Completeness_filtered"],
    var_name="Dataset",
    value_name="Completeness",
)
sns.boxplot(x="Dataset", y="Completeness", data=df_melted_completeness)
plt.title("Boxplot of Completeness (Unfiltered vs Filtered)")
plt.savefig(os.path.join(output_dir, "boxplot_completeness.png"))
plt.close()

# Boxplot: Single-Copy Comparison
plt.figure(figsize=(10, 6))
df_melted_single_copy = df.melt(
    id_vars=["Species"],
    value_vars=["Single_copy_unfiltered", "Single_copy_filtered"],
    var_name="Dataset",
    value_name="Single_Copy",
)
sns.boxplot(x="Dataset", y="Single_Copy", data=df_melted_single_copy)
plt.title("Boxplot of Single-Copy BUSCOs (Unfiltered vs Filtered)")
plt.savefig(os.path.join(output_dir, "boxplot_single_copy.png"))
plt.close()

# Violin Plot: Completeness
plt.figure(figsize=(10, 6))
sns.violinplot(x="Dataset", y="Completeness", data=df_melted_completeness, inner="quartile")
plt.title("Violin Plot of Completeness (Unfiltered vs Filtered)")
plt.savefig(os.path.join(output_dir, "violin_completeness.png"))
plt.close()ls


# Violin Plot: Single-Copy BUSCOs
plt.figure(figsize=(10, 6))
sns.violinplot(x="Dataset", y="Single_Copy", data=df_melted_single_copy, inner="quartile")
plt.title("Violin Plot of Single-Copy BUSCOs (Unfiltered vs Filtered)")
plt.savefig(os.path.join(output_dir, "violin_single_copy.png"))
plt.close()


summary_stats = df.describe()
summary_stats.to_csv(os.path.join(output_dir, "summary_statistics.csv"))

print(f"Visualizations and summary statistics saved to {output_dir}")