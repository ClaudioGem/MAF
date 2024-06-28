### Load table with pandas
import os
import numpy as np
import pandas as pd
os.chdir("/home/genomics/genomics09/Project/project5")
Variants_1 = pd.read_csv("snp_pvalue_JB1169_ed.txt", sep="\t")
Variants_1.columns
#Variants_1.shape
#print(Variants_1.head(10))
#Variants_1['Chrom'].unique()

### Filter tables including only variants present in chromosome I, II and III

file_names = ["snp_pvalue_JB1169_ed.txt", "snp_pvalue_JB1207_ed.txt", "snp_pvalue_JB22_ed.txt", "snp_pvalue_JB760_ed.txt"]

filtered_variants = [] #store filtered lines
chromosomes_of_interest = ["I", "II", "III"]

for i in range(len(file_names)):
    df = pd.read_csv(file_names[i], sep="\t") #load each file to df
    filtered_df = df[df["Chrom"].isin(chromosomes_of_interest)] #filter rows with chromosomes of interest
    filtered_variants.append(filtered_df)
print(filtered_variants)


### Merge data from all samples to create a bar plot with the number of variants per sample

import matplotlib.pyplot as plt
colors = ["steelblue", "sandybrown", "darkseagreen", "palevioletred"]

variant_counts = []

for file_name in file_names:
    df = pd.read_csv(file_name, sep="\t")
    filtered_df = df[df["Chrom"].isin(chromosomes_of_interest)]    
    variant_counts.append(len(filtered_df))
print(variant_counts)

plt.figure(figsize=(15, 10))
plt.bar(file_names, variant_counts, color=colors, edgecolor="black")
plt.ylabel("Number of Variants", fontsize=24)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.tight_layout()
plt.show()

### Create a function that transform the percentage value into minor allele frequency
###Use this function to create a new column in your data frame with
###minor allele frequency.

transformed_variance = []

def maf_from_varfreq(row):
    var_freq_raw = row["VarFreq"].rstrip('%')
    var_freq_percentage = float(var_freq_raw)
    if var_freq_percentage >= 50: #Transform VarFreq to MAF
        var_freq_percentage = 100 - var_freq_percentage
    maf = var_freq_percentage / 100
    return maf

for file_name in file_names:
    df = pd.read_csv(file_name, sep="\t")
    filtered_df = df[df["Chrom"].isin(chromosomes_of_interest)].copy()
    filtered_df.loc[:, "MAF"] = filtered_df.apply(maf_from_varfreq, axis=1)
    transformed_variance.append(filtered_df)
    output_file = f"filtered_with_maf_{file_name}"
    filtered_df.to_csv(output_file, sep="\t", index=False)
    print(f"Processed {file_name}, results saved to {output_file}")

print(transformed_variance)

### Use the MAF to plot the distribution per sample
import seaborn as sns
filtered_file_names = ["filtered_with_maf_snp_pvalue_JB1169_ed.txt", "filtered_with_maf_snp_pvalue_JB1207_ed.txt", "filtered_with_maf_snp_pvalue_JB22_ed.txt", "filtered_with_maf_snp_pvalue_JB760_ed.txt"]

for i, file_name in enumerate(filtered_file_names): # Read files, extract MAF values
    sample_name = file_name.split('_')[-2]  # sample name from file name
    df = pd.read_csv(file_name, sep='\t')
    
    plt.figure(figsize=(12, 10))
    sns.histplot(df["MAF"], bins=30, color=colors[i])
    plt.title(f"Distribution of MAF for {sample_name}", fontsize=22)
    plt.xlabel("MAF", fontsize=22)
    plt.ylabel("Counts",fontsize=22)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    if True: #only plot no saving
        plt.tight_layout()
        plt.show()
    if False: #Save plot to file        
        output_plot = f"maf_distribution_{sample_name}.png"
        plt.savefig(output_plot)
        plt.close()
        print(f"Plot saved for {sample_name} to {output_plot}")


### Set a windows size variable. For instance, 40000 bp
###Write a loop that goes along the genome by sections of your windows
###size (40000 bp) and within each window calculate the mean MAF.
###Plot values per window along the genome with a line plot

filtered_file_names = ["filtered_with_maf_snp_pvalue_JB1169_ed.txt", "filtered_with_maf_snp_pvalue_JB1207_ed.txt", "filtered_with_maf_snp_pvalue_JB22_ed.txt", "filtered_with_maf_snp_pvalue_JB760_ed.txt"]

window_size = 40000
colors = ["steelblue", "sandybrown", "darkseagreen", "palevioletred"]
chromosomes_of_interest = ["I", "II", "III"]

for i, file_name in enumerate(filtered_file_names):
    sample_name = file_name.split("_")[-2]
    df = pd.read_csv(file_name, sep="\t")
    
    for chrom in chromosomes_of_interest:
        chrom_df = df[df["Chrom"] == chrom].copy()
        chrom_df["Window"] = (chrom_df["Position"] // window_size) * window_size #calculate start position for each window using floor division
        window_means = chrom_df.groupby("Window").agg({"MAF": np.mean}).reset_index() #group by Window and aggregate MAF means for respective window

        plt.figure(figsize=(12, 10))
        sns.lineplot(data=window_means, x="Window", y="MAF", color=colors[i])
        plt.title(f"Mean MAF per {window_size} bp Window for {sample_name} - Chromosome {chrom}",fontsize=22)
        plt.xlabel("Genomic Position (Window Start)", fontsize=22)
        plt.ylabel("Mean MAF", fontsize=22)

        xticks = window_means["Window"][::20]  # every 20th window position
        xlabels = [f"{pos:,}" for pos in xticks]
        plt.xticks(ticks=xticks, labels=xlabels, fontsize=16)
        plt.yticks(fontsize=16)

        plt.tight_layout()
        plt.show()


### Your look in (i) will only consider variant sites. Improve point (i) by considering non-variant sites. 
###For this creating a vector of length 40000 (window size) with zeros, and replace zeros for the MAF value in variant sites.
###Calculate a new mean MAF per window. Produce new plot of mean MAF per window along the genome

window_size = 40000

for i, file_name in enumerate(file_names): #read files
    sample_name = file_name.split("_")[-2] 
    df = pd.read_csv(file_name, sep="\t")
    
    for chrom in chromosomes_of_interest: #select chromosomes individually
        chrom_df = df[df["Chrom"] == chrom]
        max_position = chrom_df["Position"].max()
        num_windows = (max_position // window_size) + 1 #calculate windows
        window_means = [] #empty list to be filled with means

        for window in range(num_windows): #go through each window
            start = window * window_size
            end = start + window_size
            window_data = chrom_df[(chrom_df["Position"] >= start) & (chrom_df["Position"] < end)]
            maf_vector = np.zeros(window_size) #vector of zeros
            
            for _, row in window_data.iterrows(): #fill windows with MAF at relative position
                pos_in_window = int(row["Position"] - start)
                maf_vector[pos_in_window] = row["MAF"]

            mean_maf = np.mean(maf_vector) #calculate means
            window_means.append((start, mean_maf))
        
        window_means_df = pd.DataFrame(window_means, columns=["Window", "Mean_MAF"])
        
        plt.figure(figsize=(12, 10))
        sns.lineplot(data=window_means_df, x="Window", y="Mean_MAF", color=colors[i])
        plt.title(f" Normalized mean MAF per {window_size} bp Window for {sample_name} - Chromosome {chrom}", fontsize=22)
        plt.xlabel("Genomic Position (bp)", fontsize=22)
        plt.ylabel("Mean MAF", fontsize=22)
        xticks = window_means_df["Window"][::20]
        xlabels = [f"{pos:,}" for pos in xticks]
        plt.xticks(ticks=xticks, labels=xticks, fontsize=16)
        plt.yticks(fontsize=16)
        
        plt.tight_layout()
        plt.show()
