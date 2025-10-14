DESeq2 Pipeline (Documented R Code)

# 1. Install &amp; Load DESeq2 (install only if not already installed)
if (!requireNamespace(&quot;BiocManager&quot;, quietly = TRUE))
install.packages(&quot;BiocManager&quot;)

Page 35
BiocManager::install(&quot;DESeq2&quot;, ask = FALSE) # Skip this line if already installed

library(DESeq2)

# 2. Load circRNA count data (row names = circRNA IDs, columns = sample counts)
counts_df &lt;- read.csv(&quot;circRNA_matrix.csv&quot;, row.names = 1, check.names = FALSE)

# Check the first few rows
head(counts_df)

# 3. Convert data frame to matrix format (required for DESeq2)
count_matrix &lt;- as.matrix(counts_df)

# 4. Create sample metadata table
# First 7 columns are healthy controls (DRR228774–DRR228780)
# Next 10 columns are FECD patient samples (DRR415358–DRR415367)
sample_names &lt;- colnames(count_matrix)

sample_conditions &lt;- data.frame(
row.names = sample_names,
condition = factor(
c(rep(&quot;Control&quot;, 7), rep(&quot;Patient&quot;, 10)),
levels = c(&quot;Control&quot;, &quot;Patient&quot;)
)
)

# Check that the sample conditions are correct


print(sample_conditions)

# 5. Create DESeq2 dataset object
dds &lt;- DESeqDataSetFromMatrix(
countData = count_matrix,
colData = sample_conditions,
design = ~ condition # Specifies model design formula
)
# 6. Run the full DESeq2 pipeline
dds &lt;- DESeq(dds)

# 7. Extract differential expression results (Patient vs Control)
res &lt;- results(dds, contrast = c(&quot;condition&quot;, &quot;Patient&quot;, &quot;Control&quot;))

# 8. Order results by adjusted p-value (FDR correction)
resOrdered &lt;- res[order(res$padj),]

# 9. Print summary of the differential expression analysis
summary(resOrdered)

# 10. Save complete results to CSV for downstream analysis
write.csv(
as.data.frame(resOrdered),
file = &quot;DESeq2_circRNA_.csv&quot;,
row.names = TRUE
)