# 🧠 GSE125050 AD Gene Explorer

An interactive web application designed to visualize cell-type-specific gene expression patterns in Alzheimer's disease (AD). 

This tool provides researchers with a user-friendly interface to query individual genes (e.g., *APOE*, *BECN1*, *TREM2*) and instantly evaluate their expression levels across diverse brain cell populations (neurons, astrocytes, endothelial cells, and myeloid cells) comparing AD patients with healthy controls.

## 📌 Data Source & System Attribution
* **Raw Dataset Source:** NCBI GEO ([GSE125050](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE125050)) — Cell-type-specific RNA profiles derived from the post-mortem frozen superior frontal gyrus of Alzheimer's disease (AD) patients and non-AD controls.
* **System Development & Data Processing:** Built and processed into this online explorer system by **Wendy Hsu's Lab at National Health Research Institutes (NHRI), Taiwan**.

## 🚀 Features
* **Real-time Query:** Instantly convert gene symbols to Entrez IDs via MyGene API.
* **Statistical Analysis:** Automatic computation of statistical significance (T-test) between AD and Control groups within specific cell types.
* **Publication-Ready Visualization:** High-quality bar plots with individual sample data points for robust data representation.
