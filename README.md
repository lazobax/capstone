# A Bioinformatics Overview of Publicly Available Fanconi Anemia Mutation Data

**Author:** Lazo B. Ali
**Supervisor:** Dr. Joseph Rebehmed
**Institution:** Lebanese American University, Department of Computer Science and Mathematics

## 📘 Overview

This project systematically integrates and analyzes mutation data related to **Fanconi Anemia (FA)** from two major public resources:

* **LOVD3 (Leiden Open Variation Database)**
* **NCBI ClinVar**

Despite nominal adherence to ACMG/AMP guidelines by both databases, inconsistencies in transcript versions, annotation formats, and classification schemes complicate cross-database meta-analysis. This pipeline resolves those issues through data normalization, conversion to VCF, annotation with ANNOVAR, and downstream analysis.

The final result is a unified dataset of 91,720 GRCh37-mapped variants with harmonized pathogenicity assertions and functional scores.

---

## 📂 Repository Structure

```
├── scripts/                # All processing and analysis scripts
│   ├── get.py              # Scraper for LOVD variant tables
│   ├── fixScrape.py        # Fix inconsistent LOVD headers
│   ├── parse.py            # Extract and normalize LOVD coordinates
│   ├── hgvsToVCF.py        # Convert HGVS to VCF using Biocommons
│   ├── tsvToVcf.py         # Convert TSVs to basic VCF
│   ├── removeIUPACambg.py  # Filter malformed or ambiguous variants
│   ├── replaceCHROM.py     # Replace CHROM with GRCh37 contig names
│   ├── normVCF.sh          # bcftools-based VCF normalization
│   ├── label.py            # Map consequence strings to pathogenic/benign/VUS
│   ├── plot.py             # Visualization of mutation distributions
|   ├── capstone.ipynb          # Final analysis and figures (Jupyter)
├── reqs.yaml               # environment requirements
├── reproduce.sh            # One-line reproducibility script
├── README.md               # You are here
```

---

## 🔧 How to Run

### 1. Install Dependencies

Create a conda environment with:

```bash
conda create -f reqs.yaml
conda activate reproduceCapstone
```

Make sure `bcftools` and `samtools` are in your PATH.

### 2. Run the Pipeline

Run the entire workflow with:

```bash
bash reproduce.sh
```

This will:

* Download and index the GRCh37 reference genome
* Scrape and parse LOVD variant data
* Format ClinVar and LOVD variants into unified VCFs
* Normalize and annotate them
* Merge, label, and produce analysis-ready data

### 3. Run the Analysis

After running the pipeline, open the Jupyter notebook:

```bash
jupyter notebook capstone.ipynb
```

This will generate all figures and tables used in the final report.

---

## 📈 Key Features

* **Fully automated data retrieval and formatting**
* **HGVS to VCF normalization with `biocommons.hgvs`**
* **Variant consequence harmonization (manual label mapping)**
* **ANNOVAR annotation with `refGene` + `dbNSFP v4.2a`**
* **Visualization of mutation frequency, consequence, and spatial distribution**
* **ECDF analysis of mutation reporting frequency**

---

## 📄 Citation

If you use this pipeline or dataset, please cite:

> Ali, L.B. (2025). *A Bioinformatics Overview of Publicly Available Fanconi Anemia Mutation Data*. Lebanese American University Capstone Project.

---

## 📬 Contact

For questions or collaborations, please contact:
📧 [lazo.bakhtiar@lau.edu.lb](mailto:lazo.ali@lau.edu)

