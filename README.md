# ESPRE: EccDNA Stacking Prediction for Robust Estimation

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python](https://img.shields.io/badge/Python-3.8%2B-blue)](https://www.python.org/)
[![R](https://img.shields.io/badge/R-4.0%2B-blue)](https://www.r-project.org/)
[![Platform](https://img.shields.io/badge/Platform-Linux%2FUnix-lightgrey)]()


    ███████╗███████╗██████╗ ██████╗ ███████╗
    ██╔════╝██╔════╝██╔══██╗██╔══██╗██╔════╝
    █████╗  ███████╗██████╔╝██████╔╝█████╗  
    ██╔══╝  ╚════██║██╔═══╝ ██╔══██╗██╔══╝  
    ███████╗███████║██║     ██║  ██║███████╗
    ╚══════╝╚══════╝╚═╝     ╚═╝  ╚═╝╚══════╝
    ----------------------------------------
    EccDNA Stacking Prediction for Robust Estimation
ESPRE (EccDNA Stacking Prediction for Robust Estimation) is a bioinformatics tool designed to predict Schizophrenia (SCZ) risk based on extrachromosomal circular DNA (eccDNA) features.

It employs a Stacking Ensemble Learning framework (integrating XGBoost, LightGBM, Random Forest, SVM, and Neural Networks) to analyze genomic bin features—specifically the fragment length ratio and GC content—to provide robust disease probability estimation.
## Features
- Robust Stacking Model: Combines 5 base learners and a meta-learner to reduce overfitting and improve generalization across cohorts.

- Biological Feature Extraction: Automatically calculates GC-corrected fragment length ratios and GC content in 1Mb genomic windows.

- Interpretability: Identifies and utilizes the Top 100 most influential genomic regions for prediction.

- Automated Pipeline: Seamlessly integrates bedtools, R, and Python for end-to-end analysis.
## Installation & Prerequisites
1. System Requirements
- Linux/Unix environment

- Python >= 3.8

- R >= 4.0

- bedtools (Must be in your system PATH)

2. Clone the RepositoryBash

        git clone https://github.com/xuzhe010/ESPRE.git
        cd ESPRE

3. Install Python DependenciesBash
        
        pip install -r requirements.txt

    Key python packages: `pandas`, `numpy`, `scikit-learn`, `xgboost`, `lightgbm`, `joblib`, `seaborn`

4. Install R Packages
Open your R console and run:

        Rinstall.packages(c("tidyverse", "optparse"))

5. Prepare Database (Important!)

You need to generate the 1Mb window file based on your reference genome (hg38).
        # Assuming you have hg38.fa.fai (index file)
        bedtools makewindows -g /path/to/hg38.fa.fai -w 1000000 > database/hg38_1Mb_windows.bed

## Usage
Run the main script espre.py with your eccDNA BED file and the reference genome.

        python espre.py \
            --input /path/to/your_sample.bed \
            --genome /path/to/hg38.fa \
            --output ./result_sample.csv

Arguments

|Argument|Flag|Description|
|:----|:----|:----|
Input|`-i, --input`|Path to the patient eccDNA BED file (Columns: chr start end).
Genome|`-g, --genome`|Path to the hg38 reference genome FASTA file (.fa).
Output|`-o, --output`|Path where the result CSV will be saved.
Temp Dir|--tmp_dir|(Optional) Directory for temporary files. Default: auto-generated.

## Input & Output Format
### Input Format (BED)
Standard 3-column BED file (tab-separated):

    chr1    10000   15000
    chr1    20000   20500


### Output 
Format (CSV)The software generates a CSV file with the prediction result:

|Sample_ID|Prediction|Probability|Note
|----|----|----|----|
|sample01.bed|SCZ|0.8742|Flagged High Risk
