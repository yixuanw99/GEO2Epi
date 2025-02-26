# GEO2Epi: A Workflow for Epigenetic Methylation Analysis

GEO2Epi is a modular, reproducible workflow for downloading, processing, and analyzing DNA methylation data from the Gene Expression Omnibus (GEO). The pipeline downloads GEO datasets, extracts and preprocesses raw data, normalizes IDAT files using CHAMP, and performs downstream analyses with BioLearn and MethylDetectR. It is implemented using Snakemake and leverages both Python and R scripts.

---

## Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Requirements](#requirements)
- [Installation](#installation)
- [Directory Structure](#directory-structure)
- [Pipeline Workflow](#pipeline-workflow)
- [Usage](#usage)
- [Scripts Description](#scripts-description)
- [Contributing](#contributing)
- [License](#license)
- [Contact](#contact)

---

## Overview

The GEO2Epi pipeline automates the processing of GEO methylation datasets. It downloads GEO data (specified by GSE IDs), extracts the required XML and raw tar files, generates sample sheets, normalizes IDAT files, and cleans up raw data. After preprocessing, it performs epigenetic age prediction using two different analytical modules:
- **BioLearn** – A Python-based tool running multiple methylation clocks.
- **MethylDetectR** – An R-based analysis to compute epigenetic age predictors.

By uniting these modules, GEO2Epi offers a comprehensive solution for epigenetic data analysis.

---

## Features

- **Automated Data Retrieval:** Downloads MINiML XML and supplementary raw tar files from GEO.
- **Preprocessing:** Extracts XML data, generates sample sheets, and processes IDAT files.
- **Normalization:** Uses the CHAMP R package to normalize methylation data.
- **Downstream Analysis:** Runs BioLearn and MethylDetectR to calculate various epigenetic age predictions.
- **Modular Workflow:** Built with Snakemake, allowing scalability and reproducibility.
- **Multi-language Support:** Integrates Python and R scripts for comprehensive data analysis.

---

## Requirements

- **Conda:** Create an isolated environment using the provided `environment.yml`.
- **Snakemake:** For workflow management.
- **R (v4.2.3):** With necessary packages (e.g., CHAMP, data.table).
- **Python (v3.10.9):** With packages such as pandas and biolearn.
- **UNIX Shell:** For running the provided bash scripts.

---

## Installation

1. **Clone the repository:**

   ```bash
   git clone https://github.com/yixuanw99/GEO2Epi.git
   cd GEO2Epi
   ```

2. **Create the Conda environment:**

   ```bash
   conda env create -f environment.yml
   conda activate GEO2Epi
   ```

3. **Install ChAMP from GitHub:**

   ```bash
   git clone https://github.com/YuanTian1991/ChAMPdata.git
   R CMD INSTALL ChAMPdata

   git clone https://github.com/YuanTian1991/ChAMP.git
   R CMD INSTALL ChAMP

   mv /home2/yxwu/miniconda3/envs/GEO2Epi/lib/R/library/ChAMPdata/data/probe.features.epic.rda /home2/yxwu/miniconda3/envs/GEO2Epi/lib/R/library/ChAMPdata/data/probe.features.epicv1.rda
   ```

4. **Configure the pipeline:**

   Edit the configuration file at `config/config.yaml` to set the following parameters:
   - `PROJECT_ROOT`: The absolute path to your project directory.
   - `GSE_ID_LIST`: A list of GEO Series IDs to be processed.
   - `CORES`: Number of CPU cores to use (default is 1).

---

## Directory Structure

The repository is organized as follows:

```
.
├── config/                 # Configuration files
│   └── config.yaml         # Main configuration file
├── data/                   # Data directory
│   ├── raw/                # Raw downloaded data
│   ├── processed/          # Processed normalization outputs
│   └── results/            # Analysis results
│       ├── biolearn/       # Outputs from BioLearn analysis
│       └── MethylDetectR/  # Outputs from MethylDetectR analysis
├── logs/                   # Log files from workflow execution
├── scripts/                # Analysis and helper scripts
│   ├── biolearn_run.py     # Python script for BioLearn analysis
│   ├── champ_idat2norm.R   # R script for CHAMP normalization
│   ├── download_gse.sh     # Shell script to download GEO files
│   ├── extract_miniml.sh   # Shell script to extract MINiML XML
│   ├── extract_raw.sh      # Shell script to extract raw tar files and IDATs
│   ├── MethylDetectR_run.R # R script to run MethylDetectR analysis
│   ├── sample_sheet_generator.py  # Python script to generate sample sheets
│   └── MethylDetectR_dependencies/ # Dependencies required for MethylDetectR
├── environment.yml         # Conda environment configuration file
├── Snakefile               # Snakemake workflow file
└── README.md               # Project documentation (this file)
```

---

## Pipeline Workflow

The pipeline is defined in the `Snakefile` and includes the following rules:

- **create_global_dirs:** Creates necessary directories for processed data and results.
- **download_gse:** Downloads GEO files (MINiML XML and raw tar archives) using `download_gse.sh`.
- **extract_miniml:** Extracts the MINiML XML file from the downloaded archive.
- **extract_raw:** Extracts raw data from the tar archive, moves and decompresses IDAT files.
- **generate_sample_sheet:** Parses the XML file to generate a sample sheet and determine the array type.
- **champ_idat2norm:** Normalizes IDAT data using CHAMP via an R script.
- **clean_raw:** Cleans up raw files after processing.
- **run_biolearn:** Executes the BioLearn analysis module.
- **run_MethylDetectR:** Runs the MethylDetectR analysis module.

Each step logs its progress and errors to individual log files located in the `logs/` directory.

---

## Usage

After installing the dependencies and configuring the project:

1. **Activate the environment:**

   ```bash
   conda activate GEO2Epi
   ```

2. **Run the Snakemake pipeline:**

   ```bash
   snakemake --cores 4 --verbose --keep-going --resources download_jobs=1
   ```

   Replace `4` with the number of cores you wish to allocate. Limiting download jobs with `--resources download_jobs=1` is important to avoid being banned by NCBI.

3. **Monitoring logs:**

   Check the `logs/` folder for details on the execution of each rule.

---

## Scripts Description

- **download_gse.sh:** Downloads GEO dataset files based on the provided GSE ID.
- **extract_miniml.sh:** Extracts the MINiML XML file from the downloaded tgz archive.
- **extract_raw.sh:** Extracts raw tar files, moves IDAT files, and decompresses them.
- **sample_sheet_generator.py:** Parses the XML to create a CSV sample sheet and write the array type.
- **champ_idat2norm.R:** Uses the CHAMP package to normalize the methylation data from IDAT files.
- **biolearn_run.py:** Runs multiple epigenetic prediction models from the BioLearn library.
- **MethylDetectR_run.R:** Performs epigenetic age prediction using the MethylDetectR module.

---

## Contributing

Contributions are welcome! Please follow these steps:
1. Fork the repository.
2. Create a new branch for your feature or bugfix.
3. Make your changes and ensure tests pass.
4. Submit a pull request with a detailed description of your changes.

---

## License

This project is distributed under the MIT License. See the [LICENSE](LICENSE) file for more details.

---

## Contact

For questions, suggestions, or contributions, please contact:

- **Project Lead:** [Yi-Hsuan Wu](mailto:yixuan.w99@gmail.com)
- **GitHub Issues:** Use the repository’s issue tracker for bug reports or feature requests.

---

This README provides an overview of the GEO2Epi pipeline, along with detailed installation, usage, and development instructions. For further documentation or updates, please refer to the repository’s Wiki or documentation folder.
