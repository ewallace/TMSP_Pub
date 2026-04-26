# TMSP_Pub

Analysis of protein secretion routes in fungi, focusing on the hydrophobic helix within N-terminal signal peptides and transmembrane regions.

This repository accompanies the paper:

> Protein secretion routes in fungi are mostly determined by the length of the hydrophobic helix in the signal peptide.
>
> Tristan Sones-Dykes, Edward Wallace.
>
> bioRxiv preprint
> 
> DOI: https://doi.org/10.1101/2025.07.30.667231

The analysis here was performed by:

- Tristan Sones-Dykes, [@TristanSones-Dykes](https://github.com/TristanSones-Dykes)
- Edward Wallace, [@ewallace](https://github.com/ewallace)

## Overview

The core hypothesis is that the **length of the hydrophobic helix** in a signal peptide determines whether a protein is routed through the SRP-dependent or Sec63-dependent secretion pathway. The analysis combines:

- **Phobius** predictions of signal peptides and transmembrane regions across 12+ fungal species and humans
- **DeepTMHMM** as an independent prediction method for comparison
- **Hydrophobicity scoring** using multiple scales (Kyte-Doolittle, Hessa, Rose) to compute compound hydropathy metrics
- **Experimental validation** against known SRP-dependent and Sec63-dependent proteins in *S. cerevisiae*
- **GO enrichment analysis** per species, separated by predicted SP and TM proteins
- **PSIPRED** secondary structure predictions as further validation

Key thresholds identified: ~13.5 AA helix length and compound hydropathy score of ~40 separate SRP-dependent from Sec63-dependent proteins.

## Installation

### Prerequisites

- **R** >= 4.5.2
- **Python** >= 3.12
- [**uv**](https://docs.astral.sh/uv/) (Python package manager, recommended) or pip
- **renv** (R package, installed automatically on first R session)

### R Environment

R dependencies are managed with [renv](https://rstudio.github.io/renv/). On first opening the project in R:

```r
# renv should bootstrap automatically via .Rprofile
# If not, install it manually:
install.packages("renv")

# Restore the project library from the lockfile:
renv::restore()
```

This installs all R packages at the exact versions recorded in `renv.lock`, including Bioconductor packages (Biostrings, ggtree, treeio, etc.).

### Python Environment

Python dependencies are used via `reticulate` for Phobius web scraping, DeepTMHMM parsing, and PSIPRED analysis. Set up with [uv](https://docs.astral.sh/uv/):

```bash
# Create and activate virtual environment
uv venv
source .venv/bin/activate

# Install required packages
uv pip install mechanicalsoup requests pandas numpy biopython beautifulsoup4 lxml
```

Or with pip:

```bash
python -m venv .venv
source .venv/bin/activate
pip install mechanicalsoup requests pandas numpy biopython beautifulsoup4 lxml
```

### External Tools

- **Phobius**: No local installation needed. The analysis uses the [Phobius web server](https://phobius.sbc.su.se/) via `mechanicalsoup`.
- **DeepTMHMM**: Results must be obtained separately from [DTU Health Tech](https://dtu.biolib.com/DeepTMHMM) and placed in `results/deepTMHMM/`.
- **PSIPRED**: Uses the [UCL PSIPRED API](http://bioinf.cs.ucl.ac.uk/psipred/) - no local installation needed. Set `PSIPRED_EMAIL` before submitting jobs.
- **FungiDB**: GO enrichment uses the [FungiDB](https://fungidb.org/) web API via `src/GO_analysis.sh`.

## Repository Structure

```
TMSP_Pub/
├── README.md                 # This file
├── LICENSE                   # MIT License
├── .Rprofile                 # Activates renv on R startup
├── renv.lock                 # R package versions (reproducibility)
├── .lintr                    # R linter configuration
│
├── data/                     # Input data
│   ├── Proteins/
│   │   ├── pub/              # Proteomes truncated to first 60 AA
│   │   └── full/             # Full-length proteomes (for DeepTMHMM)
│   ├── SC_SRP.txt            # Verified SRP-dependent proteins (S. cerevisiae)
│   ├── SC_non_SRP.txt        # Verified Sec63-dependent proteins
│   ├── SC_screened.txt       # Screened non-SRP proteins
│   └── scales.csv            # Hydrophobicity scales (KD, Rose, Hessa, etc.)
│
├── src/                      # Source functions and pipelines
│   ├── hydrophobicity.r      # Hydrophobicity calculations & Phobius R wrapper
│   ├── utils.r               # Visualization helpers & cladogram plotting
│   ├── deps.r                # IDE/development dependencies
│   ├── phobius.py            # Phobius web scraper (Python)
│   ├── deepTMHMM.py          # DeepTMHMM output parser (Python)
│   ├── GO_analysis.sh        # FungiDB GO enrichment shell script
│   └── psipred/              # PSIPRED secondary structure prediction
│       └── PSIPRED.py
│
├── scripts/                  # Analysis and figure generation
│   ├── pub_data.r            # Data processing pipeline (run first)
│   └── pub_figures.r         # Publication figure generation (run second)
│
└── results/                  # Analysis outputs
    ├── figures/              # Generated publication figures (PDF + PNG)
    ├── phobius/              # Phobius prediction CSVs per species
    ├── proteins/             # Protein ID lists for GO analysis
    ├── GO/                   # GO enrichment results per species
    └── deepTMHMM/            # DeepTMHMM predictions (3 species)
```

## Usage

### Running the Analysis

The analysis is split into two scripts that should be run in order:

**1. Generate processed data** (`scripts/pub_data.r`):

```r
source("scripts/pub_data.r")
```

This script:
- Runs Phobius on all species proteomes (results are cached in `results/phobius/`)
- Calculates hydrophobicity scores for *S. cerevisiae* using KD, Rose, and Hessa scales
- Writes `results/figures/SC_first_60.csv` with compound hydropathy metrics
- Generates SP/TM protein lists for GO analysis
- Runs GO enrichment via FungiDB (skips if results already exist)
- Processes additional/industry-relevant species

**2. Generate publication figures** (`scripts/pub_figures.r`):

```r
source("scripts/pub_figures.r")
```

This script uses the processed data to produce all manuscript figures:
- Scatter plots with marginal histograms (KD, Hessa, Rose scales)
- Contingency table mosaic plots with chi-squared tests
- Phylogenetic cladogram with per-species helix length distributions
- Amino acid composition parallel coordinates plot
- Phobius vs DeepTMHMM comparison scatter plot
- PSIPRED secondary structure analysis

Figures are saved as PDF and PNG (300 DPI) in `results/figures/`.

### Key Data Files

| File | Description |
|------|-------------|
| `data/scales.csv` | Hydrophobicity scales (20 AA x 8 scales) |
| `data/SC_SRP.txt` | Experimentally verified SRP-dependent proteins |
| `data/SC_non_SRP.txt` | Experimentally verified Sec63-dependent proteins |
| `results/figures/SC_first_60.csv` | Processed *S. cerevisiae* data with hydropathy scores |
| `data/Proteins/pub/proteome_table_*.txt` | Species metadata tables (filename, taxonomy ID, etc.) |

### Notes

- **Caching**: Phobius results are cached in `results/phobius/`. Delete the CSV to re-run predictions for a species.
- **Phobius web server**: The analysis depends on the Phobius web server being available. Large proteomes may require patience or batching.
- **DeepTMHMM**: Results are not generated by this pipeline. Run DeepTMHMM separately and place output in `results/deepTMHMM/<species>/`.
- **Experimental data**: Validation protein lists sourced from [Ast et al., Cell 2013](https://doi.org/10.1016/j.cell.2013.02.003).
