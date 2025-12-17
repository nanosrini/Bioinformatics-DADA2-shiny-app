#### FILE: setup.R
# Run this once to install recommended packages (requires R installed)
install_if_missing <- function(pkgs) {
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  }
}

# CRAN
install_if_missing(c('shiny','DT','zip','shinyFiles'))

# Bioconductor and DADA2
if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')
BiocManager::install(c('dada2','phyloseq','ShortRead'), ask = FALSE)

# Note: Kraken2 is an external binary. Please install Kraken2 separately:
cat("\nPlease install Kraken2 separately and ensure 'kraken2' is in your PATH.\n")
cat("Kraken2: https://github.com/DerrickWood/kraken2\n")

cat("\nSetup complete. Run the app with: shiny::runApp('app.R')\n")

```
#### FILE: README.md
# DADA2 + SILVA and Kraken2 Local Shiny App

# DADA2 Shiny App

A local, user-friendly Shiny application for NGS data analysis using the DADA2 pipeline.
The app supports single-end and paired-end FASTQ files and performs quality filtering,
denoising, ASV inference, chimera removal, and taxonomic assignment using the SILVA database.

## Features
- Local execution (Windows & macOS)
- Supports FASTQ / FASTQ.GZ
- DADA2-based ASV inference
- SILVA taxonomy assignment
- Optional Kraken2 classification
- No data upload to external servers

## Installation
```r
install.packages(c("shiny", "DT", "shinyFiles", "fs"))

## Quick start
1. Install R (>=4.0) from CRAN and optionally RStudio.
2. Clone or download this repo (place app.R, setup.R, README.md together).
3. In R, run `setup.R` to install required packages:
```r
source('setup.R')

```
## Notes and limitations
- This app expects paired files to include `_R1/_R2` or `_1/_2` in filenames for pairing. If not, please rename before upload.
- The app runs single-core (multithread disabled) for simplicity and wide compatibility.
- SILVA and Kraken2 databases are large. Please download separately before running the app.
- For very large FASTQ files, users should run the pipeline in R directly or use Docker (not provided here).

``````
4. Install Kraken2 separately if you plan to use Kraken2. Ensure `kraken2` is in your PATH.
   - Kraken2: https://github.com/DerrickWood/kraken2
5. Download SILVA (if using DADA2) and note the fasta path. Example: `SILVA_138_SSURef_NR99.fa.gz`.
   - SILVA: https://www.arb-silva.de/
6. Run the app:

```r
shiny::runApp('app.R')
```


