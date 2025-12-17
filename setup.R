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
