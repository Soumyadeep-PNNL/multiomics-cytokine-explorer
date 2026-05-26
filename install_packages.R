# ============================================================
# install_packages.R
# Run this script ONCE before launching the Shiny app
# ============================================================

# Install CRAN packages
cran_packages <- c(
  "shiny",
  "bslib",        # Modern Bootstrap themes
  "dplyr",
  "tidyr",
  "stringr",
  "ggplot2",
  "plotly",       # Interactive plots
  "DT",           # Interactive tables
  "pheatmap",     # Heatmap
  "ggrepel",      # Non-overlapping labels on volcano
  "RColorBrewer",
  "viridis",
  "scales",
  "readxl",       # Read .xlsx files
  "grid",
  "UpSetR"
)

installed <- rownames(installed.packages())
to_install <- cran_packages[!cran_packages %in% installed]

if (length(to_install) > 0) {
  message("Installing: ", paste(to_install, collapse = ", "))
  install.packages(to_install, repos = "https://cloud.r-project.org")
} else {
  message("All required packages are already installed!")
}

# Verify
missing <- cran_packages[!cran_packages %in% rownames(installed.packages())]
if (length(missing) > 0) {
  warning("The following packages failed to install: ", paste(missing, collapse=", "))
} else {
  message("\nAll packages installed successfully. You can now run the app with:")
  message("  shiny::runApp('app.R')")
  message("  or open app.R in RStudio and click 'Run App'")
}
