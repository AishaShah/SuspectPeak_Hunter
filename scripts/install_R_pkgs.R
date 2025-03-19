# R script to install required packages into a user-defined library

# Set a user-writable library path
user_lib <- "/work/FAC/FBM/CIG/nvastenh/competition_model/Aisha-Dora/R_libs"  # Customize this path if needed
dir.create(user_lib, showWarnings = FALSE, recursive = TRUE)

# Set CRAN repository
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Install BiocManager if not already installed
if (!require("BiocManager", quietly = TRUE, lib.loc = user_lib)) {
    install.packages("BiocManager", lib = user_lib)
}

# List of required CRAN and Bioconductor packages
cran_packages <- c("reshape2", "ggplot2", "tidyr", "dplyr", "patchwork", "ggh4x", "ggVennDiagram", "gridExtra")
bioconductor_packages <- c("DiffBind")

# Function to check if a package is installed
is_package_installed <- function(package) {is.element(package, installed.packages(lib.loc=user_lib)[,1])}

# Install missing CRAN packages
missing_cran <- cran_packages[!sapply(cran_packages, is_package_installed)]
if (length(missing_cran) > 0) {
    install.packages(missing_cran, lib = user_lib)
} else {
    cat("All CRAN packages are already installed.\n")
}

# Install missing Bioconductor packages
missing_bioc <- bioconductor_packages[!sapply(bioconductor_packages, is_package_installed)]
if (length(missing_bioc) > 0) {
    BiocManager::install(missing_bioc, lib = user_lib)
} else {
    cat("All Bioconductor packages are already installed.\n")
}

# Verify installed packages
cat("Installed packages in the library:\n")
installed_packages <- installed.packages(lib.loc = user_lib)
print(installed_packages[, 1])

# Try loading to check
library(dplyr, lib.loc = user_lib)
library(tidyr, lib.loc = user_lib)
library(reshape2, lib.loc = user_lib) # for melt function
library(ggplot2, lib.loc = user_lib)
library(patchwork, lib.loc = user_lib)
library(ggh4x, lib.loc = user_lib)
library(ggVennDiagram, lib.loc = user_lib)
library(gridExtra, lib.loc = user_lib)
library(DiffBind, lib.loc = user_lib)  # Ensure DiffBind is loaded correctly



