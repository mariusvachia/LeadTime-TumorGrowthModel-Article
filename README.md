# Research Code & R Package for Statistics in Medicine

**Paper Title:** A unified tumor growth model for evaluating breast cancer screening and correcting lead-time bias

This repository contains the R package (`LeadTime`) and all simulation and application scripts for the paper submitted to *Statistics in Medicine*.

---

## Repository Structure

```text
.
├── LeadTime/                        # Source code for the R package developed in the study
└── scripts_article/                 # Scripts and data used for the manuscript
    ├── simulation_study/            # Simulation code, random seeds, datasets, and fitted models/LeadTime
    └── application/                 # Script, example dataset, and R objects obtained on the real registry dataset illustrating the method
```
---

## How to Install the R Package

You can install the package directly from this GitHub repository using `devtools`:

```R
# Install devtools if not already installed
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# Install the package from GitHub
devtools::install_github("mariusvachia/LeadTime-TumorGrowthModel-article", subdir = "LeadTime")
```

