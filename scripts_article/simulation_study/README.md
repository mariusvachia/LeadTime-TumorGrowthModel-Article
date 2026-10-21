# Simulation Study

This directory contains all R scripts, random seeds, generated datasets, and fitted model/LeadTime objects required to reproduce the simulation study presented in the manuscript.

---

## Overview

The simulations address three main scientific objectives from the paper:

1. **Objective (i) — Model Selection:** Evaluating model selection performance via AIC.
2. **Objective (ii) — Model Performance:** Assessing model estimation performance under 3 scenarios.
3. **Objective (iii) — Lead-Time Corrected Survival:** Estimating lead-time corrected survival accuracy.

---

## Directory Structure

```text
simulation_study/
├── AIC/                         # Objective (i): Model selection performance
│   ├── n500/                    # Simulations with sample size N = 500
│   └── n2000/                   # Simulations with sample size N = 2000
├── Models_LeadTime/             # Objectives (ii) & (iii): Model performance & lead-time corrected survival
└── README.md                    # Documentation for the simulation pipeline
```

## Scripts and Subfolders

Within each simulation pipeline, scripts and files are organized chronologically:

- **`01_simulation*`**: Launches the simulations and handles data generation.
- **`02_estimation*`**: Fits the tumor growth models to the corresponding simulated datasets.
- **`03_.....*`**: Estimates LeadTime objects or/and aggregates fitted results, and generates paper figures/tables.
- **`seeds/`**: Stored files containing random seeds for exact reproducibility.
- **`datasets/`**: Stored generated simulated datasets stratified by scenario/age at onset distribution.
- **`est/`**: Stored the fitted tumor growth model objects by scenario/age at onset distribution.
- **`cor/`**: Stored LeadTime objects.
---
