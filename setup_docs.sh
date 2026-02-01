#!/bin/bash

# Make sure we are in the project root
cd "$(dirname "$0")"

# Create main folders
mkdir -p docs/Probability_Statistics
mkdir -p docs/Regression
mkdir -p docs/Survival_Analysis
mkdir -p docs/Clinical_Trials
mkdir -p docs/Meta_Analysis
mkdir -p docs/Infectious_Disease_Modelling
mkdir -p docs/Machine_Learning
mkdir -p docs/Bayesian_Statistics
mkdir -p docs/Reproducible_Research

# Create main markdown files in each folder
echo "# Probability & Statistics Review" > docs/Probability_Statistics/index.md
echo "# Regression" > docs/Regression/index.md
echo "# Survival Analysis" > docs/Survival_Analysis/index.md
echo "# Clinical Trials" > docs/Clinical_Trials/index.md
echo "# Meta Analysis" > docs/Meta_Analysis/index.md
echo "# Infectious Disease Modelling" > docs/Infectious_Disease_Modelling/index.md
echo "# Machine Learning" > docs/Machine_Learning/index.md
echo "# Bayesian Statistics" > docs/Bayesian_Statistics/index.md
echo "# Reproducible Research" > docs/Reproducible_Research/index.md

# Create main index page
echo "# Welcome to BioStatRx" > docs/index.md

echo "Core folders and main pages created successfully!"
