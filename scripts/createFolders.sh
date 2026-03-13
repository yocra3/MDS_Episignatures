#!/bin/bash
# -----------------------------------------------------------------------------
# Script Name:    createFolders.sh
# Description:    Creates the directory structure for the MDS_Episignatures project.
# Author:         Carlos Ruiz
# -----------------------------------------------------------------------------
# Project Directory Structure:
# .
# ├── data/              # Raw and processed data files
# ├── results/           # Analysis results and output files
# ├── scripts/           # Project scripts
# ├── figures/           # Figures
# └── docker/            # Docker files
# -----------------------------------------------------------------------------

# Create main project directories
mkdir -p data
mkdir -p results
mkdir -p scripts
mkdir -p figures
mkdir -p docker