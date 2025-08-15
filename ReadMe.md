# Hydrocarbon Property Estimator

## Purpose
This project provides a Python implementation of a group-contribution method to estimate physical properties of hydrocarbons. 
It calculates the normal boiling point (Tb in K) and liquid density at 20°C (d20 in kg/m³) from SMILES strings.
Condensed rings are not supported as per the original paper. 

The script processes a list of SMILES from an input file, computes properties, and outputs results to an Excel file with sheets for properties and group counts.

**Note:** Freezing point estimation are not yet implemented.

## Testcases
SMILES list of hydrocarbons from Table 7–14 in the paper are included in `SMILES_Testcase.txt`.

## Usage
1. Install dependencies (see `requirements.txt`).
2. Place your SMILES input in `SMILES_List.txt`.
3. Run the Python script:
```bash

python Tbd20.py
