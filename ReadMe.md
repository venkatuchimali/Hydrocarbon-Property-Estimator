# Hydrocarbon Property Estimator – Group Contribution Method

This repository implements the **Skander & Chitour (2002)** group-contribution method for estimating **boiling point (Tb)** and **liquid density at 20 °C (d20)** of pure hydrocarbons using **SMILES** representations of molecules.

It uses **RDKit** for structural parsing, SMARTS pattern matching for functional group detection, and pandas for results export.

The implemented method:
- Is adapted from:  
  **N. Skander & C.E. Chitour (2002). “A New Group-Contribution Method for the Estimation of Physical Properties of Hydrocarbons,” Oil & Gas Science and Technology – Rev. IFP, Vol. 57 (4), 369–376**.
- Distinguishes between **isomers**, ring types, aromatic substitution patterns, and branching effects using **correction groups**.
- Uses equations **(11)** for Tb and **(9)** for density from the paper.

## Features

- **Automatic group counting** via SMARTS patterns:
  - CH₃, CH₂, CH, quaternary C
  - Olefinic, acetylenic carbons
  - Ring size corrections (3–20 member rings)
  - Aromatic ring detection & ortho/meta/para & multi-substitution patterns
  - CH₃ proximity effects for branching corrections
  - cis/trans configurations

- **Calculates**:
  - ΣΔθ contributions for Tb and d20
  - Estimated **boiling point (K)**
  - Estimated **density at 20 °C (kg/m³)**

- **Outputs**:
  - `smiles_properties.xlsx` with:
    - Sheet `"Properties"` — Tb, d20, ΣΔθ values
    - Sheet `"GroupCounts"` — occurrence counts for all groups and corrections

## Installation

1. **Clone this repository:** git clone https://github.com/venkatuchimali/Hydrocarbon-Property-Estimator.git
2. **Install Dependencies:** pip install -r requirements.txt

## Usage

1. **Prepare an input file:** `smiles_list.txt` with one SMILES string per line
2. **Run the script:** python Tbd20.py
3. **View the results:**
- Output is saved in **smiles_properties.xlsx** with:
  - Sheet `Properties`: SMILES, Tb_K, d20_kg_m3
  - Sheet `GroupCounts`: raw group counts used in calculations

## Limitations

- Designed for **pure hydrocarbons only** (n-paraffins, iso-paraffins, olefins, alkynes, naphthenes, aromatics).
- No condensed polycyclic systems in current version (as per original paper's scope).
- Freezing point correlation is **not** implemented in this script (but available in paper).

## Citing

If you use this tool in your work, please cite:

> N. Skander, C.E. Chitour, "A New Group-Contribution Method for the Estimation of Physical Properties of Hydrocarbons", Oil & Gas Science and Technology – Rev. IFP, 57(4), 369–376, 2002.

## License

MIT License – feel free to use, modify, and distribute with attribution.
