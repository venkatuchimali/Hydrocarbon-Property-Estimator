from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem
import pandas as pd
import math

# ========= 1. Base groups =========
# SMARTS patterns for basic functional groups in hydrocarbons.
# These are used to identify core structural fragments.
group_smarts = {
    "CH3": "[CH3;!$(C=*);!$(C#*)]",  # Methyl group, not part of double/triple bond
    "CH2": "[CH2;!$(C=*);!$(C#*)]",  # Methylene group
    "CH": "[CH;!$(C=*);!$(C#*)]",    # Methine group
    "C_q": "[C;H0;!$(C=*);!$(C#*)]", # Quaternary carbon
    "=CH2": "[CH2;$(C=*)]",          # Vinyl methylene
    "=CH-": "[CH;$(C=*)]",           # Vinyl methine
    "=C<": "[C;H0;$(C=*)]",          # Vinyl quaternary
    "=C=": "[C;H0;$(C=*)$(C=*)]",    # Allene carbon
    "≡CH": "[CH;$(C#*)]",            # Acetylenic methine
    "≡C-": "[C;H0;$(C#*)]"           # Acetylenic quaternary
}

# Labels for CH3 proximity correction groups.
ch3_proximity_labels = [
    "C(CH3)3", "C(CH3)2", "C(CH3)",
    "C(CH3)C(CH3)", "C(CH3)C(CH3)2", "C(CH3)C(CH3)3",
    "C(CH3)2C(CH3)2", "C(CH3)2C(CH3)3"
]

# Labels for ring size corrections (3 to 20 membered rings).
ring_size_labels = [f"{n} membered ring" for n in range(3, 21)]

# Labels for aromatic substitution corrections.
aromatic_labels = [
    "Aromatic ring",
    "Ortho-substitution", "Meta-substitution", "Para-substitution",
    "Substitution 1-2-3", "Substitution 1-2-4", "Substitution 1-3-5",
    "Substitution 1-2-5", "Substitution 1-2-6", "Substitution 1-3-4",
    "Substitution 1-2-4-5", "Substitution 1-2-3-4", "Substitution 1-2-3-5",
    "Substitution 1-2-3-4-5"
]

# ========= 2. Δθ values =========
# Group contribution values (Δθ) for Tb (boiling point) and d20 (density).
# Sourced from Table 5 in the paper.
delta_theta = {
    # Main groups
    "CH3": (33.6103, 0.653292), "CH2": (6.01945, 0.0467606),
    "CH": (-22.8076, -0.57689), "C_q": (-51.0314, -1.20711),
    "=CH2": (32.9791, 0.652976), "=CH-": (6.33765, 0.0346579),
    "=C<": (-21.1492, -0.589958), "=C=": (7.94780, 0.0210506),
    "≡CH": (32.7471, 0.63950), "≡C-": (7.19358, 0.0228269),
    # CH3 proximity
    "C(CH3)3": (-1.73388, 0.0401010), "C(CH3)2": (-0.91393, 0.0199574),
    "C(CH3)": (-0.55098, 0.00470667), "C(CH3)C(CH3)": (0.71018, 0.00208435),
    "C(CH3)C(CH3)2": (1.21475, 0.00014474), "C(CH3)C(CH3)3": (1.82757, -0.00522518),
    "C(CH3)2C(CH3)2": (2.12237, -0.00461323), "C(CH3)2C(CH3)3": (3.45589, -0.0103234),
    # cis/trans
    "cis": (-0.14240, 0.00237399), "trans": (-0.20120, 0.00404111),
    # Ring sizes
    **{f"{n} membered ring": (val, dens) for n, val, dens in [
        (3, 56.7167, 1.21925), (4, 56.7215, 1.19909),
        (5, 56.7392, 1.20410), (6, 57.9497, 1.18562),
        (7, 59.9484, 1.14001), (8, 61.5173, 1.13091),
        (9, 62.4661, 1.12504), (10, 62.9563, 1.12153),
        (11, 62.5620, 1.11894), (12, 62.1962, 1.11816),
        (13, 61.8293, 1.11867), (14, 61.0686, 1.11929),
        (15, 60.5805, 1.12084), (16, 59.9864, 1.12262),
        (17, 59.2470, 1.12462), (18, 58.3214, 1.12685),
        (19, 57.6199, 1.12931), (20, 56.6837, 1.13200)]},
    # Aromatic corrections
    "Aromatic ring": (56.1560, 1.20018),
    "Ortho-substitution": (0.28767, 0.00595223),
    "Meta-substitution": (0.07109, 0.0115111),
    "Para-substitution": (0.85273, 0.0115111),
    "Substitution 1-2-3": (1.29240, 0.0124819), "Substitution 1-2-4": (0.97349, 0.0210110),
    "Substitution 1-3-5": (-0.67390, 0.0243517), "Substitution 1-2-5": (1.30243, 0.0133045),
    "Substitution 1-2-6": (1.14308, 0.0113885), "Substitution 1-3-4": (-0.35368, 0.0195087),
    "Substitution 1-2-4-5": (2.40917, 0.0313636), "Substitution 1-2-3-4": (4.09709, 0.0231769),
    "Substitution 1-2-3-5": (2.32485, 0.0299362), "Substitution 1-2-3-4-5": (7.39249, 0.0321448)
}

# ========= 3. Aromatic numbering helpers =========
# Helper function to renumber atoms in a ring starting from a specific atom.
# This ensures consistent positioning for substitution analysis.
def renumber_ring_atoms(mol, ring_atoms, start_atom):
    order = [start_atom]
    visited = {start_atom}
    current, prev = start_atom, None
    while len(order) < len(ring_atoms):
        nbrs = [n.GetIdx() for n in mol.GetAtomWithIdx(current).GetNeighbors()
                if n.GetIdx() in ring_atoms and n.GetIdx() != prev]
        for nbr in nbrs:
            if nbr not in visited:
                order.append(nbr)
                visited.add(nbr)
                prev, current = current, nbr
                break
    return order

# Analyzes substitution patterns in aromatic rings.
# Identifies positions and counts specific multi-substitutions.
# If debug=True, prints intermediate info for troubleshooting.
def analyze_aromatic_ring(mol, ring_atoms, debug=False):
    size = len(ring_atoms)
    subs_atoms = []
    for ri in ring_atoms:
        for nbr in mol.GetAtomWithIdx(ri).GetNeighbors():
            if nbr.GetIdx() not in ring_atoms and nbr.GetSymbol() != "H":
                if nbr.GetSymbol() == "C" and nbr.GetTotalNumHs() == 2:
                    if any(nn.GetIsAromatic() and nn.GetIdx() not in ring_atoms for nn in nbr.GetNeighbors()):
                        subs_atoms.append(ri)
                        break
                else:
                    subs_atoms.append(ri)
                    break
    subs_atoms = sorted(set(subs_atoms))
    if not subs_atoms:
        return {lbl:0 for lbl in aromatic_labels[1:]}, []
    bridge_atoms = []
    for ri in ring_atoms:
        for nbr in mol.GetAtomWithIdx(ri).GetNeighbors():
            if nbr.GetIdx() not in ring_atoms:
                if nbr.GetIsAromatic():
                    bridge_atoms.append(ri)
                else:
                    for nn in nbr.GetNeighbors():
                        if nn.GetIsAromatic() and nn.GetIdx() not in ring_atoms:
                            bridge_atoms.append(ri)
    bridge_atoms = sorted(set(bridge_atoms))
    start_atom = min(bridge_atoms) if bridge_atoms else min(subs_atoms)
    ring_order = renumber_ring_atoms(mol, ring_atoms, start_atom)
    subs_pos = [ring_order.index(a) + 1 for a in subs_atoms]
    if debug:
        print("Ring atoms:", ring_atoms)
        print("Renumbered order:", ring_order)
        print("Sub positions:", subs_pos)
    counts = {lbl:0 for lbl in aromatic_labels[1:]}
    if size == 6:
        for i in range(len(subs_pos)):
            for j in range(i+1, len(subs_pos)):
                d = (subs_pos[j] - subs_pos[i]) % size
                dist = min(d, size - d)
                if dist == 1: counts["Ortho-substitution"] += 1
                elif dist == 2: counts["Meta-substitution"] += 1
                elif dist == 3: counts["Para-substitution"] += 1
        sp = set(subs_pos)
        if {1, 2, 3}.issubset(sp): counts["Substitution 1-2-3"] += 1
        if {1, 2, 4}.issubset(sp): counts["Substitution 1-2-4"] += 1
        if {1, 3, 5}.issubset(sp): counts["Substitution 1-3-5"] += 1
        if {1, 2, 5}.issubset(sp): counts["Substitution 1-2-5"] += 1
        if {1, 2, 6}.issubset(sp): counts["Substitution 1-2-6"] += 1
        if {1, 3, 4}.issubset(sp): counts["Substitution 1-3-4"] += 1
        if {1, 2, 4, 5}.issubset(sp): counts["Substitution 1-2-4-5"] += 1
        if {1, 2, 3, 4}.issubset(sp): counts["Substitution 1-2-3-4"] += 1
        if {1, 2, 3, 5}.issubset(sp): counts["Substitution 1-2-3-5"] += 1
        if {1, 2, 3, 4, 5}.issubset(sp): counts["Substitution 1-2-3-4-5"] += 1
    return counts, subs_pos

# ========= 4. Group counting =========
# Counts occurrences of all groups, rings, and aromatic corrections in the molecule.
def count_groups(smiles, debug_aromatics=False):
    mol = Chem.MolFromSmiles(smiles)
    counts = {g:0 for g in group_smarts}
    assigned_atoms = set()
    all_ring_atoms = set()
    for ring in mol.GetRingInfo().AtomRings():
        all_ring_atoms.update(ring)
    ring_counts = {label:0 for label in ring_size_labels}
    for ring in mol.GetRingInfo().AtomRings():
        if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            continue
        size = len(ring)
        if 3 <= size <= 20:
            ring_counts[f"{size} membered ring"] += 1
    aromatic_counts = {label:0 for label in aromatic_labels}
    for ring in mol.GetRingInfo().AtomRings():
        if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            aromatic_counts["Aromatic ring"] += 1
            ring_corrs, _ = analyze_aromatic_ring(mol, ring, debug=debug_aromatics)
            for k, v in ring_corrs.items():
                aromatic_counts[k] += v
    Chem.Kekulize(mol, clearAromaticFlags=True)
    for gname, smarts in group_smarts.items():
        patt = Chem.MolFromSmarts(smarts)
        for match in mol.GetSubstructMatches(patt):
            if match[0] not in assigned_atoms:
                counts[gname] += 1
                assigned_atoms.add(match[0])
    is_olefin = counts['=CH2'] + counts['=CH-'] + counts['=C<'] + counts['=C='] > 0
    # CH3 proximity
    ch3_counts = {label:0 for label in ch3_proximity_labels}
    atom_ch3_count = {}
    methyl_carbons = set()
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == "C" and atom.GetTotalNumHs() == 3:
            methyl_carbons.add(atom.GetIdx())
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != "C": continue
        aid = atom.GetIdx()
        if aid in all_ring_atoms: continue
        if any(n.GetIdx() in all_ring_atoms for n in atom.GetNeighbors()): continue
        n_ch3 = sum(1 for nbr in atom.GetNeighbors() if nbr.GetIdx() in methyl_carbons)
        if n_ch3 > 0:
            atom_ch3_count[aid] = n_ch3
            if n_ch3 == 1:
                found = False
                for nbr in atom.GetNeighbors():
                    if nbr.GetSymbol() == "C" and nbr.GetIdx() not in methyl_carbons and not any(n.GetIdx() in all_ring_atoms for n in nbr.GetNeighbors()):
                        n_nbr_ch3 = sum(1 for nn in nbr.GetNeighbors() if nn.GetIdx() in methyl_carbons)
                        is_special_nbr = (nbr.GetTotalNumHs() == 2 and any(b.GetBondType() == Chem.BondType.DOUBLE for b in nbr.GetBonds()))
                        if n_nbr_ch3 >= 1 or is_special_nbr:
                            ch3_counts["C(CH3)"] += 1
                            found = True
                            break
                if not found:
                    if atom.GetTotalNumHs() == 0 and any(b.GetBondType() == Chem.BondType.DOUBLE for b in atom.GetBonds()):  # for =C< with n_ch3=1
                        ch3_counts["C(CH3)"] += 1
            elif n_ch3 == 2:
                ch3_counts["C(CH3)2"] += 1
            elif n_ch3 == 3:
                ch3_counts["C(CH3)3"] += 1
    # Additional C(CH3) for terminal =CH2 in olefins
    if counts["=CH2"] > 0:
        for atom in mol.GetAtoms():
            if atom.GetSymbol() == "C" and atom.GetTotalNumHs() == 2 and any(b.GetBondType() == Chem.BondType.DOUBLE for b in atom.GetBonds()):
                ch3_counts["C(CH3)"] += 1
    for bond in mol.GetBonds():
        a1, a2 = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if a1 in atom_ch3_count and a2 in atom_ch3_count:
            n1, n2 = atom_ch3_count[a1], atom_ch3_count[a2]
            key = (min(n1, n2), max(n1, n2))
            key_map = {
                (1,1): "C(CH3)C(CH3)",
                (1,2): "C(CH3)C(CH3)2",
                (1,3): "C(CH3)C(CH3)3",
                (2,2): "C(CH3)2C(CH3)2",
                (2,3): "C(CH3)2C(CH3)3"
            }
            if key in key_map:
                ch3_counts[key_map[key]] += 1
    counts.update(ch3_counts)
    counts.update(ring_counts)
    counts.update(aromatic_counts)
    return counts

# ========= 5. Cis/trans counting =========
# Counts cis and trans isomers in double bonds and small rings.
def count_cis_trans(smiles):
    mol = Chem.MolFromSmiles(smiles)
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    cis = trans = 0
    # For double bonds
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.DOUBLE and not bond.IsInRing():
            stereo = bond.GetStereo()
            if stereo == Chem.BondStereo.STEREOZ:
                cis += 1
            elif stereo == Chem.BondStereo.STEREOE:
                trans += 1
    # For non-aromatic rings
    mol2 = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol2, randomSeed=42)
    conf = mol2.GetConformer()
    for ring in mol.GetRingInfo().AtomRings():
        if any(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            continue  # Skip aromatic rings
        if 3 <= len(ring) <= 8:
            subs = []
            for ri in ring:
                atom = mol.GetAtomWithIdx(ri)
                for nbr in atom.GetNeighbors():
                    if nbr.GetIdx() not in ring and nbr.GetSymbol() != "H":
                        subs.append((ri, nbr.GetIdx()))
            if len(subs) == 2:
                (a1, s1), (a2, s2) = subs
                if mol.GetBondBetweenAtoms(a1, a2):  # Substituents must be adjacent
                    dih = AllChem.GetDihedralDeg(conf, s1, a1, a2, s2)
                    abs_dih = abs(dih)
                    if abs_dih < 90:
                        cis += 1
                    else:
                        trans += 1
    return {"cis": cis, "trans": trans}

# ========= 6. Property calculation =========
# Coefficients for equations from Table 5.
# Tb uses Eq. (11): θ = a + b * Σ∆θi + c (Σ∆θi)^m
# d20 uses Eq. (9): M / θ = a + b * Σ∆θi
Tb_coeffs = {"a": 2104.97, "b": -0.17744, "c": -6194.34, "m": -0.28017}
d20_coeffs = {"a": -0.391046, "b": 0.349011}

# Calculates summed contributions (S_tb, S_d), Tb, d20, and group counts.
def calc_properties(smiles, debug_aromatics=False):
    counts = count_groups(smiles, debug_aromatics)
    counts.update(count_cis_trans(smiles))
    S_tb = sum(counts.get(g, 0) * delta_theta.get(g, (0, 0))[0] for g in counts)
    S_d = sum(counts.get(g, 0) * delta_theta.get(g, (0, 0))[1] for g in counts)
    Tb = Tb_coeffs["a"] + Tb_coeffs["b"] * S_tb + Tb_coeffs["c"] * (S_tb ** Tb_coeffs["m"]) if S_tb > 0 else None
    M = Descriptors.MolWt(Chem.MolFromSmiles(smiles))
    denom = d20_coeffs["a"] + d20_coeffs["b"] * S_d
    d20 = M / denom if denom != 0 else None
    return {"S_tb": S_tb, "S_d": S_d, "Tb_K": Tb, "d20_kg_m3": d20, "counts": counts}

if __name__ == "__main__":
    input_file = "smiles_list.txt"  # TXT file with one SMILES per line
    smiles_list = [line.strip() for line in open(input_file) if line.strip()]

    properties_data = []
    group_counts_data = []

    for smi in smiles_list:
        try:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                # Invalid SMILES, append NaN
                properties_data.append({"SMILES": smi, "Tb_K": float('nan'), "d20_kg_m3": float('nan')})
                counts_row = {"SMILES": smi}
                counts_row.update({g: float('nan') for g in list(group_smarts.keys()) + ch3_proximity_labels + ring_size_labels + aromatic_labels + ["cis", "trans"]})
                group_counts_data.append(counts_row)
                continue

            result = calc_properties(smi)
            properties_data.append({
                "SMILES": smi,
                "Tb_K": result["Tb_K"],
                "d20_kg_m3": result["d20_kg_m3"]
            })
            counts_row = {"SMILES": smi}
            counts_row.update(result["counts"])
            group_counts_data.append(counts_row)
        except Exception as e:
            print(f"Error processing {smi}: {e}")
            properties_data.append({"SMILES": smi, "Tb_K": float('nan'), "d20_kg_m3": float('nan')})
            counts_row = {"SMILES": smi}
            counts_row.update({g: float('nan') for g in list(group_smarts.keys()) + ch3_proximity_labels + ring_size_labels + aromatic_labels + ["cis", "trans"]})
            group_counts_data.append(counts_row)

    df_properties = pd.DataFrame(properties_data)
    df_group_counts = pd.DataFrame(group_counts_data)

    output_file = "smiles_properties.xlsx"
    with pd.ExcelWriter(output_file) as writer:
        df_properties.to_excel(writer, sheet_name="Properties", index=False)
        df_group_counts.to_excel(writer, sheet_name="GroupCounts", index=False)

    print(f"Results saved to {output_file}")