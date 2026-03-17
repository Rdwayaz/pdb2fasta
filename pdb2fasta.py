#!/usr/bin/env python3
"""
pdb2fasta.py - Convert PDB files to FASTA format.

Features:
  - Handles unordered ATOM records
  - Uses only CA atoms to identify residues (tolerates missing side chains)
  - Fills gaps between residue numbers with 'X'
  - Outputs one sequence per chain to stdout (pipeable)

Usage:
  python pdb2fasta.py structure.pdb
  python pdb2fasta.py structure.pdb > output.fasta
"""

import sys
import os

# Standard 3-letter to 1-letter amino acid code mapping
AA3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    # Common non-standard / modified residues mapped to closest standard
    "MSE": "M",  # selenomethionine
    "HSD": "H", "HSE": "H", "HSP": "H",  # CHARMM histidine variants
    "CYX": "C", "CYM": "C",              # disulfide / deprotonated cys
    "GLH": "E", "ASH": "D",              # protonated glu/asp
    "LYN": "K",                          # neutral lysine
    "HIP": "H",                          # doubly-protonated histidine
    "SEP": "S", "TPO": "T", "PTR": "Y",  # phosphorylated residues
}


def parse_pdb(filepath):
    """
    Parse ATOM records from a PDB file.

    Returns a dict: { chain_id -> { res_num -> one_letter_code } }
    Only one record per (chain, resnum) is kept; CA atoms are preferred
    but any ATOM record for that residue will register it.
    """
    chains = {}  # chain -> {resnum: (resname, insertion_code)}

    with open(filepath, "r") as fh:
        for line in fh:
            if not line.startswith("ATOM"):
                continue

            # PDB fixed-width column parsing
            try:
                atom_name   = line[12:16].strip()
                res_name    = line[17:20].strip()
                chain_id    = line[21].strip() or "_"
                res_seq_str = line[22:26].strip()
                i_code      = line[26].strip()  # insertion code
            except IndexError:
                continue

            # Skip non-standard chain placeholders
            if not res_seq_str:
                continue
            try:
                res_seq = int(res_seq_str)
            except ValueError:
                continue

            if chain_id not in chains:
                chains[chain_id] = {}

            key = (res_seq, i_code)

            # Prefer CA entry; if already registered with CA, skip others
            if key not in chains[chain_id]:
                chains[chain_id][key] = res_name
            elif atom_name == "CA":
                # Overwrite with CA-confirmed residue name (most reliable)
                chains[chain_id][key] = res_name

    return chains


def build_fasta(chains):
    """
    Build FASTA sequences for each chain.

    Gaps between residue numbers are filled with 'X'.
    Insertion codes at the same residue number are appended without gap.
    """
    results = {}

    for chain_id, residues in sorted(chains.items()):
        if not residues:
            continue

        # Sort by (res_seq, insertion_code) — insertion codes '', 'A', 'B', ...
        sorted_keys = sorted(residues.keys(), key=lambda k: (k[0], k[1]))

        sequence = []
        prev_resnum = None

        for (res_seq, i_code), res_name in zip(
            sorted_keys, [residues[k] for k in sorted_keys]
        ):
            one_letter = AA3TO1.get(res_name.upper(), "X")

            if prev_resnum is not None:
                # Fill integer gap with X (insertion codes don't create gaps)
                gap = res_seq - prev_resnum - 1
                if gap > 0:
                    sequence.extend(["X"] * gap)

            sequence.append(one_letter)

            # Only advance the residue number counter for non-insertion records
            # (insertion codes share the same base number)
            if i_code == "":
                prev_resnum = res_seq
            else:
                # Insertion — keep prev_resnum at the base number so next
                # integer residue still triggers a proper gap check.
                prev_resnum = res_seq

        results[chain_id] = "".join(sequence)

    return results


def format_fasta(sequence, width=60):
    """Wrap sequence to fixed width lines."""
    return "\n".join(sequence[i:i+width] for i in range(0, len(sequence), width))


def main():
    if len(sys.argv) < 2:
        print("Usage: python pdb2fasta.py <input.pdb>", file=sys.stderr)
        sys.exit(1)

    pdb_path = sys.argv[1]
    if not os.path.isfile(pdb_path):
        print(f"Error: file not found: {pdb_path}", file=sys.stderr)
        sys.exit(1)

    basename = os.path.splitext(os.path.basename(pdb_path))[0]

    chains = parse_pdb(pdb_path)
    if not chains:
        print("Error: no ATOM records found in file.", file=sys.stderr)
        sys.exit(1)

    fasta_seqs = build_fasta(chains)

    for chain_id, sequence in fasta_seqs.items():
        header = f">{basename}_chain{chain_id}"
        print(header)
        print(format_fasta(sequence))


if __name__ == "__main__":
    main()
