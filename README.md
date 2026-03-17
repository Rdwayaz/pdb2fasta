pdb2fasta.py - Convert PDB files to FASTA format.

Features:
  - Handles unordered ATOM records
  - Uses only CA atoms to identify residues (tolerates missing side chains)
  - Fills gaps between residue numbers with 'X'
  - Outputs one sequence per chain to stdout (pipeable)

Usage:
  python pdb2fasta.py structure.pdb
  python pdb2fasta.py structure.pdb > output.fasta
