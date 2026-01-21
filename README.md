# Enumerate-Reaction-Mapping

This script generates unique atom-mapped reaction SMILES by enumerating atom mappings
and keeping only mappings with the minimum total add+break operations. Only edge
additions and deletions are allowed; atom identities and bond order changes are
treated as add/break based on the bond order delta.

Notes
- Input is parsed with RDKit, then kekulized and explicit hydrogens are added.
- For best results, provide kekulized (non-aromatic) SMILES where possible to avoid
  ambiguity in aromatic bond orders.

## Requirements
- Python
- RDKit
- NetworkX

## Usage

Basic:
```bash
python reaction_mapping.py "O=CCCOO>>CC=O.O=CO"
```

Debug:
```bash
python reaction_mapping.py "O=CCCOO>>CC=O.O=CO" --debug --debug-interval 10000
```

Draw reaction images (SVG):
```bash
python reaction_mapping.py "O=CCCOO>>CC=O.O=CO" --draw --draw-path outputs/rxn.svg --draw-limit 0
```

## Test Cases

Example 1 (with debug output):
```bash
python reaction_mapping.py "O=CCCOO>>CC=O.O=CO" --debug
```
Expected output:
```
[debug] atoms: 12
[debug] symbol_counts: {'O': 3, 'C': 3, 'H': 6}
[debug] group_sizes: {'O': 3, 'C': 3, 'H': 6}
[debug] estimated_mappings: 25920
[debug] done nodes=354 kept=2 drawn=0 best_total=4 best_ops=b2f2 elapsed=0.0s
[O:1]=[C:2]([C:3]([C:4]([O:5][O:6][H:12])([H:10])[H:11])([H:8])[H:9])[H:7]>>[C:3]([C:4](=[O:5])[H:11])([H:8])([H:9])[H:10].[O:1]=[C:2]([O:6][H:12])[H:7]
[O:1]=[C:2]([C:3]([C:4]([O:5][O:6][H:12])([H:10])[H:11])([H:8])[H:9])[H:7]>>[C:4](=[O:5])([O:6][H:12])[H:11].[O:1]=[C:2]([C:3]([H:8])([H:9])[H:10])[H:7]
```

Example 2:
```bash
python reaction_mapping.py "[C:1]([c:4]1[n:2]([N:5]([H:9])[H:10])[c:3]([N:6]([H:11])[H:12])[c:7]([H:13])[n:8]1)([H:14])([H:15])[H:16]>>[C:1]([C:4]1[C:7]([C:3](=[N:2][N:5]([H:9])[H:10])[N:6]([H:11])[H:12])([H:13])[N:8]1)([H:14])([H:15])[H:16]"
```
Expected output:
```
[C:1]([C:2]1=[N:8][C:7]([H:16])=[C:5]([N:6]([H:14])[H:15])[N:3]1[N:4]([H:12])[H:13])([H:9])([H:10])[H:11]>>[C:1]([C:2]1[C:7]([C:5](=[N:3][N:4]([H:12])[H:13])[N:6]([H:14])[H:15])([H:16])[N:8]1)([H:9])([H:10])[H:11]
```

## Output
- Prints unique atom-mapped reaction SMILES to stdout.
- When `--draw` is set, saves SVG files to the path given by `--draw-path`.
  If multiple mappings are saved, numbered suffixes are appended.
