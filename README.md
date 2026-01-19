# Enumerate-Reaction-Mapping

This script generates unique atom-mapped reaction SMILES by enumerating atom mappings
and filtering them with add/break bond counts. Only edge additions and deletions are
allowed; atom identities and bond order changes are treated as add/break based on the
bond order delta.

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

With add/break ranges:
```bash
python reaction_mapping.py "O=CCCOO>>CC=O.O=CO" --add-range 0,3 --del-range 0,3
```

Debug:
```bash
python reaction_mapping.py "O=CCCOO>>CC=O.O=CO" --debug --debug-interval 10000
```

Draw reaction images (SVG):
```bash
python reaction_mapping.py "O=CCCOO>>CC=O.O=CO" --draw --draw-path outputs/rxn.svg --draw-limit 0
```

## Output
- Prints unique atom-mapped reaction SMILES to stdout.
- When `--draw` is set, saves SVG files to the path given by `--draw-path`.
  If multiple mappings are saved, numbered suffixes are appended.

