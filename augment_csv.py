import argparse
import csv
import sys

from reaction_mapping import generate_atom_mapped_reactions


def main():
    parser = argparse.ArgumentParser(
        description="Augment reaction SMILES in a CSV using atom-mapped enumeration."
    )
    parser.add_argument("input_csv", help="Input CSV path.")
    parser.add_argument("output_csv", help="Output CSV path.")
    parser.add_argument(
        "--origin-idx-col",
        default="origin_idx",
        help="Column name for origin_idx (default: origin_idx).",
    )
    parser.add_argument(
        "--origin-rsmi-col",
        default="origin_rsmi",
        help="Column name for origin_rsmi (default: origin_rsmi).",
    )
    parser.add_argument(
        "--aug-idx-col",
        default="aug_idx",
        help="Column name for augmented index (default: aug_idx).",
    )
    parser.add_argument(
        "--keep-empty",
        action="store_true",
        help="Keep rows with no mappings (aug_idx=0).",
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Print per-row progress to stderr.",
    )
    args = parser.parse_args()

    with open(args.input_csv, "r", encoding="utf-8", newline="") as infile:
        reader = csv.DictReader(infile)
        if reader.fieldnames is None:
            raise ValueError("Input CSV has no header.")
        if args.origin_idx_col not in reader.fieldnames:
            raise ValueError(f"Missing column: {args.origin_idx_col}")
        if args.origin_rsmi_col not in reader.fieldnames:
            raise ValueError(f"Missing column: {args.origin_rsmi_col}")

        fieldnames = list(reader.fieldnames)
        if args.aug_idx_col not in fieldnames:
            fieldnames.append(args.aug_idx_col)

        with open(args.output_csv, "w", encoding="utf-8", newline="") as outfile:
            writer = csv.DictWriter(outfile, fieldnames=fieldnames)
            writer.writeheader()

            for row_num, row in enumerate(reader, start=1):
                rsmi = (row.get(args.origin_rsmi_col) or "").strip()
                if not rsmi:
                    if args.keep_empty:
                        out_row = dict(row)
                        out_row[args.aug_idx_col] = "0"
                        writer.writerow(out_row)
                    if args.debug:
                        print(f"[debug] row={row_num} empty rsmi", file=sys.stderr)
                    continue

                mappings = generate_atom_mapped_reactions(
                    rxn_smiles=rsmi,
                )
                if args.debug:
                    origin_idx = row.get(args.origin_idx_col, "")
                    print(
                        f"[debug] row={row_num} origin_idx={origin_idx} mappings={len(mappings)}",
                        file=sys.stderr,
                    )

                if not mappings:
                    if args.keep_empty:
                        out_row = dict(row)
                        out_row[args.origin_rsmi_col] = rsmi
                        out_row[args.aug_idx_col] = "0"
                        writer.writerow(out_row)
                    continue

                for idx, mapped in enumerate(mappings, start=1):
                    out_row = dict(row)
                    out_row[args.origin_rsmi_col] = mapped
                    out_row[args.aug_idx_col] = str(idx)
                    writer.writerow(out_row)


if __name__ == "__main__":
    main()
