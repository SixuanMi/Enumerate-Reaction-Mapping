import argparse
import math
import os
import sys
import time
from collections import Counter

import networkx as nx
from rdkit import Chem
from rdkit.Chem import AllChem, Draw

EPS = 1e-6


def parse_reaction_smiles(rxn_smiles):
    if ">>" in rxn_smiles:
        parts = rxn_smiles.split(">>")
        if len(parts) != 2:
            raise ValueError("Reaction SMILES should contain a single '>>'.")
        reactants, products = parts
    else:
        parts = rxn_smiles.split(">")
        if len(parts) != 3:
            raise ValueError("Reaction SMILES should contain '>>' or three '>' parts.")
        reactants, products = parts[0], parts[2]

    if not reactants or not products:
        raise ValueError("Reaction SMILES must include reactants and products.")

    return reactants, products


def smiles_to_mol(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Could not parse SMILES: {smiles}")
    try:
        Chem.Kekulize(mol, clearAromaticFlags=True)
        mol = Chem.AddHs(mol)
    except Exception as e:
        raise ValueError(f"Could not kekulize molecule from SMILES: {smiles}") from e
    return mol


def mol_to_graph(mol):
    graph = nx.Graph()
    for atom in mol.GetAtoms():
        graph.add_node(atom.GetIdx(), symbol=atom.GetSymbol())
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtomIdx()
        a2 = bond.GetEndAtomIdx()
        graph.add_edge(a1, a2, bond_order=bond.GetBondTypeAsDouble())
    return graph


def atom_symbol_counts(mol):
    return Counter(atom.GetSymbol() for atom in mol.GetAtoms())


def symmetry_ranks(mol):
    return Chem.CanonicalRankAtoms(mol, breakTies=False, includeAtomMaps=False)


def mapping_key(atom_mapping, reactant_ranks, product_ranks):
    counts = {}
    for react_idx, prod_idx in atom_mapping.items():
        r_rank = reactant_ranks[react_idx]
        p_rank = product_ranks[prod_idx]
        counts.setdefault(r_rank, {})
        counts[r_rank][p_rank] = counts[r_rank].get(p_rank, 0) + 1

    key_parts = []
    for r_rank in sorted(set(reactant_ranks)):
        p_counts = counts.get(r_rank, {})
        key_parts.append((r_rank, tuple(sorted(p_counts.items()))))
    return tuple(key_parts)


def bond_order_matrix(mol):
    n_atoms = mol.GetNumAtoms()
    matrix = [[0.0 for _ in range(n_atoms)] for _ in range(n_atoms)]
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtomIdx()
        a2 = bond.GetEndAtomIdx()
        order = bond.GetBondTypeAsDouble()
        matrix[a1][a2] = order
        matrix[a2][a1] = order
    return matrix


def group_atoms_by_symbol(graph):
    groups = {}
    for node, data in graph.nodes(data=True):
        groups.setdefault(data["symbol"], []).append(node)
    return groups


def estimate_mapping_count(groups):
    total = 1
    for nodes in groups.values():
        total *= math.factorial(len(nodes))
    return total


def add_break_bounds_from_delta(min_delta, max_delta):
    if max_delta <= 0:
        min_add = 0.0
        max_add = 0.0
        min_break = -max_delta
        max_break = -min_delta
    elif min_delta >= 0:
        min_add = min_delta
        max_add = max_delta
        min_break = 0.0
        max_break = 0.0
    else:
        min_add = 0.0
        max_add = max_delta
        min_break = 0.0
        max_break = -min_delta
    return min_add, max_add, min_break, max_break


def format_count(value):
    if abs(value - round(value)) < EPS:
        return str(int(round(value)))
    return f"{value:.3f}".rstrip("0").rstrip(".")


def format_best_pairs(pairs):
    if not pairs:
        return "none"
    parts = []
    for break_count, add_count in sorted(pairs):
        parts.append(f"b{format_count(break_count)}f{format_count(add_count)}")
    return ",".join(parts)


def minmax_product_orders(candidates_a, candidates_b, product_orders, allow_same):
    min_p = None
    max_p = None
    for a in candidates_a:
        row = product_orders[a]
        for b in candidates_b:
            if not allow_same and a == b:
                continue
            val = row[b]
            if min_p is None or val < min_p:
                min_p = val
            if max_p is None or val > max_p:
                max_p = val
    return min_p, max_p


def compute_remaining_min_total(
    order,
    pos,
    mapping,
    used_product,
    reactant_orders,
    product_orders,
    reactant_symbols,
    product_indices_by_symbol,
):
    unassigned = order[pos:]
    if not unassigned:
        return 0.0

    remaining_by_symbol = {}
    for symbol, indices in product_indices_by_symbol.items():
        remaining_by_symbol[symbol] = [idx for idx in indices if not used_product[idx]]

    needed = Counter(reactant_symbols[idx] for idx in unassigned)
    for symbol, count in needed.items():
        if len(remaining_by_symbol[symbol]) < count:
            return None

    symbol_pair_minmax = {}
    min_total = 0.0
    assigned = order[:pos]
    minmax_assigned = {}

    for i in unassigned:
        sym_i = reactant_symbols[i]
        candidates_i = remaining_by_symbol[sym_i]
        if not candidates_i:
            return None
        for j in assigned:
            prod_j = mapping[j]
            cache_key = (prod_j, sym_i)
            if cache_key not in minmax_assigned:
                min_p, max_p = minmax_product_orders(
                    [prod_j], candidates_i, product_orders, allow_same=False
                )
                minmax_assigned[cache_key] = (min_p, max_p)
            min_p, max_p = minmax_assigned[cache_key]
            min_delta = min_p - reactant_orders[i][j]
            max_delta = max_p - reactant_orders[i][j]
            a0, a1, b0, b1 = add_break_bounds_from_delta(min_delta, max_delta)
            min_total += a0 + b0

    for idx_a, i in enumerate(unassigned):
        sym_i = reactant_symbols[i]
        for j in unassigned[idx_a + 1 :]:
            sym_j = reactant_symbols[j]
            key = (sym_i, sym_j)
            if key not in symbol_pair_minmax:
                allow_same = sym_i != sym_j
                min_p, max_p = minmax_product_orders(
                    remaining_by_symbol[sym_i],
                    remaining_by_symbol[sym_j],
                    product_orders,
                    allow_same=allow_same,
                )
                if min_p is None:
                    return None
                symbol_pair_minmax[key] = (min_p, max_p)
                symbol_pair_minmax[(sym_j, sym_i)] = (min_p, max_p)
            min_p, max_p = symbol_pair_minmax[key]
            min_delta = min_p - reactant_orders[i][j]
            max_delta = max_p - reactant_orders[i][j]
            a0, a1, b0, b1 = add_break_bounds_from_delta(min_delta, max_delta)
            min_total += a0 + b0

    return min_total


def build_mapped_reaction_smiles(reactant_mol, product_mol, atom_mapping):
    reactant = Chem.Mol(reactant_mol)
    product = Chem.Mol(product_mol)

    for atom in reactant.GetAtoms():
        atom.SetAtomMapNum(atom.GetIdx() + 1)

    for start_idx, target_idx in atom_mapping.items():
        product.GetAtomWithIdx(target_idx).SetAtomMapNum(start_idx + 1)

    reactant_smiles = Chem.MolToSmiles(reactant, canonical=True, allHsExplicit=True)
    product_smiles = Chem.MolToSmiles(product, canonical=True, allHsExplicit=True)
    return f"{reactant_smiles}>>{product_smiles}"


def build_draw_path(base_path, index):
    base, ext = os.path.splitext(base_path)
    if not ext:
        ext = ".svg"
    if index == 0:
        return f"{base}{ext}"
    return f"{base}_{index + 1}{ext}"


def ensure_parent_dir(path):
    parent = os.path.dirname(path)
    if parent:
        os.makedirs(parent, exist_ok=True)


def save_reaction_svg(reactant_mol, product_mol, atom_mapping, output_path):
    reactant = Chem.Mol(reactant_mol)
    product = Chem.Mol(product_mol)

    for atom in reactant.GetAtoms():
        atom.SetAtomMapNum(atom.GetIdx() + 1)

    for start_idx, target_idx in atom_mapping.items():
        product.GetAtomWithIdx(target_idx).SetAtomMapNum(start_idx + 1)

    rxn = AllChem.ChemicalReaction()
    rxn.AddReactantTemplate(reactant)
    rxn.AddProductTemplate(product)
    svg = Draw.ReactionToImage(rxn, useSVG=True)
    ensure_parent_dir(output_path)
    with open(output_path, "w", encoding="utf-8") as handle:
        handle.write(svg)


def generate_atom_mapped_reactions(
    rxn_smiles,
    debug=False,
    debug_interval=100000,
    draw_path=None,
    draw_limit=0,
):
    reactant_smiles, product_smiles = parse_reaction_smiles(rxn_smiles)
    reactant_mol = smiles_to_mol(reactant_smiles)
    product_mol = smiles_to_mol(product_smiles)

    # Check atom number and formula
    if reactant_mol.GetNumAtoms() != product_mol.GetNumAtoms():
        raise ValueError("Reactant and product must have the same number of atoms.")
    if atom_symbol_counts(reactant_mol) != atom_symbol_counts(product_mol):
        raise ValueError("Reactant and product must have the same chemical formula.")

    reactant_graph = mol_to_graph(reactant_mol)
    product_graph = mol_to_graph(product_mol)
    reactant_groups = group_atoms_by_symbol(reactant_graph)
    product_groups = group_atoms_by_symbol(product_graph)
    if reactant_groups.keys() != product_groups.keys():
        return []

    reactant_symbols = [atom.GetSymbol() for atom in reactant_mol.GetAtoms()]
    product_symbols = [atom.GetSymbol() for atom in product_mol.GetAtoms()]
    product_indices_by_symbol = {}
    for idx, sym in enumerate(product_symbols):
        product_indices_by_symbol.setdefault(sym, []).append(idx)

    reactant_orders = bond_order_matrix(reactant_mol)
    product_orders = bond_order_matrix(product_mol)
    reactant_ranks = symmetry_ranks(reactant_mol)
    product_ranks = symmetry_ranks(product_mol)

    if debug:
        counts = atom_symbol_counts(reactant_mol)
        group_sizes = {symbol: len(nodes) for symbol, nodes in reactant_groups.items()}
        est = estimate_mapping_count(reactant_groups)
        print(f"[debug] atoms: {reactant_mol.GetNumAtoms()}", file=sys.stderr)
        print(f"[debug] symbol_counts: {dict(counts)}", file=sys.stderr)
        print(f"[debug] group_sizes: {group_sizes}", file=sys.stderr)
        print(f"[debug] estimated_mappings: {est}", file=sys.stderr)

    order = sorted(
        range(reactant_mol.GetNumAtoms()),
        key=lambda i: (
            len(product_indices_by_symbol[reactant_symbols[i]]),
            -reactant_graph.degree[i],
        ),
    )

    results = set()
    seen_keys = set()
    used_product = [False] * product_mol.GetNumAtoms()
    mapping = {}

    nodes = 0
    kept = 0
    drawn = 0
    best_total = None
    best_pairs = set()
    start_time = time.time()

    def backtrack(pos, add_count, break_count):
        nonlocal nodes, kept, drawn, best_total, best_pairs
        nodes += 1
        if debug and nodes % debug_interval == 0:
            elapsed = time.time() - start_time
            print(
                f"[debug] nodes={nodes} kept={kept} elapsed={elapsed:.1f}s",
                file=sys.stderr,
            )

        current_total = add_count + break_count
        if best_total is not None and current_total > best_total + EPS:
            return

        min_remaining = compute_remaining_min_total(
            order,
            pos,
            mapping,
            used_product,
            reactant_orders,
            product_orders,
            reactant_symbols,
            product_indices_by_symbol,
        )
        if min_remaining is None:
            return
        if best_total is not None and current_total + min_remaining > best_total + EPS:
            return

        if pos == len(order):
            total = current_total
            if best_total is None or total < best_total - EPS:
                best_total = total
                best_pairs = {(break_count, add_count)}
                results.clear()
                seen_keys.clear()
            elif abs(total - best_total) <= EPS:
                best_pairs.add((break_count, add_count))
            else:
                return

            key = mapping_key(mapping, reactant_ranks, product_ranks)
            if key in seen_keys:
                return
            seen_keys.add(key)
            mapped = build_mapped_reaction_smiles(reactant_mol, product_mol, mapping)
            results.add(mapped)
            kept += 1
            if draw_path and (draw_limit <= 0 or drawn < draw_limit):
                output_path = build_draw_path(draw_path, drawn)
                save_reaction_svg(reactant_mol, product_mol, mapping, output_path)
                drawn += 1
            return

        react_idx = order[pos]
        symbol = reactant_symbols[react_idx]
        candidates = product_indices_by_symbol[symbol]
        for prod_idx in candidates:
            if used_product[prod_idx]:
                continue
            delta_add = 0.0
            delta_break = 0.0
            for assigned_idx, assigned_prod in mapping.items():
                delta = product_orders[prod_idx][assigned_prod] - reactant_orders[
                    react_idx
                ][assigned_idx]
                if delta > 0:
                    delta_add += delta
                elif delta < 0:
                    delta_break += -delta

            mapping[react_idx] = prod_idx
            used_product[prod_idx] = True
            backtrack(pos + 1, add_count + delta_add, break_count + delta_break)
            used_product[prod_idx] = False
            mapping.pop(react_idx, None)

    backtrack(0, 0.0, 0.0)
    if debug:
        elapsed = time.time() - start_time
        best_total_str = "none" if best_total is None else format_count(best_total)
        best_pairs_str = format_best_pairs(best_pairs)
        print(
            f"[debug] done nodes={nodes} kept={kept} drawn={drawn} "
            f"best_total={best_total_str} best_ops={best_pairs_str} "
            f"elapsed={elapsed:.1f}s",
            file=sys.stderr,
        )
    return sorted(results)


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Generate unique atom-mapped reaction SMILES with minimal add+break count."
        )
    )
    parser.add_argument("smiles", help="Reaction SMILES, e.g. O=CCCOO>>CC=O.O=CO")
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Print debug information and progress to stderr.",
    )
    parser.add_argument(
        "--debug-interval",
        type=int,
        default=100000,
        help="Progress interval when --debug is set, default 100000.",
    )
    parser.add_argument(
        "--draw",
        action="store_true",
        help="Save reaction SVGs for accepted mappings.",
    )
    parser.add_argument(
        "--draw-path",
        default="reaction.svg",
        help="Output SVG path or base name when --draw is set.",
    )
    parser.add_argument(
        "--draw-limit",
        type=int,
        default=0,
        help="Maximum number of SVGs to write when --draw is set; 0 for all.",
    )
    args = parser.parse_args()

    results = generate_atom_mapped_reactions(
        args.smiles,
        debug=args.debug,
        debug_interval=args.debug_interval,
        draw_path=args.draw_path if args.draw else None,
        draw_limit=args.draw_limit,
    )

    for rxn in results:
        print(rxn)


if __name__ == "__main__":
    main()
