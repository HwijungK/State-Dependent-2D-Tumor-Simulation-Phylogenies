#!/usr/bin/env python

import re
from pathlib import Path

import dendropy
import numpy as np
from dendropy.calculate import treecompare
from scipy.stats import spearmanr

SIM_DIR = Path(__file__).resolve().parent
SDEVO_FILE = SIM_DIR / "data" / "summary.tree"

def load_sdevo_tree(tns: dendropy.TaxonNamespace) -> dendropy.Tree:
    """Read SDevo's NEXUS tree, strip BEAST comments, and rename tips to bare IDs."""
    text = SDEVO_FILE.read_text()
    text = re.sub(r"\[[^\]]*\]", "", text)
    text = re.sub(r"(?:center|edge)_(\d+)_\d+", r"\1", text)
    return dendropy.Tree.get(
        data=text, schema="nexus",
        rooting="force-rooted", taxon_namespace=tns,
    )


def keep_tips_and_branch_points(tree: dendropy.Tree, name: str) -> None:
    """Keep only sampled tips and ancestors where lineages split."""
    tree.suppress_unifurcations()

    while len(tree.seed_node.child_nodes()) == 1:
        child = tree.seed_node.child_nodes()[0]
        child.parent_node = None
        tree.seed_node = child
    tree.seed_node.edge_length = None

    for node in tree.preorder_node_iter():
        if node.is_leaf():
            if node.taxon is None:
                raise ValueError(f"{name}: a tip has no label.")
        else:
            if len(node.child_nodes()) < 2:
                raise ValueError(f"{name}: unary internal node remains.")
            if node.taxon is not None:
                raise ValueError(
                    f"{name}: internal node labeled {node.taxon.label} is a sampled "
                    "cell with sampled descendants; it can't be treated as a tip."
                )
    tree.encode_bipartitions()
    print(f"{name}: {len(tree.leaf_nodes())} tips, "
          f"{len(tree.internal_nodes())} branching ancestors")


def tip_labels(tree: dendropy.Tree) -> set[str]:
    return {leaf.taxon.label for leaf in tree.leaf_node_iter()}


def scaled_pairwise_distances(tree: dendropy.Tree, labels: list[str]) -> np.ndarray:
    """Pairwise tip-to-tip distances, scaled to 0-1 so time units don't matter."""
    pdm = tree.phylogenetic_distance_matrix()
    taxa = {t.label: t for t in tree.taxon_namespace}
    values = np.array([
        pdm.distance(taxa[labels[i]], taxa[labels[j]])
        for i in range(len(labels))
        for j in range(i + 1, len(labels))
    ])
    return values / values.max()


def compare_to_sdevo(sim_tree_in: dendropy.Tree) -> None:
    """Compare a sim DendroPy tree with SDevo's summary tree."""
    if not SDEVO_FILE.exists():
        raise FileNotFoundError(f"SDevo tree not found: {SDEVO_FILE}")

    tns = dendropy.TaxonNamespace()
    sdevo_tree = load_sdevo_tree(tns)


    #work on a copy so notebook's tree isn't modified
    sim_tree = sim_tree_in.clone(depth=2)
    for leaf in sim_tree.leaf_node_iter():
        leaf.taxon = sim_tree.taxon_namespace.require_taxon(label=leaf.label)
        leaf.label = None
    sim_tree.migrate_taxon_namespace(tns)

    keep_tips_and_branch_points(sim_tree, "sim")
    keep_tips_and_branch_points(sdevo_tree, "SDevo")

    sim_tips, sdevo_tips = tip_labels(sim_tree), tip_labels(sdevo_tree)
    if sim_tips != sdevo_tips:
        raise ValueError(
            f"Tip labels differ. Only in sim: {sorted(sim_tips - sdevo_tips)}; "
            f"only in SDevo: {sorted(sdevo_tips - sim_tips)}"
        )
    n_tips = len(sim_tips)

    rf = treecompare.symmetric_difference(sim_tree, sdevo_tree)
    false_pos, false_neg = treecompare.false_positives_and_negatives(sim_tree, sdevo_tree)
    max_rf = 2 * (n_tips - 2)

    print(f"Tips compared: {n_tips}")
    print(f"RF distance: {rf}  (normalized: {rf / max_rf:.3f})")

    report_clade_differences(sim_tree, sdevo_tree)

    print(f"Clades only in SDevo: {false_pos}")
    print(f"Clades in sim missed by SDevo: {false_neg}")

    labels = sorted(sim_tips)
    sim_d = scaled_pairwise_distances(sim_tree, labels)
    sdevo_d = scaled_pairwise_distances(sdevo_tree, labels)
    rho, p_value = spearmanr(sim_d, sdevo_d)
    print(f"Pairwise-depth Spearman rho: {rho:.3f} (p = {p_value:.2g})")
    print(f"Pairwise-depth mean abs error (scaled): {np.mean(np.abs(sim_d - sdevo_d)):.3f}")



#these tell which groupings of cells differ btwn simulated and sdevo tree
def clade_sets(tree):
    """Return each internal clade as a frozenset of tip labels."""
    """essentially return family (descendant tips) of every internal node except root"""
    n_tips = len(tree.leaf_nodes())
    clades = set()
    for node in tree.preorder_node_iter():
        if node.is_leaf() or node is tree.seed_node:
            continue
        tips = frozenset(leaf.taxon.label for leaf in node.leaf_iter())
        if 1 < len(tips) < n_tips:
            clades.add(tips)
    return clades


def report_clade_differences(sim_tree, sdevo_tree):
    """compares two trees using clade_sets"""
    sim_clades = clade_sets(sim_tree)
    sdevo_clades = clade_sets(sdevo_tree)

    only_sdevo = sdevo_clades - sim_clades
    only_sim = sim_clades - sdevo_clades

    #print each clade's size, its sorted cell IDS (smallest first)
    print(f"Clades only in SDevo: {len(only_sdevo)}")
    for c in sorted(only_sdevo, key=len):
        print(f"  size {len(c)}: {sorted(c)}")

    print(f"Clades in sim missed by SDevo: {len(only_sim)}")
    for c in sorted(only_sim, key=len):
        print(f"  size {len(c)}: {sorted(c)}")

