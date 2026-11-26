"""Reconstruct and plot cell ancestry from a completed Eden simulation log."""

#!/usr/bin/env python

import re
from pathlib import Path
from typing import Literal

import dendropy
import matplotlib.pyplot as plt
import pandas as pd
from Bio import Phylo
from Bio.Phylo.BaseTree import Clade, Tree
from matplotlib.axes import Axes
from matplotlib.figure import Figure

import matplotlib as mpl
from matplotlib.lines import Line2D


def read_sim_log(path: str | Path) -> pd.DataFrame:
    """Read a simulation TSV, preserving cell and parent IDs as strings."""
    return pd.read_csv(
        path,
        sep="\t",
        dtype={
            "id": "string",
            "parent_id": "string",
            "birth_date": "int64",
            "death_date": "int64",
            "states": "string",
        },
    )


def build_dendropy_tree(
    df: pd.DataFrame,
    branch_length: Literal["time", "mutation"] = "time",
) -> dendropy.Tree:
    """Build a rooted tree with cell-state annotations on every node.

    ``df`` must contain the sampled cells and all ancestors back to one root.
    Only final survivors are labeled with cell IDs; they have the maximum
    ``death_date`` in the completed log. Every node records its logged state
    in a ``type`` annotation. Include these in Newick exports by passing
    ``suppress_annotations=False`` to ``tree.as_string`` or ``tree.write``.
    ``branch_length="time"`` uses each cell's ``death_date - birth_date``;
    nodes therefore represent the end of each cell's lifetime, and final
    survivors align at the simulation end time (relative to the root's birth).
    ``branch_length="mutation"`` uses one unit per generation, including the
    root, matching the simulator's generation count. It does not count observed
    sequence mutations. The root's incoming edge represents its own lifetime.
    """
    if branch_length not in ("time", "mutation"):
        raise ValueError("branch_length must be 'time' or 'mutation'.")

    required_columns = {"id", "parent_id", "birth_date", "death_date", "states"}
    missing_columns = required_columns - set(df.columns)
    if missing_columns:
        raise ValueError(f"Missing required columns: {', '.join(sorted(missing_columns))}.")
    if df.empty:
        raise ValueError("The simulation log is empty.")
    if df["id"].isna().any() or not df["id"].is_unique:
        raise ValueError("Each id must be unique and nonmissing.")
    if not df["states"].isin(("edge", "center")).all():
        raise ValueError("Each individual must have a states value of 'edge' or 'center'.")
    if df[["birth_date", "death_date"]].isna().any().any():
        raise ValueError("Each cell must have a birth_date and death_date.")
    if df["death_date"].lt(df["birth_date"]).any():
        raise ValueError("A cell's death_date cannot precede its birth_date.")

    root_ids = df.loc[df["parent_id"].isna(), "id"].tolist()
    if len(root_ids) != 1:
        raise ValueError("Expected exactly one root with a blank parent_id.")
    missing_parents = set(df["parent_id"].dropna()) - set(df["id"])
    if missing_parents:
        raise ValueError(f"Missing parent IDs: {', '.join(map(str, missing_parents))}.")

    # Validate connectivity before DendroPy can reject a cyclic attachment.
    children = {cell_id: [] for cell_id in df["id"]}
    for row in df.itertuples(index=False):
        if pd.notna(row.parent_id):
            children[row.parent_id].append(row.id)
    connected_ids = set()
    pending = [root_ids[0]]
    while pending:
        cell_id = pending.pop()
        connected_ids.add(cell_id)
        pending.extend(children[cell_id])
    if len(connected_ids) != len(df):
        raise ValueError("Some individuals do not connect to the root; check for cycles.")

    survivor_ids = set(df.loc[df["death_date"].eq(df["death_date"].max()), "id"])
    # Create every node first so the input row order does not matter.
    nodes = {}
    for row in df.itertuples(index=False):
        node = dendropy.Node(
            label=str(row.id) if row.id in survivor_ids else None,
            edge_length=row.death_date - row.birth_date if branch_length == "time" else 1,
        )
        node.annotations.add_new("type", row.states)
        nodes[row.id] = node
    for parent_id, child_ids in children.items():
        for child_id in child_ids:
            nodes[parent_id].add_child(nodes[child_id])

    return dendropy.Tree(seed_node=nodes[root_ids[0]], is_rooted=True)


def to_biopython_tree(tree: dendropy.Tree) -> Tree:
    """Convert to Biopython clades, preserving labels, lengths, and types.

    Direct conversion preserves numeric ID labels and zero-length edges.
    ``type`` annotations become Newick comments such as ``&type=edge``, which
    survive saving and loading with ``Bio.Phylo.write`` and ``Bio.Phylo.read``.
    The returned tree is ladderized for plotting; the input is unchanged.
    """
    clades = {
        node: Clade(
            name=node.label or (node.taxon.label if node.taxon else None),
            branch_length=node.edge_length,
        )
        for node in tree.preorder_node_iter()
    }
    for node, clade in clades.items():
        node_type = node.annotations.get_value("type")
        clade.comment = f"&type={node_type}" if node_type is not None else None
        clade.clades = [clades[child] for child in node.child_node_iter()]

    plot_tree = Tree(root=clades[tree.seed_node], rooted=tree.is_rooted)
    plot_tree.ladderize()
    return plot_tree


def plot_sim_tree(
    tree: Tree,
    *,
    ax: Axes | None = None,
    birth_rates: dict[str, float] | None = None,
    vmin: float | None = None,
    vmax: float | None = None,
) -> tuple[Figure, Axes]:
    """Plot state markers for survivors and their branching common ancestors.

    Ancestors receive a marker where at least two child lineages contain final
    survivors. Ancestor IDs stay hidden, and ancestors along a single retained
    lineage stay unmarked. Only survivors should have node names, as set by
    ``build_dendropy_tree``. Types come from each clade's Newick comment
    (``&type=edge`` or ``&type=center``); no simulation dataframe is needed.

    A new figure scales with the number of tips, or callers can supply their
    own axes. Returns ``(fig, ax)`` without calling ``plt.show()`` so callers
    can customize, display, or save the figure.
    """
    state_colors = {"edge": "#E69F00", "center": "#0072B2"}
    types_by_clade = {}
    for clade in tree.find_clades():
        match = re.search(r"(?:^&|,)\s*type\s*=\s*([^,]+)", getattr(clade, "comment", None) or "")
        node_type = match.group(1).strip().strip("\"'") if match else None
        if node_type not in state_colors:
            raise ValueError("Each node must have a type annotation of 'edge' or 'center'.")
        types_by_clade[clade] = node_type


    rates_by_clade = {}
    if birth_rates:
        for clade in tree.find_clades(order="postorder"):
            if clade.name is not None:
                rates_by_clade[clade] = birth_rates[clade.name]
            else:
                tip_rates = [birth_rates[t.name] for t in clade.get_terminals()
                             if t.name in birth_rates]
                rates_by_clade[clade] = sum(tip_rates) / len(tip_rates)

    use_rate = bool(rates_by_clade)
    if use_rate:
        norm = mpl.colors.Normalize(
            vmin=vmin if vmin is not None else min(rates_by_clade.values()),
            vmax=vmax if vmax is not None else max(rates_by_clade.values()),
        )
        cmap = plt.get_cmap("viridis")
        for clade, r in rates_by_clade.items():
            clade.color = mpl.colors.to_hex(cmap(norm(r)))
            clade.width = 1.8
    else:
        tree.root.color = "#237A8A"
        tree.root.width = 1.8

    if ax is None:
        fig, ax = plt.subplots(figsize=(10, max(4, 0.28 * tree.count_terminals())))
    else:
        fig = ax.figure

    Phylo.draw(
        tree,
        axes=ax,
        do_show=False,
        show_confidence=False,
        label_func=lambda clade: clade.name,
    )

    # Match Phylo.draw's coordinates, including its fallback for all-zero lengths.
    x_positions = tree.depths()
    if not max(x_positions.values()):
        x_positions = tree.depths(unit_branch_lengths=True)
    y_positions = {
        tip: index for index, tip in enumerate(tree.get_terminals(), start=1)
    }
    has_survivor = {}
    marked_clades = []
    for clade in tree.find_clades(order="postorder"):
        is_survivor = clade.name is not None
        surviving_branches = sum(has_survivor[child] for child in clade.clades)
        has_survivor[clade] = is_survivor or surviving_branches > 0
        if clade.clades:
            y_positions[clade] = (
                y_positions[clade.clades[0]] + y_positions[clade.clades[-1]]
            ) / 2
        if is_survivor or surviving_branches >= 2:
            marked_clades.append(clade)

    markers = {"edge": "o", "center": "s"}
    for state, marker in markers.items():
        clades = [c for c in marked_clades if types_by_clade[c] == state]
        if not clades:
            continue
        xs = [x_positions[c] for c in clades]
        ys = [y_positions[c] for c in clades]
        if use_rate:
            ax.scatter(xs, ys, c=[rates_by_clade[c] for c in clades],
                       cmap=cmap, norm=norm, marker=marker, s=36,
                       edgecolors="white", linewidths=0.6, zorder=4)
        else:
            ax.scatter(xs, ys, marker=marker, s=36, color=state_colors[state],
                       edgecolors="white", linewidths=0.6, zorder=4)

    if use_rate:
        sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
        fig.colorbar(sm, ax=ax, pad=0.02, shrink=0.4,
                     label="Mean birth rate (generation / tick)")
        handles = [Line2D([], [], marker=m, linestyle="", color="gray", label=s)
                   for s, m in markers.items()]
        ax.legend(handles=handles, title="Cell state", loc="upper right")
    else:
        ax.legend(title="Cell state", loc="upper right")
    return fig, ax
