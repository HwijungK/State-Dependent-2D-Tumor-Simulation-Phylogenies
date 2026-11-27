#!/usr/bin/env python

"""Write final-world plots and a sampled ancestry tree without opening a GUI."""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import pandas as pd

from reconstruct_sim_tree import read_sim_log, build_dendropy_tree, to_biopython_tree, plot_sim_tree


def visualize(simulation_log: Path, ancestry_log: Path, rates_file: Path, output_dir: Path) -> None:
    columns = ["id", "death_date", "states", "position_x", "position_y"]
    cells = pd.read_csv(simulation_log, sep="\t", usecols=columns, dtype={"id": "string"})
    cells = cells.loc[cells["death_date"].eq(cells["death_date"].max())].copy()
    rates = pd.read_csv(rates_file, sep="\t", dtype={"id": "string"})
    cells = cells.merge(rates, on="id")
    output_dir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(8, 7))
    for state, color in {"edge": "#E69F00", "center": "#0072B2"}.items():
        group = cells.loc[cells["states"].eq(state)]
        ax.scatter(group["position_x"], group["position_y"], c=color, s=35, label=state)
    ax.legend(title="Cell state", loc="upper left", bbox_to_anchor=(1.02, 1))
    ax.set_title(f"Final world ({len(cells)} surviving cells)")
    finish_world_plot(fig, ax, output_dir / "world_states.png")

    fig, ax = plt.subplots(figsize=(8, 7))
    points = ax.scatter(
        cells["position_x"], cells["position_y"],
        c=cells["lineage_mean_birthrate"], cmap="viridis", s=35,
    )
    fig.colorbar(points, ax=ax, label="Lineage mean birth rate (generation / tick)")
    ax.set_title("Final survivor lineage mean birth rates")
    finish_world_plot(fig, ax, output_dir / "world_lineage_mean_birthrate.png")

    tree = build_dendropy_tree(read_sim_log(ancestry_log))
    tree.write(path=str(output_dir / "sampled_tree.nwk"), schema="newick", suppress_annotations=False)
    rates_by_id = dict(zip(rates["id"], rates["lineage_mean_birthrate"]))
    fig, ax = plot_sim_tree(to_biopython_tree(tree), birth_rates=rates_by_id)
    legend = ax.get_legend()
    ax.legend(handles=legend.legend_handles, title="Cell state", loc="upper left", bbox_to_anchor=(1.02, 1))
    fig.axes[1].set_ylabel("Lineage mean birth rate\n(generation / tick)", fontsize=10)
    ax.set_xlabel("Logged simulation time")
    ax.set_title("Sampled ancestry colored by lineage mean birth rate")
    fig.savefig(output_dir / "sampled_tree.svg", bbox_inches="tight")
    fig.savefig(output_dir / "sampled_tree.png", dpi=140, bbox_inches="tight")
    plt.close(fig)


def finish_world_plot(fig, ax, output: Path) -> None:
    ax.set(xlabel="position_x", ylabel="position_y", aspect="equal")
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    fig.savefig(output, dpi=160, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--simulation-log", type=Path, required=True)
    parser.add_argument("--ancestry-log", type=Path, required=True)
    parser.add_argument("--rates", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    visualize(args.simulation_log, args.ancestry_log, args.rates, args.output_dir)


if __name__ == "__main__":
    main()
