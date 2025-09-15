#!/usr/bin/env python3
"""
Polymer Analysis Runner

This module provides functions to run complete conformational analysis
on LAMMPS simulation data and return structured results for plotting
and further analysis.
"""

import sys
import os
import numpy as np
from typing import Dict, Any, Tuple, Optional
import matplotlib.pyplot as plt

# Add the project root to the path
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(os.path.dirname(os.path.dirname(current_dir)))
sys.path.insert(0, project_root)

from src.tools.analysis import read_simulation_data, calculate_radius_of_gyration, calculate_persistence_length, calculate_diffusion_coefficient, calculate_pq


def run_complete_analysis(datafile_path: str, dump_pattern: str) -> Dict[str, Any]:
    """
    Run complete conformational analysis on LAMMPS simulation data.

    Parameters:
    - datafile_path: Path to LAMMPS data file
    - dump_pattern: Glob pattern for dump files (e.g., "coord/dump.*.txt")

    Returns:
    - results: Dictionary containing all analysis results
    """
    results = {"metadata": {"datafile_path": datafile_path, "dump_pattern": dump_pattern, "success": False, "error": None}, "system_info": {}, "analysis_results": {}}

    try:
        # Step 1: Read simulation data
        polymer_atoms, all_atoms, bonds, polymer_trajectory, box = read_simulation_data(datafile_path, dump_pattern)

        # Determine polymer type
        polymer_type = "linear"
        if "brush" in datafile_path:
            polymer_type = "brush"
        elif "star" in datafile_path:
            polymer_type = "star"
        elif "ring" in datafile_path:
            polymer_type = "ring"

        # Store system information
        results["system_info"] = {
            "n_atoms": len(all_atoms["ids"]),
            "n_polymer_atoms": len(polymer_atoms["ids"]),
            "n_bonds": len(bonds),
            "atom_types": list(set(all_atoms["types"])),
            "box_dimensions": box["lengths"].tolist(),
            "box_angles": box["angles"].tolist(),
            "polymer_type": polymer_type,
        }

        # Store trajectory information
        results["trajectory_info"] = {
            "n_frames": len(polymer_trajectory["polymer_atom_positions"]),
            "step_range": [int(polymer_trajectory["steps"][0]), int(polymer_trajectory["steps"][-1])],
            "step_size": int(polymer_trajectory["steps"][1] - polymer_trajectory["steps"][0]) if len(polymer_trajectory["steps"]) > 1 else 1,
        }

        # Step 2: Run conformation analysis

        # rg
        results["analysis_results"]["radius_of_gyration_trajectory"] = calculate_radius_of_gyration(polymer_trajectory["polymer_atom_positions"]).tolist()

        results["analysis_results"]["mean_radius_of_gyration"] = np.mean(results["analysis_results"]["radius_of_gyration_trajectory"])

        # persistence length
        lp_traj = []
        lp_chains_traj = []
        print("polymer_type", polymer_type)
        for frame_positions in polymer_trajectory["polymer_atom_positions"]:
            lp = calculate_persistence_length(frame_positions, bonds, polymer_trajectory["polymer_atom_ids"], polymer_type, polymer_atoms["types"])
            if isinstance(lp, dict):
                lp_traj.append(lp["overall"])
                lp_chains_traj.append(lp["chains"])
            else:
                lp_traj.append(lp)
                lp_chains_traj.append([lp])
        results["analysis_results"]["persistence_length_trajectory"] = lp_traj
        results["analysis_results"]["persistence_length_chains_trajectory"] = lp_chains_traj
        results["analysis_results"]["mean_persistence_length"] = np.mean([lp for lp in lp_traj if lp is not None and not np.isnan(lp)])

        # diffusion
        D, msd, dts = calculate_diffusion_coefficient(polymer_trajectory["polymer_atom_positions"], polymer_trajectory["steps"], 0.01 * 1000)  # matching lammps script
        results["analysis_results"]["msd"] = msd.tolist() if msd is not None else None
        results["analysis_results"]["msd_delt"] = 0.01 * 1000
        results["analysis_results"]["diffusion_coefficient"] = float(D)

        # scattering
        # 1. use Rg to find appropriate q range
        Rg = results["analysis_results"]["mean_radius_of_gyration"]
        q_value = np.logspace(-1, 1, 100) * 2 * np.pi / Rg
        # 2. calculate P(q)
        all_pq = calculate_pq(polymer_trajectory["polymer_atom_positions"], q_value)
        results["analysis_results"]["scattering_q"] = q_value.tolist()
        results["analysis_results"]["scattering_pq"] = all_pq.tolist()
        results["analysis_results"]["mean_scattering_pq"] = np.mean(all_pq, axis=0).tolist()

        # TODO: need to test the scattering calculation

        # Mark as successful
        results["metadata"]["success"] = True

    except Exception as e:
        results["metadata"]["error"] = str(e)
        results["metadata"]["success"] = False

    return results


def save_analysis_results(results: Dict[str, Any], output_file: str) -> None:
    """
    Save analysis results to a JSON file.

    Parameters:
    - results: Analysis results dictionary
    - output_file: Path to output JSON file
    """
    import json

    # Convert numpy arrays to lists for JSON serialization
    def convert_to_serializable(obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, (np.float64, np.float32)):
            return float(obj)
        elif isinstance(obj, (np.int64, np.int32)):
            return int(obj)
        elif isinstance(obj, np.bool_):
            return bool(obj)
        elif isinstance(obj, dict):
            return {k: convert_to_serializable(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [convert_to_serializable(item) for item in obj]
        elif isinstance(obj, tuple):
            return tuple(convert_to_serializable(item) for item in obj)
        else:
            return obj

    serializable_results = convert_to_serializable(results)

    with open(output_file, "w") as f:
        json.dump(serializable_results, f, indent=2)


def load_analysis_results(input_file: str) -> Dict[str, Any]:
    """
    Load analysis results from a JSON file.

    Parameters:
    - input_file: Path to input JSON file

    Returns:
    - results: Analysis results dictionary
    """
    import json

    with open(input_file, "r") as f:
        results = json.load(f)

    return results


def plot_analysis_results(data_dir, conformation_analysis_json: str):
    """Load analysis results and create plots."""
    import json

    # Load results
    with open(conformation_analysis_json, "r") as f:
        conformation_results = json.load(f)

    if not conformation_results["metadata"]["success"]:
        print(f"Analysis failed: {conformation_results['metadata']['error']}")
        return

    # 1 plot conformation results
    rg_traj = conformation_results["analysis_results"].get("radius_of_gyration_trajectory", None)
    mean_rg = conformation_results["analysis_results"].get("mean_radius_of_gyration", None)
    lp_traj = conformation_results["analysis_results"].get("persistence_length_trajectory", None)
    mean_lp = conformation_results["analysis_results"].get("mean_persistence_length", None)
    msd = conformation_results["analysis_results"].get("msd", None)
    msd_delt = conformation_results["analysis_results"].get("msd_delt", None)  # single float number
    t = np.arange(len(rg_traj)) * msd_delt if rg_traj is not None and msd_delt is not None else None

    print("msd_delt", msd_delt)

    # Create figure with subplots
    fig, axes = plt.subplots(1, 3, figsize=(3.3 * 2, 3.3 * 0.7))
    # fig.suptitle("Polymer Analysis Results", fontsize=16)

    # Plot 1: Radius of gyration trajectory
    if rg_traj is not None:
        axes[0].plot(t, rg_traj, label="tajectory")
        if mean_rg is not None:
            axes[0].axhline(mean_rg, color="r", linestyle="--", label=f"mean: {mean_rg:.3f}")
        axes[0].set_xlabel(r"$t$", fontsize=9)
        axes[0].set_ylabel(r"$R_g$", fontsize=9)
        axes[0].legend(fontsize=9)
        axes[0].grid(True, alpha=0.3)
        axes[0].tick_params(axis="both", which="both", direction="in", labelsize=7)

    # Plot 2: Persistence length trajectory (if available)
    if lp_traj is not None and any(lp is not None for lp in lp_traj):
        axes[1].plot(t, lp_traj, label="trajectory", color="g")
        if mean_lp is not None:
            axes[1].axhline(mean_lp, color="r", linestyle="--", label=f"mean: {mean_lp:.3f}")
        axes[1].set_xlabel(r"$t$", fontsize=9)
        axes[1].set_ylabel(r"$l_p$", fontsize=9)
        axes[1].legend(fontsize=9)
        axes[1].grid(True, alpha=0.3)
        axes[1].tick_params(axis="both", which="both", direction="in", labelsize=7)
    else:
        axes[1].text(0.5, 0.5, "No persistence length data available", transform=axes[1].transAxes, ha="center", va="center", fontsize=10)
        axes[1].axis("off")

    # Plot 3: diffusion
    if msd is not None and t is not None:
        axes[2].plot(t[: len(msd)], msd, "o", mfc="None", color="m", label="MSD")
        D = conformation_results["analysis_results"].get("diffusion_coefficient", "N/A")
        if D != "N/A":
            axes[2].plot(t[: len(msd)], t[: len(msd)] * D * 6, "--", color="gray", label="6D*t")
        axes[2].text(0.5, 0.9, f"D ~ {D:.2e}", transform=axes[2].transAxes, ha="center", va="center", fontsize=9)
        axes[2].legend(fontsize=9)
        axes[2].set_xlabel(r"$t$", fontsize=9)
        axes[2].set_ylabel(r"MSD", fontsize=9)
        axes[2].grid(True, alpha=0.3)
        axes[2].tick_params(axis="both", which="both", direction="in", labelsize=7)
    else:
        axes[2].text(0.5, 0.5, "No MSD data available", transform=axes[2].transAxes, ha="center", va="center", fontsize=9)
        axes[2].axis("off")

    plt.tight_layout(pad=0.5)
    plt.savefig(f"{data_dir}/conformation_analysis.png")
    plt.close()


def main():
    """Main function demonstrating the analysis workflow."""

    # Example usage - replace with your actual file paths
    data_dir = f"{project_root}/data/test/system_linear_polymer_30"
    datafile_path = f"{project_root}/data/test/system_linear_polymer_30.data"


    data_dir = f"{project_root}/data/test/system_brush_polymer_50_0.3_10"
    datafile_path = f"{project_root}/data/test/system_brush_polymer_50_0.3_10.data"


    data_dir = f"{project_root}/data/test/system_star_polymer_8_4"
    datafile_path = f"{project_root}/data/test/system_star_polymer_8_4.data"

    dump_pattern = os.path.join(data_dir, "coord", "dump.*.txt")

    print("Polymer Analysis Runner")
    print("=" * 50)

    # Run complete analysis
    results = run_complete_analysis(datafile_path, dump_pattern)

    # Optionally save results
    if results["metadata"]["success"]:
        output_file = f"{data_dir}/analysis_results.json"
        save_analysis_results(results, output_file)
        print(f"\n💾 Results saved to {output_file}")
    else:
        print(f"Analysis failed: {results['metadata']['error']}")

    # Plot results
    if results["metadata"]["success"]:
        conformation_analysis_json = f"{data_dir}/analysis_results.json"
        plot_analysis_results(data_dir, conformation_analysis_json)
        print(f"\n📊 Plots saved to {data_dir}/conformation_analysis.png")
    else:
        print(f"Analysis failed: {results['metadata']['error']}")

    return results


if __name__ == "__main__":
    main()
