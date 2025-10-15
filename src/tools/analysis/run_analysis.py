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

from src.tools.analysis import read_simulation_data, calculate_radius_of_gyration, calculate_persistence_length, calculate_diffusion_coefficient, calculate_pq, calculate_gr


def analyze_conformation(atoms, bonds, trajectory, box, polymer_type="linear"):
    """
    Analyze polymer conformation from simulation data structures.

    Parameters:
    - atoms: Polymer atoms dict
    - bonds: Bonds array
    - trajectory: Trajectory dict
    - box: Box dict
    - polymer_type: Type of polymer

    Returns:
    - results: Dict with analysis results
    """
    from src.tools.analysis.calculate_conformation import calculate_end_to_end_distance

    results = {}

    # Radius of gyration
    rg_traj = calculate_radius_of_gyration(trajectory['polymer_atom_positions'])
    results['radius_of_gyration'] = rg_traj
    results['mean_radius_of_gyration'] = np.mean(rg_traj)

    # End-to-end distance
    e2e_traj = []
    for frame_positions in trajectory['polymer_atom_positions']:
        e2e = calculate_end_to_end_distance(frame_positions, trajectory['polymer_atom_ids'], bonds, polymer_type, atoms['types'])
        e2e_traj.append(e2e)
    results['end_to_end_distance'] = e2e_traj
    valid_e2e = [e for e in e2e_traj if e is not None]
    results['mean_end_to_end_distance'] = np.mean(valid_e2e) if valid_e2e else None

    # Persistence length
    lp_traj = []
    for frame_positions in trajectory['polymer_atom_positions']:
        lp = calculate_persistence_length(frame_positions, bonds, trajectory['polymer_atom_ids'], polymer_type, atoms['types'])
        if isinstance(lp, dict):
            lp_traj.append(lp["overall"])
        else:
            lp_traj.append(lp)
    results['persistence_length'] = lp_traj
    results['mean_persistence_length'] = np.mean([lp for lp in lp_traj if lp is not None and not np.isnan(lp)])

    # Diffusion coefficient
    D, msd, dts = calculate_diffusion_coefficient(trajectory['polymer_atom_positions'], trajectory['steps'], 0.01 * 1000)
    results['diffusion_coefficient'] = float(D)
    results['msd'] = msd.tolist() if msd is not None else None

    # Scattering
    Rg = results['mean_radius_of_gyration']
    q_value = np.logspace(-1, 1, 100) * 2 * np.pi / Rg
    all_pq = calculate_pq(trajectory['polymer_atom_positions'], q_value)
    results['scattering_q'] = q_value.tolist()
    results['scattering_pq'] = all_pq.tolist()
    results['mean_scattering_pq'] = np.mean(all_pq, axis=0).tolist()

    return results


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

    data_dir = os.path.dirname(datafile_path)

    try:
        from src.tools.analysis.calculate_conformation import calculate_end_to_end_distance

        # Step 1: Read simulation data
        polymer_atoms, all_atoms, bonds, polymer_trajectory, box, metadata = read_simulation_data(datafile_path, dump_pattern)

        # Determine polymer type
        polymer_type = "linear"
        if "brush" in datafile_path:
            polymer_type = "brush"
        elif "star" in datafile_path:
            polymer_type = "star"
        elif "ring" in datafile_path:
            polymer_type = "ring"
        elif "dendrimer" in datafile_path:
            polymer_type = "dendrimer"

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

        n_frames = results["trajectory_info"]["n_frames"]
        # rg
        results["analysis_results"]["radius_of_gyration_trajectory"] = calculate_radius_of_gyration(polymer_trajectory["polymer_atom_positions"]).tolist()

        #results["analysis_results"]["mean_radius_of_gyration"] = np.mean(results["analysis_results"]["radius_of_gyration_trajectory"])
        results["analysis_results"]["mean_radius_of_gyration"] = np.mean(results["analysis_results"]["radius_of_gyration_trajectory"][int(n_frames/2):])

        # end-to-end distance
        end_to_end_traj = []
        end_to_end_chains_traj = []
        for frame_positions in polymer_trajectory["polymer_atom_positions"]:
            e2e = calculate_end_to_end_distance(frame_positions, polymer_trajectory["polymer_atom_ids"], bonds, polymer_type, polymer_atoms["types"])
            if isinstance(e2e, dict):
                end_to_end_traj.append(e2e["overall"])
                end_to_end_chains_traj.append(e2e["chains"])
            else:
                end_to_end_traj.append(e2e)
                end_to_end_chains_traj.append([e2e])
        results["analysis_results"]["end_to_end_distance_trajectory"] = end_to_end_traj
        results["analysis_results"]["end_to_end_distance_chains_trajectory"] = end_to_end_chains_traj
        valid_e2e_second_half = [e for i, e in enumerate(end_to_end_traj) if i >= int(n_frames/2) and e is not None and not np.isnan(e)]
        results["analysis_results"]["mean_end_to_end_distance"] = np.mean(valid_e2e_second_half) if valid_e2e_second_half else None

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
        valid_lp_second_half = [lp for i, lp in enumerate(lp_traj) if i >= int(n_frames/2) and lp is not None and not np.isnan(lp)]
        results["analysis_results"]["mean_persistence_length"] = np.mean(valid_lp_second_half) if valid_lp_second_half else None

        # diffusion
        half_idx = int(n_frames / 2)
        positions_second_half = polymer_trajectory["polymer_atom_positions"][half_idx:]
        steps_second_half = polymer_trajectory["steps"][half_idx:]
        D, msd, dts = calculate_diffusion_coefficient(positions_second_half, steps_second_half, 0.01 * 1000)  # matching lammps script
        results["analysis_results"]["msd"] = msd.tolist() if msd is not None else None
        results["analysis_results"]["msd_delt"] = 0.01 * 1000
        results["analysis_results"]["diffusion_coefficient"] = float(D)

        # scattering
        # 1. use Rg2 to find appropriate q range
        Rg2 = results["analysis_results"]["mean_radius_of_gyration"]
        Rg = np.sqrt(Rg2)
        q_value = np.logspace(-1.5, 1.2, 80) * 2 * np.pi / Rg
        # 2. calculate P(q)
        all_pq = calculate_pq(polymer_trajectory["polymer_atom_positions"], q_value)
        results["analysis_results"]["scattering_q"] = q_value.tolist()
        results["analysis_results"]["scattering_pq"] = all_pq.tolist()
        results["analysis_results"]["mean_scattering_pq"] = np.mean(all_pq[int(n_frames/2):], axis=0).tolist()

        # 3. calculate g(r)
        box_length = box["lengths"][0]  # assuming cubic box
        r_bins = np.linspace(0, box_length / 4, 100)
        all_gr = calculate_gr(polymer_trajectory["polymer_atom_positions"], r_bins, box_length)
        r_mid = (r_bins[:-1] + r_bins[1:]) / 2
        results["analysis_results"]["radial_distribution_r"] = r_mid.tolist()
        results["analysis_results"]["radial_distribution_gr"] = all_gr.tolist()
        results["analysis_results"]["mean_radial_distribution_gr"] = np.mean(all_gr[int(n_frames/2):], axis=0).tolist()

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
    e2e_traj = conformation_results["analysis_results"].get("end_to_end_distance_trajectory", None)
    mean_e2e = conformation_results["analysis_results"].get("mean_end_to_end_distance", None)
    msd = conformation_results["analysis_results"].get("msd", None)
    msd_delt = conformation_results["analysis_results"].get("msd_delt", None)  # single float number
    scattering_q = conformation_results["analysis_results"].get("scattering_q", None)
    mean_scattering_pq = conformation_results["analysis_results"].get("mean_scattering_pq", None)
    radial_distribution_r = conformation_results["analysis_results"].get("radial_distribution_r", None)
    mean_radial_distribution_gr = conformation_results["analysis_results"].get("mean_radial_distribution_gr", None)
    t = np.arange(len(rg_traj)) * msd_delt if rg_traj is not None and msd_delt is not None else None

    print("msd_delt", msd_delt)
    print(f"Data availability - rg_traj: {rg_traj is not None}, msd: {msd is not None}, scattering_q: {scattering_q is not None}, radial_distribution_r: {radial_distribution_r is not None}")

    # Check if essential data is available
    if rg_traj is None or msd is None or scattering_q is None or mean_scattering_pq is None or radial_distribution_r is None or mean_radial_distribution_gr is None:
        print("❌ Essential data missing for plotting")
        return

    # Determine available data
    has_persistence = lp_traj is not None and any(lp is not None for lp in lp_traj)
    has_e2e = e2e_traj is not None and any(e is not None for e in e2e_traj)

    print(f"Optional data - persistence: {has_persistence}, e2e: {has_e2e}")

    # Create figure with 3x2 layout
    fig, axes = plt.subplots(3, 2, figsize=(3.3, 3.3 * 1.5))
    plot_positions = [(0, 0), (0, 1), (1, 0), (1, 1), (2, 0), (2, 1)]
    empty_positions = []

    print(f"Using 3x2 layout with {len(plot_positions)} potential plots")

    try:
        # Plot 1: Radius of gyration trajectory (always available)
        pos = plot_positions[0]
        axes[pos].plot(t, rg_traj, lw=1)
        if mean_rg is not None:
            axes[pos].axhline(mean_rg, color="r", linestyle="--", label=f"mean: {mean_rg:.3f}")
        axes[pos].set_xlabel(r"$t$", fontsize=9, labelpad=0)
        axes[pos].set_ylabel(r"$R_g^2$", fontsize=9, labelpad=0)
        axes[pos].legend(fontsize=7, frameon=True, ncol=1, columnspacing=0.5, handlelength=1, handletextpad=0.2)
        axes[pos].grid(True, alpha=0.3)
        axes[pos].tick_params(axis="both", which="both", direction="in", labelsize=7, pad=1)
        axes[pos].set_title("Radius of gyration", fontsize=9, pad=3)

        # Plot 2: Diffusion/MSD (always available)
        pos = plot_positions[1]
        axes[pos].plot(t[: len(msd)], msd, "o", mfc="None", color="m", lw=1)
        D = conformation_results["analysis_results"].get("diffusion_coefficient", "N/A")
        if D != "N/A":
            axes[pos].plot(t[: len(msd)], t[: len(msd)] * D * 6, "--", color="gray", label="6D*t"+"\n"+f"D ~ {D:.2e}", lw=1)
        axes[pos].legend(fontsize=7, frameon=True, ncol=1, columnspacing=0.5, handlelength=1, handletextpad=0.2)
        axes[pos].set_xlabel(r"$\Delta t$", fontsize=9, labelpad=0)
        axes[pos].set_ylabel(r"MSD", fontsize=9, labelpad=0)
        axes[pos].grid(True, alpha=0.3)
        axes[pos].tick_params(axis="both", which="both", direction="in", labelsize=7, pad=1)
        axes[pos].set_title("Diffusion", fontsize=9, pad=3)

        # Plot 3: Scattering (P(q)) (always available)
        pos = plot_positions[2]
        axes[pos].loglog(scattering_q, mean_scattering_pq, 'b-', lw=1)
        axes[pos].set_xlabel(r"$q$", fontsize=9, labelpad=0)
        axes[pos].set_ylabel(r"$P(q)$", fontsize=9, labelpad=0)
        axes[pos].grid(True, alpha=0.3)
        axes[pos].tick_params(axis="both", which="both", direction="in", labelsize=7, pad=1)
        axes[pos].set_title("Scattering Function", fontsize=9, pad=3)

        # Plot 4: Radial distribution function g(r) (always available)
        pos = plot_positions[3]
        axes[pos].plot(radial_distribution_r, mean_radial_distribution_gr, 'r-', lw=1)
        axes[pos].set_xlabel(r"$r$", fontsize=9, labelpad=0)
        axes[pos].set_ylabel(r"$g(r)$", fontsize=9, labelpad=0)
        axes[pos].grid(True, alpha=0.3)
        axes[pos].tick_params(axis="both", which="both", direction="in", labelsize=7, pad=1)
        axes[pos].set_title("Radial Distribution", fontsize=9, pad=3)

        # Plot optional data if available
        plot_idx = 4
        if has_persistence and plot_idx < len(plot_positions):
            pos = plot_positions[plot_idx]
            axes[pos].plot(t, lp_traj, color="g", lw=1)
            if mean_lp is not None:
                axes[pos].axhline(mean_lp, color="r", linestyle="--", label=f"mean: {mean_lp:.3f}")
            axes[pos].set_xlabel(r"$t$", fontsize=9, labelpad=0)
            axes[pos].set_ylabel(r"$l_p$", fontsize=9, labelpad=0)
            axes[pos].legend(fontsize=7, frameon=True, ncol=1, columnspacing=0.5, handlelength=1, handletextpad=0.2)
            axes[pos].grid(True, alpha=0.3)
            axes[pos].tick_params(axis="both", which="both", direction="in", labelsize=7, pad=1)
            axes[pos].set_title("Persistence length", fontsize=9, pad=3)
            plot_idx += 1
        else:
            empty_positions.append(plot_positions[plot_idx])
            plot_idx += 1

        if has_e2e and plot_idx < len(plot_positions):
            pos = plot_positions[plot_idx]
            valid_indices = [i for i, e in enumerate(e2e_traj) if e is not None]
            valid_t = [t[i] for i in valid_indices] if t is not None else valid_indices
            valid_e2e = [e2e_traj[i] for i in valid_indices]
            axes[pos].plot(valid_t, valid_e2e, color="orange", lw=1)
            if mean_e2e is not None:
                axes[pos].axhline(mean_e2e, color="r", linestyle="--", label=f"mean: {mean_e2e:.3f}")
            axes[pos].set_xlabel(r"$t$", fontsize=9, labelpad=0)
            axes[pos].set_ylabel(r"$R_{ee}$", fontsize=9, labelpad=0)
            axes[pos].legend(fontsize=7, frameon=True, ncol=1, columnspacing=0.5, handlelength=1, handletextpad=0.2)
            axes[pos].grid(True, alpha=0.3)
            axes[pos].tick_params(axis="both", which="both", direction="in", labelsize=7, pad=1)
            axes[pos].set_title("End-to-end distance", fontsize=9, pad=3)
        else:
            empty_positions.append(plot_positions[plot_idx])

        # Turn off empty subplots
        for pos in empty_positions:
            axes[pos].axis("off")

        plt.tight_layout(pad=0.1)
        plt.savefig(f"{data_dir}/conformation_analysis.png", dpi=500, bbox_inches='tight')
        plt.close()
        print(f"✅ Plot saved successfully to {data_dir}/conformation_analysis.png")

    except Exception as e:
        print(f"❌ Error during plotting: {e}")
        import traceback
        traceback.print_exc()


def main():
    """Main function demonstrating the analysis workflow."""

    # Example usage - replace with your actual file paths
    # data_dir = f"{project_root}/data/test/system_linear_polymer"
    # datafile_path = f"{project_root}/data/test/system_linear_polymer.data"

    # data_dir = f"{project_root}/data/test/system_brush_polymer"
    # datafile_path = f"{project_root}/data/test/system_brush_polymer.data"

    # Use the ring polymer data that exists
    data_dir = f"{project_root}/data/test/brush_workflow_50_0.2_10_20250919_160223/system_brush"
    datafile_path = f"{project_root}/data/test/brush_workflow_50_0.2_10_20250919_160223/system_brush.data"

    print(f"Looking for data in: {data_dir}")
    print(f"Data file: {datafile_path}")

    if not os.path.exists(datafile_path):
        print(f"❌ Data file not found: {datafile_path}")
        return None

    dump_pattern = os.path.join(data_dir, "coord", "dump.*.txt")
    print(f"Dump pattern: {dump_pattern}")

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
        print(f"❌ Analysis failed: {results['metadata']['error']}")
        return results

    # Plot results
    if results["metadata"]["success"]:
        conformation_analysis_json = f"{data_dir}/analysis_results.json"
        if os.path.exists(conformation_analysis_json):
            plot_analysis_results(data_dir, conformation_analysis_json)
            print(f"\n📊 Plots saved to {data_dir}/conformation_analysis.png")
        else:
            print(f"❌ Analysis results file not found: {conformation_analysis_json}")
    else:
        print(f"❌ Analysis failed: {results['metadata']['error']}")

    return results


if __name__ == "__main__":
    main()
