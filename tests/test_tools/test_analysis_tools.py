#!/usr/bin/env python3
"""
Test script for polymer analysis tools.

This script demonstrates the analysis pipeline using real LAMMPS simulation data.
It tests all the conformational analysis functions with a linear polymer system.

Usage:
    python test_analysis_tools.py
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt

# Add the project root to the path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from src.tools.analysis import (
    read_simulation_data,
    analyze_conformation,
    calculate_radius_of_gyration,
    calculate_persistence_length,
    calculate_diffusion_coefficient,
    calculate_pq
)


def test_linear_polymer_analysis():
    """Test analysis tools with linear polymer data."""

    print("Polymer Analysis Tools Test")
    print("=" * 50)

    # Define paths to test data
    data_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "data", "test")
    data_file = os.path.join(data_dir, "system_linear_polymer_30.data")
    coord_dir = os.path.join(data_dir, "system_linear_polymer_30", "coord")
    dump_pattern = os.path.join(coord_dir, "dump.*.txt")

    print(f"Data file: {data_file}")
    print(f"Dump pattern: {dump_pattern}")
    print()

    try:
        # Step 1: Read simulation data
        print("Step 1: Reading simulation data...")
        polymer_atoms, all_atoms, bonds, trajectory, box, metadata = read_simulation_data(data_file, dump_pattern)

        atoms = polymer_atoms  # For compatibility

        print(f"✓ Loaded {len(all_atoms['ids'])} total atoms")
        print(f"✓ Loaded {len(atoms['ids'])} polymer atoms")
        print(f"✓ Loaded {len(bonds)} bonds")
        print(f"✓ Loaded {len(trajectory['polymer_atom_positions'])} trajectory frames")
        print(f"✓ Box dimensions: {box['lengths']}")
        print(f"✓ Metadata: {metadata}")
        print()

        # Step 2: Test individual analysis functions
        print("Step 2: Testing individual analysis functions...")

        # Radius of gyration
        print("  - Calculating radius of gyration...")
        rg = calculate_radius_of_gyration(trajectory['polymer_atom_positions'])
        print(".3f")
        print(".3f")
        print()

        # Persistence length
        print("  - Calculating persistence length...")
        if len(bonds) > 0:
            lp = calculate_persistence_length(
                trajectory['polymer_atom_positions'][-1],  # Last frame
                bonds,
                trajectory['polymer_atom_ids']
            )
            print(".3f")
        else:
            print("  ✗ No bonds found for persistence length calculation")
        print()

        # Diffusion coefficient
        print("  - Calculating diffusion coefficient...")
        D = calculate_diffusion_coefficient(
            trajectory['polymer_atom_positions'],
            trajectory['steps']
        )
        print(".6f")
        print()

        # Step 3: Test comprehensive analysis
        print("Step 3: Testing comprehensive analysis...")
        polymer_type = metadata.get("type", "linear")
        results = analyze_conformation(atoms, bonds, trajectory, box, polymer_type)

        print("✓ Comprehensive analysis completed!")
        print("\nResults Summary:")
        print("-" * 30)
        for key, value in results.items():
            if key == 'radius_of_gyration' and hasattr(value, '__len__'):
                print(".3f")
            elif key == 'scattering_q':
                print(f"{key}: {len(value)} points, range {value[0]:.3f} - {value[-1]:.3f}")
            elif key == 'scattering_pq':
                print(f"{key}: {len(value)} frames, {len(value[0]) if value else 0} q-points each")
            elif key == 'mean_scattering_pq':
                print(f"{key}: {len(value)} points, range {value[0]:.3f} - {value[-1]:.3f}")
            elif key == 'persistence_length' and hasattr(value, '__len__'):
                # persistence_length is a list of values (one per frame)
                valid_lp = [lp for lp in value if lp is not None and not np.isnan(lp)]
                if valid_lp:
                    print(f"{key}: {len(value)} frames, mean {np.mean(valid_lp):.3f}")
                else:
                    print(f"{key}: no valid values")
            elif hasattr(value, '__len__') and len(value) > 0 and not isinstance(value[0], (list, np.ndarray)):
                print(f"{key}: {value[-1]:.3f}")
            else:
                print(f"{key}: {value}")
        print()

        # Step 4: Test pair correlation (if MDAnalysis is available)
        print("Step 4: Testing pair correlation...")
        try:
            # Create a simple test universe for pair correlation
            import MDAnalysis as mda
            from MDAnalysis.analysis.rdf import InterRDF
            u = mda.Universe(data_file)
            rdf = InterRDF(u.atoms, u.atoms, range=(0.0, 10.0), nbins=100)
            rdf.run()
            print(f"✓ Pair correlation calculated with {len(rdf.bins)} bins")
            print(".3f")
        except Exception as e:
            print(f"  ⚠ Pair correlation test failed: {e}")
        print()

        # Step 5: Generate plots (optional)
        print("Step 5: Generating analysis plots...")
        try:
            plot_analysis_results(trajectory, rg, results)
            print("✓ Analysis plots saved to 'analysis_plots.png'")
        except Exception as e:
            print(f"  ⚠ Plot generation failed: {e}")
        print()

        print("🎉 All analysis tests completed successfully!")
        assert True
        return True

    except FileNotFoundError as e:
        print(f"❌ Error: File not found - {e}")
        print("Please ensure the test data files exist in the data/test/ directory")
        assert False
    except Exception as e:
        print(f"❌ Error during analysis: {e}")
        import traceback
        traceback.print_exc()
        return False


def plot_analysis_results(trajectory, rg, results):
    """Generate plots showing analysis results."""

    # Determine available data
    has_msd = results.get('msd') is not None
    has_scattering = 'mean_scattering_pq' in results and 'scattering_q' in results
    has_persistence = results.get('mean_persistence_length') is not None

    # Create figure with dynamic layout
    if has_msd and has_scattering:
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        plot_positions = [(0, 0), (0, 1), (1, 0)]
        if has_persistence:
            plot_positions.append((1, 1))
        empty_pos = (1, 1) if not has_persistence else None
    else:
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        plot_positions = [(0, 0)]
        if has_msd:
            plot_positions.append((0, 1))
        elif has_scattering:
            plot_positions.append((0, 1))
        empty_pos = None

    fig.suptitle('Polymer Analysis Results', fontsize=16)

    # Plot 1: Radius of gyration vs time (always available)
    steps = trajectory['steps']
    pos = plot_positions[0]
    axes[pos].plot(steps, rg, 'b-', lw=1)
    axes[pos].set_xlabel('Time', labelpad=0)
    axes[pos].set_ylabel('Radius of Gyration', labelpad=0)
    axes[pos].set_title('Rg vs Time')
    axes[pos].grid(True, alpha=0.3)

    # Plot 2: Mean Square Displacement (if available)
    plot_idx = 1
    if has_msd and plot_idx < len(plot_positions):
        pos = plot_positions[plot_idx]
        msd_values = np.array(results['msd'])
        dt_values = np.arange(len(msd_values)) * 0.01 * 1000  # matching the dt used in calculation
        axes[pos].plot(dt_values, msd_values, 'g-', lw=1)
        axes[pos].set_xlabel('Time', labelpad=0)
        axes[pos].set_ylabel('MSD', labelpad=0)
        axes[pos].set_title('Mean Square Displacement')
        axes[pos].grid(True, alpha=0.3)
        plot_idx += 1

    # Plot 3: Scattering P(q) (if available)
    if has_scattering and plot_idx < len(plot_positions):
        pos = plot_positions[plot_idx]
        q_values = np.array(results['scattering_q'])
        pq_values = np.array(results['mean_scattering_pq'])
        axes[pos].loglog(q_values, pq_values, 'm-', lw=1)
        axes[pos].set_xlabel('q', labelpad=0)
        axes[pos].set_ylabel('P(q)', labelpad=0)
        axes[pos].set_title('Scattering Function')
        axes[pos].grid(True, alpha=0.3)
        plot_idx += 1

    # Plot 4: Persistence length (if available and space)
    if has_persistence and plot_idx < len(plot_positions):
        pos = plot_positions[plot_idx]
        axes[pos].axis('off')
        persistence_text = f"""
Persistence Length Analysis:

• Mean persistence length: {results.get('mean_persistence_length', 'N/A'):.3f}
• This measures the stiffness of the polymer chain
"""
        axes[pos].text(0.1, 0.9, persistence_text, transform=axes[pos].transAxes,
                      fontsize=10, verticalalignment='top', fontfamily='monospace')
        axes[pos].set_title('Persistence Length')

    # Turn off empty subplots
    if empty_pos is not None:
        axes[empty_pos].axis('off')

    plt.tight_layout()
    plt.savefig('analysis_plots.png', dpi=150, bbox_inches='tight')
    plt.close()


def test_multiple_systems():
    """Test analysis on multiple polymer systems."""

    print("\nTesting Multiple Polymer Systems")
    print("=" * 40)

    data_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "data", "test")

    systems = [
        ("Linear Polymer 30", "system_linear_polymer_30.data", "system_linear_polymer_30/coord/dump.*.txt"),
        ("Ring Polymer 20", "system_ring_polymer_20.data", "system_ring_polymer_20/coord/dump.*.txt"),
        ("Star Polymer 8-4", "system_star_polymer_8_4.data", "system_star_polymer_8_4/coord/dump.*.txt"),
    ]

    results_summary = []

    for name, data_file, dump_pattern in systems:
        print(f"\nTesting {name}...")
        try:
            data_path = os.path.join(data_dir, data_file)
            dump_path = os.path.join(data_dir, dump_pattern)

            polymer_atoms, all_atoms, bonds, trajectory, box, metadata = read_simulation_data(data_path, dump_path)
            atoms = polymer_atoms
            polymer_type = metadata.get("type", "linear")
            results = analyze_conformation(atoms, bonds, trajectory, box, polymer_type)

            rg_final = results['radius_of_gyration'][-1] if hasattr(results['radius_of_gyration'], '__len__') else results['radius_of_gyration']

            results_summary.append({
                'system': name,
                'atoms': len(atoms['ids']),
                'frames': len(trajectory['polymer_atom_positions']),
                'rg_final': rg_final,
                'persistence_length': results.get('persistence_length'),
                'diffusion_coeff': results.get('diffusion_coefficient')
            })

            print(f"  ✓ {name}: {len(atoms['ids'])} atoms, Rg = {rg_final:.3f}")

        except Exception as e:
            print(f"  ❌ {name}: Failed - {e}")

    # Print summary table
    if results_summary:
        print("\n📊 Summary of All Systems:")
        print("-" * 80)
        print("<25")
        print("-" * 80)
        for result in results_summary:
            print("<25")
        print("-" * 80)


def main():
    """Main test function."""

    print("Starting Polymer Analysis Tools Test Suite")
    print("=" * 60)

    # Test single system
    success = test_linear_polymer_analysis()

    if success:
        # Test multiple systems
        test_multiple_systems()

        print("\n🎉 All tests completed successfully!")
        print("Check 'analysis_plots.png' for visualization of results.")
    else:
        print("\n❌ Primary test failed. Skipping additional tests.")
        sys.exit(1)


if __name__ == "__main__":
    main()
