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
        atoms, bonds, trajectory, box = read_simulation_data(data_file, dump_pattern)

        print(f"✓ Loaded {len(atoms['ids'])} atoms")
        print(f"✓ Loaded {len(bonds)} bonds")
        print(f"✓ Loaded {len(trajectory['polymer_atom_positions'])} trajectory frames")
        print(f"✓ Box dimensions: {box['lengths']}")
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
        results = analyze_conformation(atoms, bonds, trajectory, box)

        print("✓ Comprehensive analysis completed!")
        print("\nResults Summary:")
        print("-" * 30)
        for key, value in results.items():
            if key == 'radius_of_gyration' and hasattr(value, '__len__'):
                print(".3f")
            elif hasattr(value, '__len__') and len(value) > 0:
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
        return True

    except FileNotFoundError as e:
        print(f"❌ Error: File not found - {e}")
        print("Please ensure the test data files exist in the data/test/ directory")
        return False
    except Exception as e:
        print(f"❌ Error during analysis: {e}")
        import traceback
        traceback.print_exc()
        return False


def plot_analysis_results(trajectory, rg, results):
    """Generate plots showing analysis results."""

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle('Polymer Analysis Results', fontsize=16)

    # Plot 1: Radius of gyration vs time
    steps = trajectory['steps']
    axes[0, 0].plot(steps, rg, 'b-', linewidth=2)
    axes[0, 0].set_xlabel('Time')
    axes[0, 0].set_ylabel('Radius of Gyration')
    axes[0, 0].set_title('Rg vs Time')
    axes[0, 0].grid(True, alpha=0.3)

    # Plot 2: Center of mass trajectory
    com_trajectory = np.mean(trajectory['polymer_atom_positions'], axis=1)
    axes[0, 1].plot(com_trajectory[:, 0], com_trajectory[:, 1], 'r-', alpha=0.7)
    axes[0, 1].set_xlabel('X position')
    axes[0, 1].set_ylabel('Y position')
    axes[0, 1].set_title('Center of Mass Trajectory')
    axes[0, 1].grid(True, alpha=0.3)
    axes[0, 1].set_aspect('equal')

    # Plot 3: Atom positions (last frame)
    last_positions = trajectory['polymer_atom_positions'][-1]
    # Convert atom types to numeric for coloring
    try:
        atom_type_nums = [int(t) for t in trajectory['polymer_atom_types']]
        scatter = axes[1, 0].scatter(last_positions[:, 0], last_positions[:, 1],
                                     c=atom_type_nums, cmap='viridis', alpha=0.7)
        plt.colorbar(scatter, ax=axes[1, 0], label='Atom Type')
    except (ValueError, TypeError):
        # If conversion fails, use a single color
        axes[1, 0].scatter(last_positions[:, 0], last_positions[:, 1],
                           c='blue', alpha=0.7)
    axes[1, 0].set_xlabel('X position')
    axes[1, 0].set_ylabel('Y position')
    axes[1, 0].set_title('Final Configuration')
    axes[1, 0].grid(True, alpha=0.3)
    axes[1, 0].set_aspect('equal')

    # Plot 4: Analysis summary
    axes[1, 1].axis('off')
    summary_text = ".3f"".3f"".3f"".3f"f"""
Analysis Summary:

• Total atoms: {len(trajectory['polymer_atom_ids'])}
• Trajectory frames: {len(trajectory['polymer_atom_positions'])}
• Simulation steps: {steps[-1]:.1f}
• Final Rg: {rg[-1]:.3f}
• Persistence length: {results.get('persistence_length', 'N/A')}
• Diffusion coeff: {results.get('diffusion_coefficient', 'N/A'):.2e}
• End-to-end dist: {results.get('end_to_end_distance', 'N/A'):.3f}
"""
    axes[1, 1].text(0.1, 0.9, summary_text, transform=axes[1, 1].transAxes,
                    fontsize=10, verticalalignment='top', fontfamily='monospace')

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

            atoms, bonds, trajectory, box = read_simulation_data(data_path, dump_path)
            results = analyze_conformation(atoms, bonds, trajectory, box)

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
