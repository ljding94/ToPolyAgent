#!/usr/bin/env python3
"""
Full Workflow Script for Polymer Simulation and Analysis

This script runs the complete workflow for polymer simulations:
1. Generate polymer configuration (linear, ring, brush, star, or dendrimer)
2. Pack solvent around the polymer
3. Run LAMMPS molecular dynamics simulation
4. Analyze conformational properties
5. Generate analysis plots

The script creates a timestamped output directory containing:
- polymer_<type>.data: Initial polymer configuration
- system_<type>.data: System with solvent
- simulation outputs (dump files, log, final config)
- analysis_results.json: Detailed analysis results
- conformation_analysis.png: Analysis plots
- workflow_summary.json: Summary of the entire workflow

Usage:
    python test_full_workflow.py

Customization:
    The script defines type-specific config dictionaries in main(). Each config includes:
    - 'polymer_type': "linear", "ring", "brush", "star", or "dendrimer"
    - Common params: 'box_size', 'solvent_density', 'run_steps', 'thermostat', 'interaction_params'
    - 'polymer_params': dict of type-specific parameters

    To test a specific polymer type, modify the test_configs dict or run selectively.

Expected runtime: ~5-15 minutes depending on parameters and system performance
"""

import os
import sys
import json
from datetime import datetime

# Add project paths to sys.path


def find_project_root(start_path):
    current_path = start_path
    while current_path != os.path.dirname(current_path):
        if os.path.exists(os.path.join(current_path, 'setup.py')):
            return current_path
        current_path = os.path.dirname(current_path)
    raise ValueError("Could not find project root with setup.py")


project_root = find_project_root(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, project_root)
sys.path.insert(0, os.path.join(project_root, 'src'))


# Import tools
from tools.config_gen.generate_polymer import (
    generate_brush_polymer_config,
    generate_linear_polymer_config,
    generate_ring_polymer_config,
    generate_star_polymer_config,
    generate_dendrimer_config
)
from tools.config_gen.pack_solvent import pack_solvent
from tools.config_gen.plot_config import plot_configuration
from tools.simulation.run_lammps import run_lammps
from tools.analysis.run_analysis import run_complete_analysis, plot_analysis_results, save_analysis_results


def run_polymer_workflow(config, output_dir=None):
    """
    Run the complete polymer workflow using a config dict.

    Parameters:
    - config: dict with keys:
        - 'polymer_type': str ("linear", "ring", "brush", "star", "dendrimer")
        - Common params: 'box_size' (float, default 40.0), 'solvent_density' (float, default 0.1),
                         'run_steps' (int, default 50000), 'thermostat' (str, default "langevin"),
                         'interaction_params' (dict, default {"pp": 1.0, "ss": 0.8, "sp": 0.5})
        - 'polymer_params': dict of type-specific params (e.g., for 'brush': {'backbone_length': int, ...})
    - output_dir: str, output directory (auto-generated if None)

    Returns:
    - dict: workflow results and paths
    """

    print("=" * 60)
    print("POLYMER FULL WORKFLOW")
    print("=" * 60)
    print(f"Started at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print()

    # Extract and validate common params with defaults
    polymer_type = config['polymer_type']
    box_size = config.get('box_size', 40.0)
    solvent_density = config.get('solvent_density', 0.1)
    run_steps = config.get('run_steps', 50000)
    thermostat = config.get('thermostat', "langevin")
    interaction_params = config.get('interaction_params', {"pp": 1.0, "ss": 0.8, "sp": 0.5})

    # Extract type-specific params
    polymer_params = config.get('polymer_params', {})
    required_params = _get_required_params(polymer_type)  # Helper function, defined below
    for param in required_params:
        if param not in polymer_params:
            raise ValueError(f"Missing required param '{param}' for polymer_type '{polymer_type}'")

    # Create output directory (customize name based on type and key params)
    if output_dir is None:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        dir_suffix = _get_dir_suffix(polymer_type, polymer_params)  # Helper function, defined below
        output_dir = os.path.join(project_root, 'data', 'test', f"{polymer_type}_workflow_{dir_suffix}_{timestamp}")

    output_dir = os.path.abspath(output_dir)
    os.makedirs(output_dir, exist_ok=True)
    os.chdir(output_dir)

    print(f"Output directory: {os.path.abspath(output_dir)}")
    print()

    workflow_results = {
        "parameters": config,  # Store full config for reference
        "steps": {},
        "output_paths": {}
    }

    try:
        # Step 1: Generate polymer configuration
        print("Step 1: Generating polymer configuration...")
        print(f"  Polymer type: {polymer_type}")

        if polymer_type == "brush":
            print(f"  Backbone length: {polymer_params['backbone_length']}")
            print(f"  Grafting density: {polymer_params['grafting_density']}")
            print(f"  Side chain length: {polymer_params['side_chain_length']}")
            polymer_file = generate_brush_polymer_config(
                backbone_length=polymer_params['backbone_length'],
                grafting_density=polymer_params['grafting_density'],
                side_chain_length=polymer_params['side_chain_length'],
                box_size=box_size
            )
        elif polymer_type == "linear":
            print(f"  Chain length: {polymer_params['chain_length']}")
            polymer_file = generate_linear_polymer_config(
                chain_length=polymer_params['chain_length'],
                box_size=box_size
            )
        elif polymer_type == "ring":
            print(f"  Ring length: {polymer_params['ring_length']}")
            polymer_file = generate_ring_polymer_config(
                ring_length=polymer_params['ring_length'],
                box_size=box_size
            )
        elif polymer_type == "star":
            print(f"  Arm length: {polymer_params['arm_length']}")
            print(f"  Number of arms: {polymer_params['num_arms']}")
            polymer_file = generate_star_polymer_config(
                arm_length=polymer_params['arm_length'],
                num_arms=polymer_params['num_arms'],
                box_size=box_size
            )
        elif polymer_type == "dendrimer":
            print(f"  Generations: {polymer_params['generations']}")
            print(f"  Branching factor: {polymer_params['branching_factor']}")
            print(f"  Spacer: {polymer_params['spacer']}")
            polymer_file = generate_dendrimer_config(
                generations=polymer_params['generations'],
                branching_factor=polymer_params['branching_factor'],
                spacer=polymer_params['spacer'],
                box_size=box_size
            )
        else:
            raise ValueError(f"Unsupported polymer type: {polymer_type}")

        print(f"  ✓ Polymer configuration saved to: {polymer_file}")
        workflow_results["steps"]["config_generation"] = "success"
        workflow_results["output_paths"]["polymer_config"] = polymer_file

        # Plot polymer configuration
        print("  Plotting polymer configuration...")
        plot_configuration(polymer_file)
        polymer_plot_path = os.path.splitext(polymer_file)[0] + '.png'
        workflow_results["output_paths"]["polymer_plot"] = polymer_plot_path

        # Step 2: Pack solvent
        print("\nStep 2: Packing solvent around polymer...")
        print(f"  Solvent density: {solvent_density}")

        system_file = pack_solvent(
            polymer_datafile=polymer_file,
            solvent_density=solvent_density,
            box_size=box_size
        )

        print(f"  ✓ System configuration saved to: {system_file}")
        workflow_results["steps"]["solvent_packing"] = "success"
        workflow_results["output_paths"]["system_config"] = system_file

        # Plot system configuration
        print("  Plotting system configuration...")
        plot_configuration(system_file)
        system_plot_path = os.path.splitext(system_file)[0] + '.png'
        workflow_results["output_paths"]["system_plot"] = system_plot_path

        # Step 3: Run LAMMPS simulation
        print("\nStep 3: Running LAMMPS simulation...")
        print(f"  Simulation steps: {run_steps}")
        print(f"  Thermostat: {thermostat}")
        print(f"  Interaction parameters: {interaction_params}")

        dump_path = system_file.replace(".data", "")
        simulation_results = run_lammps(
            dump_path=dump_path,
            datafile_path=system_file,
            thermostat=thermostat,
            interaction_params=interaction_params,
            run_steps=run_steps
        )

        print("  ✓ Simulation completed")
        print(f"  ✓ Dump files: {simulation_results['dump_files']}")
        print(f"  ✓ Final config: {simulation_results['final_config']}")
        print(f"  ✓ Log file: {simulation_results['log']}")

        workflow_results["steps"]["simulation"] = "success"
        workflow_results["output_paths"]["simulation"] = simulation_results

        # Plot final configuration
        print("  Plotting final configuration...")
        final_config_path = simulation_results['final_config']

        # Plot with solvent
        plot_configuration(final_config_path, plot_solvent=True)
        final_plot_with_solvent = os.path.splitext(final_config_path)[0] + '.png'

        # Plot without solvent for better polymer visualization
        plot_configuration(final_config_path, plot_solvent=False)
        final_plot_without_solvent = os.path.splitext(final_config_path)[0] + '_nosolvent.png'

        workflow_results["output_paths"]["final_plot"] = final_plot_with_solvent
        workflow_results["output_paths"]["final_plot_nosolvent"] = final_plot_without_solvent

        # Step 4: Analyze results
        print("\nStep 4: Analyzing simulation results...")

        datafile_path = system_file
        dump_pattern = simulation_results['dump_files']

        analysis_results = run_complete_analysis(
            datafile_path=datafile_path,
            dump_pattern=dump_pattern
        )

        if analysis_results["metadata"]["success"]:
            print("  ✓ Analysis completed successfully")

            # Determine simulation subfolder (where final_state.data is saved)
            simulation_subfolder = os.path.dirname(simulation_results['final_config'])

            # Save analysis results in the simulation subfolder
            analysis_file = os.path.join(simulation_subfolder, "analysis_results.json")
            save_analysis_results(analysis_results, analysis_file)

            print(f"  ✓ Analysis results saved to: {analysis_file}")
            workflow_results["steps"]["analysis"] = "success"
            workflow_results["output_paths"]["analysis"] = analysis_file

            # Print key results
            print("\n  Key Analysis Results:")
            print(".3f")
            if "mean_persistence_length" in analysis_results["analysis_results"]:
                print(".3f")
            if "diffusion_coefficient" in analysis_results["analysis_results"]:
                print(".2e")

        else:
            print(f"  ✗ Analysis failed: {analysis_results['metadata']['error']}")
            workflow_results["steps"]["analysis"] = "failed"
            workflow_results["error"] = analysis_results["metadata"]["error"]

        # Step 5: Generate plots
        print("\nStep 5: Generating analysis plots...")

        if analysis_results["metadata"]["success"]:
            # Use simulation subfolder for plotting
            simulation_subfolder = os.path.dirname(simulation_results['final_config'])
            plot_analysis_results(simulation_subfolder, analysis_file)
            plot_file = os.path.join(simulation_subfolder, "conformation_analysis.png")
            print(f"  ✓ Plots saved to: {plot_file}")
            workflow_results["steps"]["plotting"] = "success"
            workflow_results["output_paths"]["plots"] = plot_file
        else:
            print("  ✗ Plotting skipped due to analysis failure")
            workflow_results["steps"]["plotting"] = "skipped"

        # Save workflow summary
        workflow_file = os.path.join(output_dir, "workflow_summary.json")
        with open(workflow_file, 'w') as f:
            json.dump(workflow_results, f, indent=2)

        print(f"\n✓ Workflow summary saved to: {workflow_file}")

        print("\n" + "=" * 60)
        print("WORKFLOW COMPLETED SUCCESSFULLY!")
        print("=" * 60)
        print(f"Completed at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print(f"Total output directory: {os.path.abspath(output_dir)}")

        workflow_results["status"] = "completed"

    except Exception as e:
        print(f"\n✗ Workflow failed with error: {str(e)}")
        workflow_results["status"] = "failed"
        workflow_results["error"] = str(e)

        # Save error summary
        workflow_file = os.path.join(output_dir, "workflow_summary.json")
        with open(workflow_file, 'w') as f:
            json.dump(workflow_results, f, indent=2)

    return workflow_results


# Helper: Get required params for each type
def _get_required_params(polymer_type):
    if polymer_type == "brush":
        return ['backbone_length', 'grafting_density', 'side_chain_length']
    elif polymer_type == "linear":
        return ['chain_length']
    elif polymer_type == "ring":
        return ['ring_length']
    elif polymer_type == "star":
        return ['arm_length', 'num_arms']
    elif polymer_type == "dendrimer":
        return ['generations', 'branching_factor', 'spacer']
    else:
        raise ValueError(f"Unsupported polymer type: {polymer_type}")


# Helper: Get suffix for output_dir naming
def _get_dir_suffix(polymer_type, polymer_params):
    if polymer_type == "brush":
        return f"{polymer_params['backbone_length']}_{polymer_params['grafting_density']}_{polymer_params['side_chain_length']}"
    elif polymer_type == "linear":
        return f"{polymer_params['chain_length']}"
    elif polymer_type == "ring":
        return f"{polymer_params['ring_length']}"
    elif polymer_type == "star":
        return f"{polymer_params['arm_length']}_{polymer_params['num_arms']}"
    elif polymer_type == "dendrimer":
        return f"{polymer_params['generations']}_{polymer_params['branching_factor']}"
    return ""


def main():
    """Main function to run polymer workflows with type-specific test configs."""

    # Define test configs for each type (only include relevant params)
    test_configs = {
        "linear": {
            "polymer_type": "linear",
            "box_size": 20.0,
            "solvent_density": 0.3,
            "run_steps": 10000,
            #"thermostat": "langevin",
            "thermostat": "nosehoover",
            "interaction_params": {"pp": 0.3, "ss": 0.3, "sp": 1.5},
            "polymer_params": {
                "chain_length": 50
            }
        },
        "ring": {
            "polymer_type": "ring",
            "box_size": 20.0,
            "solvent_density": 0.3,
            "run_steps": 10000,
            "thermostat": "langevin",
            "interaction_params": {"pp": 0.3, "ss": 0.3, "sp": 1.5}, #
            "polymer_params": {
                "ring_length": 60
            }
        },
        "brush": {
            "polymer_type": "brush",
            "box_size": 20.0,
            "solvent_density": 0.3,
            "run_steps": 10000,
            "thermostat": "langevin",
            "interaction_params": {"pp": 0.3, "ss": 0.3, "sp": 1.5},
            "polymer_params": {
                "backbone_length": 50,
                "grafting_density": 0.2,
                "side_chain_length": 10
            }
        },
        "star": {
            "polymer_type": "star",
            "box_size": 20.0,
            "solvent_density": 0.4,
            "run_steps": 20000,
            "thermostat": "langevin",
            "interaction_params": {"pp": 0.3, "ss": 0.3, "sp": 1.5},
            "polymer_params": {
                "arm_length": 10,
                "num_arms": 4
            }
        },
        "dendrimer": {
            "polymer_type": "dendrimer",
            "box_size": 20.0,
            "solvent_density": 0.2,
            "run_steps": 20000,
            "thermostat": "langevin",
            "interaction_params": {"pp": 0.3, "ss": 0.3, "sp": 1.5},
            "polymer_params": {
                "generations": 3,
                "branching_factor": 3,
                "spacer": 4
            }
        }
    }

    # Run selected workflows (e.g., all for full testing, or just one)
    results = {}
    for polymer_type, config in test_configs.items():
        if polymer_type != "dendrimer":
            continue
        print(f"\nRunning test workflow for {polymer_type.upper()}")
        print("Parameters:")
        for key, value in config.items():
            print(f"  {key}: {value}")
        results[polymer_type] = run_polymer_workflow(config)

    # Optionally, summarize all results (e.g., check statuses)
    print("\nTest Summary:")
    for polymer_type, res in results.items():
        status = res.get('status', 'unknown')
        print(f"  {polymer_type}: {status}")

    return results


if __name__ == "__main__":
    main()
