#!/usr/bin/env python3
"""
Test script to run wrappers with real data and save outputs to data/test/wrapper/

This script tests the CrewAI wrapper tools with preset configurations for different polymer types.
Similar to test_full_workflow.py but using the wrapper tools instead of direct function calls.

The script defines TEST_CONFIGS with type-specific configurations. Each config includes:
- 'polymer_type': "linear", "ring", "brush", "star", or "dendrimer"
- Common params: 'box_size', 'solvent_density', 'run_steps', 'thermostat', 'interaction_params'
- 'polymer_params': dict of type-specific parameters

Usage:
    python test_wrappers_real.py

To test specific types, modify the selected_types list in main().
Expected runtime: ~2-5 minutes per polymer type depending on parameters.
"""

import os
import sys
from datetime import datetime

# Add project root to path
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(os.path.dirname(current_dir))
sys.path.insert(0, project_root)

from src.wrappers.config_wrappers import (
    GenerateLinearPolymerTool,
    GenerateRingPolymerTool,
    GenerateBrushPolymerTool,
    GenerateStarPolymerTool,
    GenerateDendrimerTool,
    PackSolventTool
)
from src.wrappers.sim_wrappers import RunLammpsTool
from src.wrappers.analysis_wrappers import ComprehensiveAnalysisTool, PlotAnalysisTool


# Preset configurations for different polymer types
TEST_CONFIGS = {
    "linear": {
        "polymer_type": "linear",
        "box_size": 20.0,
        "solvent_density": 0.3,
        "run_steps": 10000,
        "thermostat": "langevin",
        "interaction_params": {"pp": 0.3, "ss": 0.3, "sp": 1.5},
        "polymer_params": {
            "chain_length": 30
        }
    },
    "ring": {
        "polymer_type": "ring",
        "box_size": 20.0,
        "solvent_density": 0.3,
        "run_steps": 10000,
        "thermostat": "langevin",
        "interaction_params": {"pp": 0.3, "ss": 0.3, "sp": 1.5},
        "polymer_params": {
            "chain_length": 30
        }
    },
    "brush": {
        "polymer_type": "brush",
        "box_size": 20.0,
        "solvent_density": 0.2,
        "run_steps": 10000,
        "thermostat": "langevin",
        "interaction_params": {"pp": 0.3, "ss": 0.3, "sp": 1.5},
        "polymer_params": {
            "backbone_length": 20,
            "grafting_density": 0.2,
            "side_chain_length": 5
        }
    },
    "star": {
        "polymer_type": "star",
        "box_size": 20.0,
        "solvent_density": 0.4,
        "run_steps": 15000,
        "thermostat": "langevin",
        "interaction_params": {"pp": 0.3, "ss": 0.3, "sp": 1.5},
        "polymer_params": {
            "arm_length": 8,
            "num_arms": 4
        }
    },
    "dendrimer": {
        "polymer_type": "dendrimer",
        "box_size": 20.0,
        "solvent_density": 0.2,
        "run_steps": 15000,
        "thermostat": "langevin",
        "interaction_params": {"pp": 0.3, "ss": 0.3, "sp": 1.5},
        "polymer_params": {
            "generation": 2,
            "branching_factor": 3
        }
    }
}


def run_wrapper_test(config, output_dir=None):
    """
    Run the complete wrapper test workflow using a config dict.

    Parameters:
    - config: dict with polymer configuration
    - output_dir: str, output directory (auto-generated if None)

    Returns:
    - dict: test results and paths
    """

    polymer_type = config['polymer_type']
    box_size = config.get('box_size', 20.0)
    solvent_density = config.get('solvent_density', 0.3)
    run_steps = config.get('run_steps', 10000)
    thermostat = config.get('thermostat', "langevin")
    interaction_params = config.get('interaction_params', {"pp": 0.5, "ss": 0.5, "sp": 0.5})
    polymer_params = config.get('polymer_params', {})

    print(f"\nTesting {polymer_type.upper()} polymer with wrappers")
    print(f"Parameters: {config}")

    # Create output directory (similar to test_full_workflow.py)
    if output_dir is None:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        # Create a suffix based on key parameters
        if polymer_type == "linear":
            suffix = f"{polymer_params['chain_length']}"
        elif polymer_type == "ring":
            suffix = f"{polymer_params['chain_length']}"
        elif polymer_type == "brush":
            suffix = f"{polymer_params['backbone_length']}_{polymer_params['grafting_density']}_{polymer_params['side_chain_length']}"
        elif polymer_type == "star":
            suffix = f"{polymer_params['arm_length']}_{polymer_params['num_arms']}"
        elif polymer_type == "dendrimer":
            suffix = f"{polymer_params['generation']}_{polymer_params['branching_factor']}"
        else:
            suffix = "unknown"
        output_dir = os.path.join(project_root, 'data', 'test', 'wrapper', f"{polymer_type}_workflow_{suffix}_{timestamp}")

    output_dir = os.path.abspath(output_dir)
    os.makedirs(output_dir, exist_ok=True)
    os.chdir(output_dir)

    print(f"Output directory: {output_dir}")

    test_results = {
        "config": config,
        "steps": {},
        "output_paths": {}
    }

    try:
        # Step 1: Generate polymer configuration
        print("\nStep 1: Generating polymer configuration...")
        tool_map = {
            "linear": GenerateLinearPolymerTool(),
            "ring": GenerateRingPolymerTool(),
            "brush": GenerateBrushPolymerTool(),
            "star": GenerateStarPolymerTool(),
            "dendrimer": GenerateDendrimerTool()
        }
        polymer_tool = tool_map[polymer_type]

        # Call with specific parameters
        if polymer_type == "linear":
            polymer_path = polymer_tool._run(
                chain_length=polymer_params["chain_length"],
                box_size=box_size,
                output_dir=output_dir
            )
        elif polymer_type == "ring":
            polymer_path = polymer_tool._run(
                chain_length=polymer_params["chain_length"],
                box_size=box_size,
                output_dir=output_dir
            )
        elif polymer_type == "brush":
            polymer_path = polymer_tool._run(
                backbone_length=polymer_params["backbone_length"],
                grafting_density=polymer_params["grafting_density"],
                side_chain_length=polymer_params["side_chain_length"],
                box_size=box_size,
                output_dir=output_dir
            )
        elif polymer_type == "star":
            polymer_path = polymer_tool._run(
                arm_length=polymer_params["arm_length"],
                num_arms=polymer_params["num_arms"],
                box_size=box_size,
                output_dir=output_dir
            )
        elif polymer_type == "dendrimer":
            polymer_path = polymer_tool._run(
                generation=polymer_params["generation"],
                branching_factor=polymer_params["branching_factor"],
                box_size=box_size,
                output_dir=output_dir
            )

        print(f"Generated polymer: {polymer_path}")
        test_results["output_paths"]["polymer"] = polymer_path

        # Polymer plot is automatically generated
        polymer_plot_path = polymer_path.replace('.data', '_nosolvent.png')
        print(f"Plotted polymer config: {polymer_plot_path}")
        test_results["output_paths"]["polymer_plot"] = polymer_plot_path

        # Step 2: Pack solvent
        print("\nStep 2: Packing solvent...")
        solvent_tool = PackSolventTool()
        system_path = solvent_tool._run(polymer_path, solvent_density, box_size, output_dir)
        print(f"Packed solvent: {system_path}")
        test_results["output_paths"]["system"] = system_path

        # System plot is automatically generated
        system_plot_path = system_path.replace('.data', '.png')
        print(f"Plotted system config: {system_plot_path}")
        test_results["output_paths"]["system_plot"] = system_plot_path

        # Step 3: Run simulation
        print("\nStep 3: Running simulation...")
        sim_tool = RunLammpsTool()
        sim_result = sim_tool._run(system_path, thermostat, interaction_params, run_steps)
        print(f"Simulation completed: {sim_result}")
        test_results["simulation_result"] = sim_result

        # Final plots are handled by analysis workflow
        final_config = sim_result["final_config"]
        final_plot_no_solvent = final_config.replace('.data', '_nosolvent.png')
        final_plot_with_solvent = final_config.replace('.data', '.png')
        print(f"Plotted final configs: {final_plot_no_solvent}, {final_plot_with_solvent}")
        test_results["output_paths"]["final_plot_no_solvent"] = final_plot_no_solvent
        test_results["output_paths"]["final_plot_with_solvent"] = final_plot_with_solvent

        # Step 4: Run analysis
        print("\nStep 4: Running analysis...")
        dump_pattern = sim_result["dump_files"]
        analysis_tool = ComprehensiveAnalysisTool()
        results_path = analysis_tool._run(system_path, dump_pattern)
        print(f"Analysis results: {results_path}")
        test_results["output_paths"]["analysis_results"] = results_path

        # Plot analysis results
        plot_analysis_tool = PlotAnalysisTool()
        plot_path = plot_analysis_tool._run(results_path)
        print(f"Analysis plot: {plot_path}")
        test_results["output_paths"]["analysis_plot"] = plot_path

        test_results["status"] = "success"

    except Exception as e:
        print(f"Test failed: {e}")
        test_results["status"] = "failed"
        test_results["error"] = str(e)

    return test_results


def main():
    """Main function to run wrapper tests with different polymer types."""

    # Run tests for selected polymer types
    results = {}
    selected_types = ["linear", "brush"]  # Test a few types for speed

    for polymer_type in selected_types:
        if polymer_type in TEST_CONFIGS:
            config = TEST_CONFIGS[polymer_type]
            results[polymer_type] = run_wrapper_test(config)

    # Print summary
    print("\n" + "="*60)
    print("WRAPPER TEST SUMMARY")
    print("="*60)
    for polymer_type, res in results.items():
        status = res.get('status', 'unknown')
        print(f"  {polymer_type}: {status}")
        if status == "failed":
            print(f"    Error: {res.get('error', 'Unknown')}")

    return results


if __name__ == "__main__":
    main()
