from crewai.tools.base_tool import BaseTool
from src.wrappers.config_wrappers import (
    GenerateLinearPolymerTool,
    GenerateRingPolymerTool,
    GenerateBrushPolymerTool,
    GenerateStarPolymerTool,
    GenerateDendrimerTool
)
from src.wrappers.sim_wrappers import RunLammpsTool
from src.wrappers.analysis_wrappers import ComprehensiveAnalysisTool, PlotAnalysisTool
from src.utils.constants import DEFAULT_INTERACTION_PARAMS
import os
from src.utils.logging import setup_logger

logger = setup_logger()


def run_full_workflow(config):
    """
    Run the complete polymer simulation workflow using a config dict.

    Parameters:
    - config: dict with workflow configuration containing:
        - output_dir: str, the absolute path to the output directory
        - polymer_type: str ('linear'/'ring'/'brush'/'star'/'dendrimer')
        - box_size: float (default 20.0)
        - solvent_density: float (default 0.3)
        - run_steps: int (default 10000)
        - thermostat: str (default 'langevin')
        - interaction_params: dict (default {'pp': 0.5, 'ss': 0.5, 'sp': 0.5})
        - polymer_params: dict with type-specific parameters

    Returns:
    - dict: workflow results with output paths
    """

    # The output directory is now required and passed in the config
    final_output_dir = config.get('output_dir')
    if not final_output_dir:
        raise ValueError("'output_dir' must be provided in the configuration.")

    os.makedirs(final_output_dir, exist_ok=True)
    logger.info(f"Workflow output directory: {final_output_dir}")

    # Extract other config parameters
    polymer_type = config['polymer_type']
    box_size = config.get('box_size', 20.0)
    solvent_density = config.get('solvent_density', 0.3)
    run_steps = config.get('run_steps', 10000)
    thermostat = config.get('thermostat', 'langevin')
    interaction_params = config.get('interaction_params', {'pp': 0.5, 'ss': 0.5, 'sp': 0.5})
    polymer_params = config.get('polymer_params', {})

    # Validate thermostat
    valid_thermostats = ['nosehoover', 'langevin']
    if thermostat not in valid_thermostats:
        raise ValueError(f"Invalid thermostat: {thermostat}. Must be one of {valid_thermostats}")

    # Save simulation parameters to JSON file for record keeping
    import json
    params_file = os.path.join(final_output_dir, 'simulation_parameters.json')
    sim_params = {
        'polymer_type': polymer_type,
        'polymer_params': polymer_params,
        'box_size': box_size,
        'solvent_density': solvent_density,
        'simulation': {
            'run_steps': run_steps,
            'thermostat': thermostat,
            'temperature': 1.0,  # LJ units
            'timestep': 0.01  # LJ units
        },
        'interactions': {
            'polymer_polymer': interaction_params.get('pp', DEFAULT_INTERACTION_PARAMS['pp']),
            'solvent_solvent': interaction_params.get('ss', DEFAULT_INTERACTION_PARAMS['ss']),
            'polymer_solvent': interaction_params.get('sp', DEFAULT_INTERACTION_PARAMS['sp']),
            'description': 'Lennard-Jones interaction parameters (epsilon values in LJ units)'
        }
    }
    with open(params_file, 'w') as f:
        json.dump(sim_params, f, indent=2)
    logger.info(f"Saved simulation parameters to {params_file}")

    results = {
        'config': config,
        'output_dir': final_output_dir,
        'steps': {},
        'output_paths': {},
        'status': 'started'
    }
    try:
        # Validate polymer type
        valid_types = ["linear", "ring", "brush", "star", "dendrimer"]
        if polymer_type not in valid_types:
            raise ValueError(f"Unsupported polymer_type: {polymer_type}. Must be one of {valid_types}")

        # Step 1: Generate polymer configuration and pack solvent
        logger.info("Step 1: Generating polymer configuration and packing solvent")
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
            system_path = polymer_tool._run(
                chain_length=polymer_params["chain_length"],
                solvent_density=solvent_density,
                box_size=box_size,
                output_dir=final_output_dir
            )
        elif polymer_type == "ring":
            system_path = polymer_tool._run(
                chain_length=polymer_params["chain_length"],
                solvent_density=solvent_density,
                box_size=box_size,
                output_dir=final_output_dir
            )
        elif polymer_type == "brush":
            system_path = polymer_tool._run(
                backbone_length=polymer_params["backbone_length"],
                grafting_density=polymer_params["grafting_density"],
                side_chain_length=polymer_params["side_chain_length"],
                solvent_density=solvent_density,
                box_size=box_size,
                output_dir=final_output_dir
            )
        elif polymer_type == "star":
            system_path = polymer_tool._run(
                arm_length=polymer_params["arm_length"],
                num_arms=polymer_params["num_arms"],
                solvent_density=solvent_density,
                box_size=box_size,
                output_dir=final_output_dir
            )
        elif polymer_type == "dendrimer":
            system_path = polymer_tool._run(
                generation=polymer_params["generation"],
                branching_factor=polymer_params["branching_factor"],
                solvent_density=solvent_density,
                box_size=box_size,
                output_dir=final_output_dir
            )

        results["output_paths"]["system"] = system_path
        results["steps"]["configuration_generation"] = "completed"

        # Polymer-only plot is automatically generated by the generation tool
        polymer_plot_path = system_path.replace('.data', '_nosolvent.png')
        results["output_paths"]["polymer_plot"] = polymer_plot_path

        # System plot is automatically generated by the generation tool
        system_plot_path = system_path.replace('.data', '.png')
        results["output_paths"]["system_plot"] = system_plot_path

        # Step 2: Run simulation
        logger.info("Step 2: Running simulation")
        sim_tool = RunLammpsTool()
        sim_results = sim_tool._run(
            datafile_path=system_path,
            thermostat=thermostat,
            interaction_params=interaction_params,
            run_steps=run_steps,
            output_dir=final_output_dir
        )
        results["simulation_result"] = sim_results
        results["steps"]["simulation"] = "completed"

        # Final configuration plots are handled by the analysis workflow
        final_config_path = sim_results["final_config"]
        final_plot_no_solvent_path = os.path.join(final_output_dir, 'final_state_nosolvent.png')
        final_plot_with_solvent_path = os.path.join(final_output_dir, 'final_state.png')
        results["output_paths"]["final_plot_no_solvent"] = final_plot_no_solvent_path
        results["output_paths"]["final_plot_with_solvent"] = final_plot_with_solvent_path

        # Generate final configuration plots
        from src.tools.config_gen.plot_config import plot_configuration
        plot_configuration(final_config_path, plot_solvent=False, output_dir=final_output_dir)
        plot_configuration(final_config_path, plot_solvent=True, output_dir=final_output_dir)

        # Step 3: Perform analysis
        logger.info("Step 3: Performing analysis")
        analysis_tool = ComprehensiveAnalysisTool()
        results_path = analysis_tool._run(
            datafile_path=system_path,
            dump_pattern=sim_results["dump_files"]
        )
        results["output_paths"]["analysis_results"] = results_path
        results["steps"]["analysis"] = "completed"

        # Plot analysis results
        plot_analysis_tool = PlotAnalysisTool()
        analysis_plot_path = plot_analysis_tool._run(results_path)
        results["output_paths"]["analysis_plot"] = analysis_plot_path
        results["output_paths"]["simulation_parameters"] = params_file
        results["status"] = "completed"

    except Exception as e:
        results["status"] = "failed"
        results["error"] = str(e)
        logger.error(f"Workflow failed: {e}")
        raise

    return results


class FullWorkflowTool(BaseTool):
    name: str = "FullWorkflow"
    description: str = """Runs the complete polymer simulation workflow from configuration to analysis.

    This tool accepts individual parameters that define the simulation configuration.

    Required inputs:
    - output_dir (str): The absolute path where all output files will be stored
    - polymer_type (str): Type of polymer ('linear', 'ring', 'brush', 'star', or 'dendrimer')
    - polymer_params_json (str): JSON string with type-specific polymer parameters:
        * For 'linear'/'ring': {"chain_length": int}
        * For 'brush': {"backbone_length": int, "grafting_density": float, "side_chain_length": int}
        * For 'star': {"arm_length": int, "num_arms": int}
        * For 'dendrimer': {"generation": int, "branching_factor": int}

    Optional inputs:
    - box_size (float): Size of simulation box (default 20.0)
    - solvent_density (float): Density of solvent molecules (default 0.5)
    - run_steps (int): Number of simulation steps (default 20000)
    - thermostat (str): Thermostat type ('nosehoover' or 'langevin. default 'langevin')
    - interaction_params (dict): Lennard-Jones interaction parameters (default: {'pp': 0.3, 'ss': 0.3, 'sp': 1.5})
        * "pp": polymer-polymer interaction epsilon
        * "ss": solvent-solvent interaction epsilon
        * "sp" or "ps": solvent-polymer interaction epsilon
        * Good solvent: sp > pp (e.g., {'pp': 0.3, 'ss': 0.3, 'sp': 1.5})
        * Ideal solvent: sp ≈ pp (e.g., {'pp': 0.3, 'ss': 0.3, 'sp': 1.0})
        * Poor solvent: sp < pp (e.g., {'pp': 0.3, 'ss': 0.3, 'sp': 0.5})

    Outputs: A dictionary with workflow results including all output paths and status."""

    def _run(
        self,
        output_dir: str,
        polymer_type: str,
        polymer_params_json: str,
        box_size: float = 20.0,
        solvent_density: float = 0.3,
        run_steps: int = 10000,
        thermostat: str = 'langevin',
        interaction_params_json: str = '{"pp": 0.2, "ss": 0.2, "sp": 1.5}'
    ) -> dict:
        import json
        try:
            # Parse JSON strings
            polymer_params = json.loads(polymer_params_json)
            interaction_params = json.loads(interaction_params_json)

            # Build configuration dictionary
            config = {
                'output_dir': output_dir,
                'polymer_type': polymer_type,
                'polymer_params': polymer_params,
                'box_size': box_size,
                'solvent_density': solvent_density,
                'run_steps': run_steps,
                'thermostat': thermostat,
                'interaction_params': interaction_params
            }

            # Validate thermostat
            valid_thermostats = ['nosehoover', 'langevin']
            if thermostat not in valid_thermostats:
                raise ValueError(f"Invalid thermostat: {thermostat}. Must be one of {valid_thermostats}")

            # Execute workflow
            return run_full_workflow(config)
        except json.JSONDecodeError as e:
            logger.error(f"Failed to parse JSON parameters: {e}")
            raise ValueError(f"Invalid JSON format in parameters: {e}")
        except Exception as e:
            logger.error(f"Full workflow tool failed: {e}")
            raise
