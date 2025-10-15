from crewai.tools.base_tool import BaseTool
from src.tools.simulation.run_lammps import run_lammps
from src.tools.analysis.run_analysis import run_complete_analysis, plot_analysis_results, save_analysis_results
from src.tools.config_gen.plot_config import plot_configuration
from src.utils.constants import DEFAULT_INTERACTION_PARAMS, DEFAULT_RUN_STEPS, DEFAULT_THERMOSTAT
from src.utils.logging import setup_logger
import os
import json

logger = setup_logger()


class RunLammpsTool(BaseTool):
    name: str = "RunLammpsSimulation"
    description: str = """Executes LAMMPS simulation using static script.

    Inputs:
    - datafile_path (str): Path to the LAMMPS data file
    - thermostat (str): Thermostat type - 'langevin' or 'nosehoover' (default: 'langevin')
    - interaction_params (dict): Lennard-Jones interaction parameters (default: {'pp': 0.2, 'ss': 0.2, 'sp': 1.5})
        * "pp": polymer-polymer interaction epsilon
        * "ss": solvent-solvent interaction epsilon
        * "sp": solvent-polymer interaction epsilon
        * Good solvent: sp > pp (e.g., {'pp': 0.3, 'ss': 0.3, 'sp': 1.5})
        * Ideal solvent: sp ≈ pp (e.g., {'pp': 0.3, 'ss': 0.3, 'sp': 1.0})
        * Poor solvent: sp < pp (e.g., {'pp': 0.3, 'ss': 0.3, 'sp': 0.5})
    - run_steps (int): Number of simulation steps (default: 20000)
    - output_dir (str): Directory to save output files

    Outputs: Dict with keys: dump_files (str), final_config (str), final_polymer_config (str), log (str), simulation_parameters (str)."""

    def _run(self, datafile_path: str, thermostat: str = DEFAULT_THERMOSTAT, interaction_params: dict = None, run_steps: int = DEFAULT_RUN_STEPS, output_dir: str = None) -> dict:

        print("interaction_params", interaction_params)

        if interaction_params is None:
            print("==============================")
            print("using detault interaction params")
            print("==============================")
            interaction_params = DEFAULT_INTERACTION_PARAMS

        if output_dir is None:
            # If no output_dir is provided, use the directory of the datafile_path
            output_dir = os.path.dirname(datafile_path)

        # The dump_path should be inside the specified output_dir
        dump_path = os.path.join(output_dir, os.path.basename(datafile_path).replace(".data", ""))

        try:
            results = run_lammps(dump_path, datafile_path, thermostat, interaction_params, run_steps)

            # Save simulation parameters to a JSON file for record keeping
            params_file = os.path.join(output_dir, 'simulation_parameters.json')

            # Read existing parameters if the file exists (from config phase)
            existing_params = {}
            if os.path.exists(params_file):
                try:
                    with open(params_file, 'r') as f:
                        existing_params = json.load(f)
                except Exception as e:
                    logger.warning(f"Could not read existing parameters: {e}")

            # Update with simulation parameters
            sim_params = existing_params.copy()
            sim_params['simulation'] = {
                'run_steps': run_steps,
                'thermostat': thermostat,
                'temperature': 1.0,  # LJ units
                'timestep': 0.01  # LJ units
            }
            sim_params['interactions'] = {
                'polymer_polymer': interaction_params.get('pp', DEFAULT_INTERACTION_PARAMS['pp']),
                'solvent_solvent': interaction_params.get('ss', DEFAULT_INTERACTION_PARAMS['ss']),
                'polymer_solvent': interaction_params.get('sp', DEFAULT_INTERACTION_PARAMS['sp']),
                'description': 'Lennard-Jones interaction parameters (epsilon values in LJ units)'
            }

            with open(params_file, 'w') as f:
                json.dump(sim_params, f, indent=2)
            logger.info(f"Saved simulation parameters to {params_file}")

            # Add parameters file to results
            results['simulation_parameters'] = params_file

            logger.info(f"Simulation completed for {datafile_path}")
            return results
        except Exception as e:
            logger.error(f"Simulation failed: {e}")
            raise


class SimulationTool(BaseTool):
    name: str = "CompleteSimulationPipeline"
    description: str = """Executes the complete simulation pipeline: LAMMPS simulation, final configuration plotting, trajectory analysis, and analysis result plotting.

    Inputs:
    - datafile_path (str): Path to the LAMMPS data file
    - thermostat (str): Thermostat type - 'langevin' or 'nosehoover' (default: 'langevin')
    - interaction_params (dict): Lennard-Jones interaction parameters (default: {'pp': 0.3, 'ss': 0.3, 'sp': 1.5})
        * "pp": polymer-polymer interaction epsilon
        * "ss": solvent-solvent interaction epsilon
        * "sp" or "ps": solvent-polymer interaction epsilon
        * Good solvent: sp > pp (e.g., {'pp': 0.3, 'ss': 0.3, 'sp': 1.5})
        * Ideal solvent: sp ≈ pp (e.g., {'pp': 0.3, 'ss': 0.3, 'sp': 1.0})
        * Poor solvent: sp < pp (e.g., {'pp': 0.3, 'ss': 0.3, 'sp': 0.5})
    - run_steps (int): Number of simulation steps (default: 20000)
    - output_dir (str): Directory to save all output files

    Outputs: Dict with keys: simulation_results (dict), final_config_plots (str), analysis_results_path (str), analysis_plot_path (str)."""

    def _run(self, datafile_path: str, thermostat: str = DEFAULT_THERMOSTAT, interaction_params: dict = None, run_steps: int = DEFAULT_RUN_STEPS, output_dir: str = None) -> dict:

        if interaction_params is None:
            interaction_params = DEFAULT_INTERACTION_PARAMS

        if output_dir is None:
            output_dir = os.path.dirname(datafile_path)

        # Step 1: Run LAMMPS simulation
        logger.info("Starting LAMMPS simulation...")
        dump_path = os.path.join(output_dir, os.path.basename(datafile_path).replace(".data", ""))
        simulation_results = run_lammps(dump_path, datafile_path, thermostat, interaction_params, run_steps)

        # Save simulation parameters
        params_file = os.path.join(output_dir, 'simulation_parameters.json')
        existing_params = {}
        if os.path.exists(params_file):
            try:
                with open(params_file, 'r') as f:
                    existing_params = json.load(f)
            except Exception as e:
                logger.warning(f"Could not read existing parameters: {e}")

        sim_params = existing_params.copy()
        sim_params['simulation'] = {
            'run_steps': run_steps,
            'thermostat': thermostat,
            'temperature': 1.0,
            'timestep': 0.01
        }
        sim_params['interactions'] = {
            'polymer_polymer': interaction_params.get('pp', DEFAULT_INTERACTION_PARAMS['pp']),
            'solvent_solvent': interaction_params.get('ss', DEFAULT_INTERACTION_PARAMS['ss']),
            'polymer_solvent': interaction_params.get('sp', DEFAULT_INTERACTION_PARAMS['sp']),
            'description': 'Lennard-Jones interaction parameters (epsilon values in LJ units)'
        }

        with open(params_file, 'w') as f:
            json.dump(sim_params, f, indent=2)
        simulation_results['simulation_parameters'] = params_file
        logger.info(f"Simulation completed, parameters saved to {params_file}")

        # Step 2: Plot final configurations
        logger.info("Plotting final configurations...")
        final_config_path = simulation_results.get('final_config')
        if final_config_path:
            plot_configuration(final_config_path, plot_solvent=True, interactive=False, output_dir=output_dir)
            plot_configuration(final_config_path, plot_solvent=False, interactive=False, output_dir=output_dir)
            system_plot = os.path.join(output_dir, "final_state.png")
            polymer_plot = os.path.join(output_dir, "final_state_nosolvent.png")
            final_config_plots = f"Final configuration plots: {system_plot}, {polymer_plot}"
            logger.info(f"Final configuration plots generated: {system_plot}, {polymer_plot}")
        else:
            final_config_plots = "No final configuration available for plotting"
            logger.warning("No final configuration path found in simulation results")

        # Step 3: Run analysis
        logger.info("Running trajectory analysis...")
        dump_pattern = simulation_results.get('dump_files')
        if dump_pattern:
            analysis_results = run_complete_analysis(datafile_path, dump_pattern)
            if analysis_results["metadata"]["success"]:
                analysis_results_path = os.path.join(output_dir, "analysis_results.json")
                save_analysis_results(analysis_results, analysis_results_path)
                logger.info(f"Analysis completed, results saved to {analysis_results_path}")
            else:
                error_msg = analysis_results["metadata"]["error"]
                logger.error(f"Analysis failed: {error_msg}")
                raise Exception(f"Analysis failed: {error_msg}")
        else:
            analysis_results_path = None
            logger.warning("No dump files available for analysis")

        # Step 4: Plot analysis results
        logger.info("Plotting analysis results...")
        if analysis_results_path:
            plot_analysis_results(output_dir, analysis_results_path)
            analysis_plot_path = os.path.join(output_dir, "conformation_analysis.png")
            logger.info(f"Analysis plot generated at {analysis_plot_path}")
        else:
            analysis_plot_path = "No analysis results available for plotting"
            logger.warning("No analysis results to plot")

        # Return comprehensive results
        return {
            'simulation_results': simulation_results,
            'final_config_plots': final_config_plots,
            'analysis_results_path': analysis_results_path,
            'analysis_plot_path': analysis_plot_path
        }
