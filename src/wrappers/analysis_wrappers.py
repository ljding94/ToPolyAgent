from crewai.tools.base_tool import BaseTool
from src.tools.analysis.run_analysis import run_complete_analysis, plot_analysis_results, save_analysis_results
import os
from src.utils.logging import setup_logger
from src.tools.config_gen.plot_config import plot_configuration

logger = setup_logger()


class ComprehensiveAnalysisTool(BaseTool):
    name: str = "ComprehensiveAnalysis"
    description: str = """Runs full analysis on simulation data.
    Inputs: datafile_path (str), dump_pattern (str), output_dir (str, optional).
    Outputs: Path to analysis_results.json (str)."""

    def _run(self, datafile_path: str, dump_pattern: str, output_dir: str = None) -> str:
        try:
            results = run_complete_analysis(datafile_path, dump_pattern)
            if results["metadata"]["success"]:
                # Save analysis results to specified output_dir or same directory as datafile
                if output_dir is None:
                    output_dir = os.path.dirname(datafile_path)
                results_path = os.path.join(output_dir, "analysis_results.json")
                save_analysis_results(results, results_path)
                logger.info(f"Analysis completed, results at {results_path}")
                return results_path
            else:
                error_msg = results["metadata"]["error"]
                logger.error(f"Analysis failed: {error_msg}")
                raise Exception(f"Analysis failed: {error_msg}")
        except Exception as e:
            logger.error(f"Analysis failed: {e}")
            raise


class PlotAnalysisTool(BaseTool):
    name: str = "PlotAnalysis"
    description: str = """Plots analysis results from a JSON file.
    Inputs: analysis_json_path (str), output_dir (str, optional).
    Outputs: Path to the generated plot (str)."""

    def _run(self, analysis_json_path: str, output_dir: str = None) -> str:
        try:
            # Save plot to specified output_dir or same directory as analysis file
            if output_dir is None:
                output_dir = os.path.dirname(analysis_json_path)
            plot_analysis_results(output_dir, analysis_json_path)
            plot_path = os.path.join(output_dir, "conformation_analysis.png")
            logger.info(f"Analysis plot generated at {plot_path}")
            return plot_path
        except Exception as e:
            logger.error(f"Plotting failed: {e}")
            raise


class PlotFinalConfigurationsTool(BaseTool):
    name: str = "PlotFinalConfigurations"
    description: str = """Plots final configurations from simulation (both with and without solvent).
    Inputs: final_config_path (str), output_dir (str, optional).
    Outputs: Paths to the generated plots (str)."""

    def _run(self, final_config_path: str, output_dir: str = None) -> str:
        try:
            # Determine where to save plots
            if output_dir is None:
                output_dir = os.path.dirname(final_config_path)

            # Derive polymer config path from system config path
            polymer_config_path = final_config_path  # Use the same file, but plot without solvent

            # Plot full system final configuration
            plot_configuration(final_config_path, plot_solvent=True, interactive=False, output_dir=output_dir)

            # Plot polymer-only final configuration (using same data file but without solvent)
            plot_configuration(polymer_config_path, plot_solvent=False, interactive=False, output_dir=output_dir)

            # Get plot paths
            system_plot = os.path.join(output_dir, "final_state.png")
            polymer_plot = os.path.join(output_dir, "final_state_nosolvent.png")

            logger.info(f"Final configuration plots generated: {system_plot}, {polymer_plot}")
            return f"Final configuration plots: {system_plot}, {polymer_plot}"
        except Exception as e:
            logger.error(f"Final configuration plotting failed: {e}")
            raise
