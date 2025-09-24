from crewai.tools.base_tool import BaseTool
from src.tools.analysis.run_analysis import run_complete_analysis, plot_analysis_results, save_analysis_results
import os
from src.utils.logging import setup_logger

logger = setup_logger()


class ComprehensiveAnalysisTool(BaseTool):
    name: str = "ComprehensiveAnalysis"
    description: str = """Runs full analysis on simulation data.
    Inputs: datafile_path (str), dump_pattern (str).
    Outputs: Path to analysis_results.json (str)."""

    def _run(self, datafile_path: str, dump_pattern: str) -> str:
        try:
            results = run_complete_analysis(datafile_path, dump_pattern)
            if results["metadata"]["success"]:
                # Save analysis results in the same directory as the simulation outputs
                # (system_xxx/ folder where final_state.data is located)
                base_name = os.path.splitext(os.path.basename(datafile_path))[0]  # e.g., "system_linear"
                analysis_dir = os.path.join(os.path.dirname(datafile_path), base_name)
                os.makedirs(analysis_dir, exist_ok=True)

                results_path = os.path.join(analysis_dir, "analysis_results.json")
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
    Inputs: analysis_json_path (str).
    Outputs: Path to the generated plot (str)."""

    def _run(self, analysis_json_path: str) -> str:
        try:
            # analysis_json_path should be in system_xxx/analysis_results.json
            # So the plot should also go to system_xxx/
            data_dir = os.path.dirname(analysis_json_path)
            plot_analysis_results(data_dir, analysis_json_path)
            plot_path = os.path.join(data_dir, "conformation_analysis.png")
            logger.info(f"Analysis plot generated at {plot_path}")
            return plot_path
        except Exception as e:
            logger.error(f"Plotting failed: {e}")
            raise
