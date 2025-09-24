from crewai.tools.base_tool import BaseTool
from src.tools.simulation.run_lammps import run_lammps
from src.utils.constants import DEFAULT_INTERACTION_PARAMS, DEFAULT_RUN_STEPS, DEFAULT_THERMOSTAT
from src.utils.logging import setup_logger

logger = setup_logger()


class RunLammpsTool(BaseTool):
    name: str = "RunLammpsSimulation"
    description: str = """Executes LAMMPS simulation using static script.
    Inputs: datafile_path (str), thermostat (str='langevin'/'nosehoover'), interaction_params (dict={'pp':0.5, 'ss':0.5, 'sp':0.5}), run_steps (int=100000).
    Outputs: Dict with 'dump_files' (str), 'final_config' (str), 'log' (str)."""

    def _run(self, datafile_path: str, thermostat: str = DEFAULT_THERMOSTAT, interaction_params: dict = None, run_steps: int = DEFAULT_RUN_STEPS) -> dict:
        if interaction_params is None:
            interaction_params = DEFAULT_INTERACTION_PARAMS
        dump_path = datafile_path.replace(".data", "")
        try:
            results = run_lammps(dump_path, datafile_path, thermostat, interaction_params, run_steps)
            logger.info(f"Simulation completed for {datafile_path}")
            return results
        except Exception as e:
            logger.error(f"Simulation failed: {e}")
            return {"error": str(e)}
