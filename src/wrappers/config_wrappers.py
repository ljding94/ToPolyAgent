from crewai.tools.base_tool import BaseTool
from src.tools.config_gen.generate_polymer import (
    generate_linear_polymer_config,
    generate_ring_polymer_config,
    generate_brush_polymer_config,
    generate_star_polymer_config,
    generate_dendrimer_config
)
from src.tools.config_gen.pack_solvent import pack_solvent
from src.tools.config_gen.plot_config import plot_configuration
from src.utils.constants import DEFAULT_BOX_SIZE
from src.utils.logging import setup_logger

logger = setup_logger()


class PolymerGeneratorTool(BaseTool):
    name: str = "GeneratePolymerConfig"
    description: str = """Generates LAMMPS datafile for specified polymer topology.
    Inputs: polymer_type (str: 'linear'/'ring'/'brush'/'star'/'dendrimer'), params (dict: topology-specific, e.g., {'chain_length':30} for linear), box_size (float=50.0).
    Outputs: Path to polymer datafile (str)."""

    def _run(self, polymer_type: str, params: dict, box_size: float = DEFAULT_BOX_SIZE) -> str:
        try:
            if polymer_type == "linear":
                path = generate_linear_polymer_config(params.get("chain_length", 30), box_size)
            elif polymer_type == "ring":
                path = generate_ring_polymer_config(params.get("chain_length", 30), box_size)
            elif polymer_type == "brush":
                path = generate_brush_polymer_config(
                    params.get("backbone_length", 50),
                    params.get("grafting_density", 0.3),
                    params.get("side_chain_length", 10),
                    box_size
                )
            elif polymer_type == "star":
                path = generate_star_polymer_config(
                    params.get("arm_length", 10),
                    params.get("num_arms", 4),
                    box_size
                )
            elif polymer_type == "dendrimer":
                path = generate_dendrimer_config(
                    params.get("generation", 3),
                    params.get("branching_factor", 3),
                    box_size
                )
            else:
                raise ValueError(f"Invalid polymer_type: {polymer_type}")
            logger.info(f"Generated {polymer_type} config at {path}")
            return path
        except Exception as e:
            logger.error(f"Failed to generate polymer config: {e}")
            raise


class PackSolventTool(BaseTool):
    name: str = "PackSolvent"
    description: str = """Packs solvent into the polymer system.
    Inputs: polymer_file (str: path to polymer datafile), solvent_density (float: 0-1), box_size (float=50.0).
    Outputs: Path to system datafile (str)."""

    def _run(self, polymer_file: str, solvent_density: float, box_size: float = DEFAULT_BOX_SIZE) -> str:
        try:
            path = pack_solvent(polymer_file, solvent_density, box_size)
            logger.info(f"Packed solvent into {polymer_file}, output at {path}")
            return path
        except Exception as e:
            logger.error(f"Failed to pack solvent: {e}")
            raise


class PlotConfigTool(BaseTool):
    name: str = "PlotConfig"
    description: str = """Plots the configuration of the system.
    Inputs: datafile_path (str: path to datafile), show_solvent (bool=True).
    Outputs: Path to plot image (str)."""

    def _run(self, datafile_path: str, show_solvent: bool = True) -> str:
        try:
            path = plot_configuration(datafile_path, plot_solvent=show_solvent)
            logger.info(f"Plotted config for {datafile_path} at {path}")
            return path
        except Exception as e:
            logger.error(f"Failed to plot config: {e}")
            raise
