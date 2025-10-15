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
import json
import os

logger = setup_logger()


def save_polymer_params(output_dir: str, polymer_type: str, polymer_params: dict, box_size: float, solvent_density: float = None):
    """Helper function to save polymer configuration parameters to JSON file."""
    params_file = os.path.join(output_dir, 'simulation_parameters.json')

    # Read existing parameters if the file exists
    existing_params = {}
    if os.path.exists(params_file):
        try:
            with open(params_file, 'r') as f:
                existing_params = json.load(f)
        except Exception as e:
            logger.warning(f"Could not read existing parameters: {e}")

    # Update with polymer configuration parameters
    sim_params = existing_params.copy()
    sim_params['polymer_type'] = polymer_type
    sim_params['polymer_params'] = polymer_params
    sim_params['box_size'] = box_size
    if solvent_density is not None:
        sim_params['solvent_density'] = solvent_density

    with open(params_file, 'w') as f:
        json.dump(sim_params, f, indent=2)
    logger.info(f"Saved polymer configuration parameters to {params_file}")


class GenerateLinearPolymerTool(BaseTool):
    name: str = "GenerateLinearPolymer"
    description: str = """Generates a LAMMPS datafile for a linear polymer and packs it with solvent.
    Inputs: chain_length (int: >=2, default 20), solvent_density (float: 0-1, default 0.3), box_size (float=20.0), output_dir (str: directory to save files, required).
    Outputs: Path to the system datafile with solvent (str)."""

    def _run(self, chain_length: int, solvent_density: float, box_size: float = DEFAULT_BOX_SIZE, output_dir: str = None) -> str:
        if chain_length < 2:
            raise ValueError("chain_length must be at least 2")
        if not (0 <= solvent_density <= 1):
            raise ValueError("solvent_density must be between 0 and 1")
        if output_dir is None:
            raise ValueError("output_dir is required")
        try:
            # Generate polymer config
            polymer_path = generate_linear_polymer_config(chain_length, box_size, output_dir)

            # Plot polymer-only config
            plot_configuration(polymer_path, plot_solvent=False, output_dir=output_dir)

            # Pack solvent
            system_path = pack_solvent(polymer_path, solvent_density, box_size, output_dir)

            # Plot system with solvent
            plot_configuration(system_path, plot_solvent=True, output_dir=output_dir)

            # Save all parameters
            save_polymer_params(
                output_dir,
                'linear',
                {'chain_length': chain_length},
                box_size,
                solvent_density
            )

            logger.info(f"Generated linear polymer system with solvent at {system_path}")
            return system_path
        except Exception as e:
            logger.error(f"Failed to generate linear polymer system: {e}")
            raise


class GenerateRingPolymerTool(BaseTool):
    name: str = "GenerateRingPolymer"
    description: str = """Generates a LAMMPS datafile for a ring polymer and packs it with solvent.
    Inputs: chain_length (int: number of beads in ring, >=3, default 30), solvent_density (float: 0-1, default 0.3), box_size (float=20.0), output_dir (str: directory to save files, required).
    Outputs: Path to the system datafile with solvent (str)."""

    def _run(self, chain_length: int, solvent_density: float, box_size: float = DEFAULT_BOX_SIZE, output_dir: str = None) -> str:
        if chain_length < 3:
            raise ValueError("chain_length must be at least 3 for a ring polymer")
        if not (0 <= solvent_density <= 1):
            raise ValueError("solvent_density must be between 0 and 1")
        if output_dir is None:
            raise ValueError("output_dir is required")
        try:
            # Generate polymer config
            polymer_path = generate_ring_polymer_config(chain_length, box_size, output_dir)

            # Plot polymer-only config
            plot_configuration(polymer_path, plot_solvent=False, output_dir=output_dir)

            # Pack solvent
            system_path = pack_solvent(polymer_path, solvent_density, box_size, output_dir)

            # Plot system with solvent
            plot_configuration(system_path, plot_solvent=True, output_dir=output_dir)

            # Save all parameters
            save_polymer_params(
                output_dir,
                'ring',
                {'chain_length': chain_length},
                box_size,
                solvent_density
            )

            logger.info(f"Generated ring polymer system with solvent at {system_path}")
            return system_path
        except Exception as e:
            logger.error(f"Failed to generate ring polymer system: {e}")
            raise


class GenerateBrushPolymerTool(BaseTool):
    name: str = "GenerateBrushPolymer"
    description: str = """Generates a LAMMPS datafile for a brush polymer and packs it with solvent.
    Inputs: backbone_length (int: >=5, default 20), grafting_density (float: 0-1, default 0.3), side_chain_length (int: >=2, default 5), solvent_density (float: 0-1, default 0.3), box_size (float=20.0), output_dir (str: required).
    Outputs: Path to the system datafile with solvent (str)."""

    def _run(self, backbone_length: int, grafting_density: float, side_chain_length: int, solvent_density: float, box_size: float = DEFAULT_BOX_SIZE, output_dir: str = None) -> str:
        if backbone_length < 5:
            raise ValueError("backbone_length must be at least 5")
        if not (0 <= grafting_density <= 1):
            raise ValueError("grafting_density must be between 0 and 1")
        if side_chain_length < 2:
            raise ValueError("side_chain_length must be at least 2")
        if not (0 <= solvent_density <= 1):
            raise ValueError("solvent_density must be between 0 and 1")
        if output_dir is None:
            raise ValueError("output_dir is required")
        try:
            # Generate polymer config
            polymer_path = generate_brush_polymer_config(backbone_length, grafting_density, side_chain_length, box_size, output_dir)

            # Plot polymer-only config
            plot_configuration(polymer_path, plot_solvent=False, output_dir=output_dir)

            # Pack solvent
            system_path = pack_solvent(polymer_path, solvent_density, box_size, output_dir)

            # Plot system with solvent
            plot_configuration(system_path, plot_solvent=True, output_dir=output_dir)

            # Save all parameters
            save_polymer_params(
                output_dir,
                'brush',
                {'backbone_length': backbone_length, 'grafting_density': grafting_density, 'side_chain_length': side_chain_length},
                box_size,
                solvent_density
            )

            logger.info(f"Generated brush polymer system with solvent at {system_path}")
            return system_path
        except Exception as e:
            logger.error(f"Failed to generate brush polymer system: {e}")
            raise


class GenerateStarPolymerTool(BaseTool):
    name: str = "GenerateStarPolymer"
    description: str = """Generates a LAMMPS datafile for a star polymer and packs it with solvent.
    Inputs: arm_length (int: >=5, default 10), num_arms (int: >=3, default 6), solvent_density (float: 0-1, default 0.3), box_size (float=20.0), output_dir (str: required).
    Outputs: Path to the system datafile with solvent (str)."""

    def _run(self, arm_length: int, num_arms: int, solvent_density: float, box_size: float = DEFAULT_BOX_SIZE, output_dir: str = None) -> str:
        if arm_length < 5:
            raise ValueError("arm_length must be at least 5")
        if num_arms < 3:
            raise ValueError("num_arms must be at least 3")
        if not (0 <= solvent_density <= 1):
            raise ValueError("solvent_density must be between 0 and 1")
        if output_dir is None:
            raise ValueError("output_dir is required")
        try:
            # Generate polymer config
            polymer_path = generate_star_polymer_config(arm_length, num_arms, box_size, output_dir)

            # Plot polymer-only config
            plot_configuration(polymer_path, plot_solvent=False, output_dir=output_dir)

            # Pack solvent
            system_path = pack_solvent(polymer_path, solvent_density, box_size, output_dir)

            # Plot system with solvent
            plot_configuration(system_path, plot_solvent=True, output_dir=output_dir)

            # Save all parameters
            save_polymer_params(
                output_dir,
                'star',
                {'arm_length': arm_length, 'num_arms': num_arms},
                box_size,
                solvent_density
            )

            logger.info(f"Generated star polymer system with solvent at {system_path}")
            return system_path
        except Exception as e:
            logger.error(f"Failed to generate star polymer system: {e}")
            raise


class GenerateDendrimerTool(BaseTool):
    name: str = "GenerateDendrimer"
    description: str = """Generates a LAMMPS datafile for a dendrimer and packs it with solvent.
    Inputs: generation (int: >=1, default 3), branching_factor (int: >=2, default 2), solvent_density (float: 0-1, default 0.3), box_size (float=20.0), output_dir (str: required).
    Outputs: Path to the system datafile with solvent (str)."""

    def _run(self, generation: int, branching_factor: int, solvent_density: float, box_size: float = DEFAULT_BOX_SIZE, output_dir: str = None) -> str:
        if generation < 1:
            raise ValueError("generation must be at least 1")
        if branching_factor < 2:
            raise ValueError("branching_factor must be at least 2")
        if not (0 <= solvent_density <= 1):
            raise ValueError("solvent_density must be between 0 and 1")
        if output_dir is None:
            raise ValueError("output_dir is required")
        try:
            # Generate polymer config
            polymer_path = generate_dendrimer_config(generation, branching_factor, 5, box_size, output_dir)

            # Plot polymer-only config
            plot_configuration(polymer_path, plot_solvent=False, output_dir=output_dir)

            # Pack solvent
            system_path = pack_solvent(polymer_path, solvent_density, box_size, output_dir)

            # Plot system with solvent
            plot_configuration(system_path, plot_solvent=True, output_dir=output_dir)

            # Save all parameters
            save_polymer_params(
                output_dir,
                'dendrimer',
                {'generation': generation, 'branching_factor': branching_factor},
                box_size,
                solvent_density
            )

            logger.info(f"Generated dendrimer system with solvent at {system_path}")
            return system_path
        except Exception as e:
            logger.error(f"Failed to generate dendrimer system: {e}")
            raise


class PackSolventTool(BaseTool):
    name: str = "PackSolvent"
    description: str = """Packs solvent into the polymer system. This tool is deprecated and should not be used directly. The polymer generation tools now handle solvent packing.
    Inputs: polymer_file (str: path to polymer datafile), solvent_density (float: 0-1, default 0.5), box_size (float=20.0), output_dir (str: directory to save files).
    Outputs: Path to system datafile (str)."""

    def _run(self, polymer_file: str, solvent_density: float, box_size: float = DEFAULT_BOX_SIZE, output_dir: str = None) -> str:
        logger.warning("PackSolventTool is deprecated. Solvent packing is now integrated into polymer generation tools.")
        try:
            path = pack_solvent(polymer_file, solvent_density, box_size, output_dir)
            # Update the parameters file with solvent density
            if output_dir:
                params_file = os.path.join(output_dir, 'simulation_parameters.json')
                if os.path.exists(params_file):
                    try:
                        with open(params_file, 'r') as f:
                            sim_params = json.load(f)
                        sim_params['solvent_density'] = solvent_density
                        with open(params_file, 'w') as f:
                            json.dump(sim_params, f, indent=2)
                    except Exception as e:
                        logger.warning(f"Could not update solvent density: {e}")
            # Automatically plot the system configuration with solvent
            plot_path = plot_configuration(path, plot_solvent=True, output_dir=output_dir)
            logger.info(f"Packed solvent into {polymer_file}, system at {path} and plotted at {plot_path}")
            return path
        except Exception as e:
            logger.error(f"Failed to pack solvent: {e}")
            raise
