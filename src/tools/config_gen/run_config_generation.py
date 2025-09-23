import os
import sys

# Add the src directory to the path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..'))

from src.tools.config_gen.generate_polymer import (
    generate_linear_polymer_config,
    generate_ring_polymer_config,
    generate_brush_polymer_config,
    generate_star_polymer_config,
    generate_dendrimer_config
)
from src.tools.config_gen.pack_solvent import pack_solvent
from src.tools.config_gen.plot_config import plot_configuration


def run_config_generation(polymer_type, output_dir=None, plot=True, interactive=False, **kwargs):
    """
    Generate polymer configuration, pack solvent, and optionally plot the result.

    Parameters:
    - polymer_type: str, type of polymer ('linear', 'ring', 'brush', 'star', 'dendrimer')
    - output_dir: str, output directory (current dir if None)
    - plot: bool, whether to generate plots
    - interactive: bool, whether to show interactive plot (default: False)
    - **kwargs: additional parameters depending on polymer type

    Returns:
    - dict: paths to generated files
    """

    if output_dir is None:
        output_dir = os.getcwd()
    else:
        os.makedirs(output_dir, exist_ok=True)

    results = {}

    # Extract common parameters
    solvent_density = kwargs.pop('solvent_density', 0.1)
    box_size = kwargs.pop('box_size', 50.0)

    # Generate polymer configuration
    print(f"Generating {polymer_type} polymer configuration...")

    if polymer_type == 'linear':
        polymer_file = generate_linear_polymer_config(**kwargs)
    elif polymer_type == 'ring':
        polymer_file = generate_ring_polymer_config(**kwargs)
    elif polymer_type == 'brush':
        polymer_file = generate_brush_polymer_config(**kwargs)
    elif polymer_type == 'star':
        polymer_file = generate_star_polymer_config(**kwargs)
    elif polymer_type == 'dendrimer':
        polymer_file = generate_dendrimer_config(**kwargs)
    else:
        raise ValueError(f"Unknown polymer type: {polymer_type}")

    results['polymer_config'] = polymer_file
    print(f"Polymer configuration saved to: {polymer_file}")

    if plot:
        plot_configuration(polymer_file, interactive=interactive)
        polymer_plot = os.path.splitext(polymer_file)[0] + '.png'
        results['polymer_plot'] = polymer_plot
        print(f"Polymer plot saved to: {polymer_plot}")

    # Pack solvent
    print("Packing solvent...")

    system_file = pack_solvent(
        polymer_datafile=polymer_file,
        solvent_density=solvent_density,
        box_size=box_size
    )

    results['system_config'] = system_file
    print(f"System configuration saved to: {system_file}")

    if plot:
        plot_configuration(system_file, interactive=interactive)
        system_plot = os.path.splitext(system_file)[0] + '.png'
        results['system_plot'] = system_plot
        print(f"System plot saved to: {system_plot}")

    return results


def generate_linear_config(chain_length, box_size=50.0, solvent_density=0.1, output_dir=None, plot=True, interactive=False):
    """
    Generate linear polymer configuration with solvent and plot.

    Parameters:
    - chain_length: int, number of beads in the chain
    - box_size: float, simulation box size
    - solvent_density: float, number density of solvent beads
    - output_dir: str, output directory
    - plot: bool, whether to generate plots

    Returns:
    - dict: paths to generated files
    """
    return run_config_generation(
        'linear',
        output_dir=output_dir,
        plot=plot,
        interactive=interactive,
        chain_length=chain_length,
        box_size=box_size,
        solvent_density=solvent_density
    )


def generate_ring_config(ring_length, box_size=50.0, solvent_density=0.1, output_dir=None, plot=True, interactive=False):
    """
    Generate ring polymer configuration with solvent and plot.

    Parameters:
    - ring_length: int, number of beads in the ring
    - box_size: float, simulation box size
    - solvent_density: float, number density of solvent beads
    - output_dir: str, output directory
    - plot: bool, whether to generate plots

    Returns:
    - dict: paths to generated files
    """
    return run_config_generation(
        'ring',
        output_dir=output_dir,
        plot=plot,
        interactive=interactive,
        ring_length=ring_length,
        box_size=box_size,
        solvent_density=solvent_density
    )


def generate_brush_config(backbone_length, grafting_density, side_chain_length, box_size=50.0, solvent_density=0.1, output_dir=None, plot=True, interactive=False):
    """
    Generate brush polymer configuration with solvent and plot.

    Parameters:
    - backbone_length: int, number of beads in backbone
    - grafting_density: float, probability of grafting side chains
    - side_chain_length: int, number of beads in each side chain
    - box_size: float, simulation box size
    - solvent_density: float, number density of solvent beads
    - output_dir: str, output directory
    - plot: bool, whether to generate plots

    Returns:
    - dict: paths to generated files
    """
    return run_config_generation(
        'brush',
        output_dir=output_dir,
        plot=plot,
        interactive=interactive,
        backbone_length=backbone_length,
        grafting_density=grafting_density,
        side_chain_length=side_chain_length,
        box_size=box_size,
        solvent_density=solvent_density
    )


def generate_star_config(arm_length, num_arms, box_size=50.0, solvent_density=0.1, output_dir=None, plot=True, interactive=False):
    """
    Generate star polymer configuration with solvent and plot.

    Parameters:
    - arm_length: int, number of beads in each arm
    - num_arms: int, number of arms
    - box_size: float, simulation box size
    - solvent_density: float, number density of solvent beads
    - output_dir: str, output directory
    - plot: bool, whether to generate plots

    Returns:
    - dict: paths to generated files
    """
    return run_config_generation(
        'star',
        output_dir=output_dir,
        plot=plot,
        interactive=interactive,
        arm_length=arm_length,
        num_arms=num_arms,
        box_size=box_size,
        solvent_density=solvent_density
    )


def run_dendrimer_config(generations, branching_factor, box_size=50.0, solvent_density=0.1, output_dir=None, plot=True, interactive=False):
    """
    Generate dendrimer configuration with solvent and plot.

    Parameters:
    - generations: int, number of generations
    - branching_factor: int, branching factor
    - box_size: float, simulation box size
    - solvent_density: float, number density of solvent beads
    - output_dir: str, output directory
    - plot: bool, whether to generate plots

    Returns:
    - dict: paths to generated files
    """
    return run_config_generation(
        'dendrimer',
        output_dir=output_dir,
        plot=plot,
        interactive=interactive,
        generations=generations,
        branching_factor=branching_factor,
        box_size=box_size,
        solvent_density=solvent_density
    )


if __name__ == "__main__":
    # Example usage - generate a brush polymer configuration
    print("Running example configuration generation...")

    # Create output directory
    output_dir = "example_config_output"
    os.makedirs(output_dir, exist_ok=True)

    # Generate brush polymer configuration
    results = generate_brush_config(
        backbone_length=50,
        grafting_density=0.2,
        side_chain_length=10,
        box_size=30.0,
        solvent_density=0.1,
        output_dir=output_dir,
        plot=True,
        interactive=True
    )

    print("\nConfiguration generation completed!")
    print("Generated files:")
    for key, path in results.items():
        print(f"  {key}: {path}")

    print(f"\nOutput directory: {os.path.abspath(output_dir)}")
    print("Check the directory for .data files and .png plots.")
