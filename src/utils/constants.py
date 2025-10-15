# Default constants for ToPolyAgent
DEFAULT_BOX_SIZE = 20.0
DEFAULT_BOND_LENGTH = 1.0
DEFAULT_INTERACTION_PARAMS = {
    "pp": 0.2,  # polymer-polymer
    "ss": 0.2,  # solvent-solvent
    "sp": 1.5   # solvent-polymer
}
DEFAULT_RUN_STEPS = 20000
DEFAULT_THERMOSTAT = "langevin"

# Default polymer parameters
DEFAULT_POLYMER_PARAMS = {
    "linear": {
        "chain_length": 20
    },
    "ring": {
        "chain_length": 30
    },
    "brush": {
        "backbone_length": 20,
        "grafting_density": 0.3,
        "side_chain_length": 5
    },
    "star": {
        "arm_length": 15,
        "num_arms": 6
    },
    "dendrimer": {
        "generation": 3,
        "branching_factor": 3,
        "spacer": 5
    }
}
