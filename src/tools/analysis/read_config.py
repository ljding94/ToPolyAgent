import MDAnalysis as mda
import numpy as np
import glob
from collections import defaultdict


def identify_polymer_atoms(atoms, bonds):
    """
    Identify atoms that belong to polymer chains based on bond connectivity.

    Parameters:
    - atoms: Dictionary with atom information
    - bonds: Array of bond tuples (atom1_id, atom2_id, bond_type)

    Returns:
    - polymer_atom_ids: Set of atom IDs that belong to polymers
    - polymer_mask: Boolean mask for polymer atoms
    """
    if len(bonds) == 0:
        # If no bonds, assume all atoms are polymer atoms
        return set(atoms['ids']), np.ones(len(atoms['ids']), dtype=bool)

    # Create adjacency list for bond connectivity
    adjacency = defaultdict(list)
    for bond in bonds:
        # Convert string IDs to integers for consistency
        atom1_id = int(bond[0]) if isinstance(bond[0], str) else bond[0]
        atom2_id = int(bond[1]) if isinstance(bond[1], str) else bond[1]
        adjacency[atom1_id].append(atom2_id)
        adjacency[atom2_id].append(atom1_id)

    # Find all atoms connected by bonds (polymer atoms)
    visited = set()
    polymer_atoms = set()

    def dfs(atom_id):
        """Depth-first search to find connected components"""
        if atom_id in visited:
            return
        visited.add(atom_id)
        polymer_atoms.add(atom_id)
        for neighbor in adjacency[atom_id]:
            dfs(neighbor)

    # Start DFS from each unvisited atom that has bonds
    for bond in bonds:
        atom1_id = int(bond[0]) if isinstance(bond[0], str) else bond[0]
        atom2_id = int(bond[1]) if isinstance(bond[1], str) else bond[1]
        if atom1_id not in visited:
            dfs(atom1_id)
        if atom2_id not in visited:
            dfs(atom2_id)

    # Create boolean mask for polymer atoms
    atom_ids_array = np.array(atoms['ids'])
    polymer_mask = np.isin(atom_ids_array, list(polymer_atoms))

    return polymer_atoms, polymer_mask


def read_lammps_data(data_file):
    """
    Read LAMMPS data file and extract atom and bond information for polymers only.

    Parameters:
    - data_file: Path to LAMMPS data file

    Returns:
    - polymer_atoms: Dictionary with polymer atom information (id, type, position, etc.)
    - bonds: List of bond tuples (atom1_id, atom2_id, bond_type)
    - box: Box dimensions
    """
    u = mda.Universe(data_file, format='DATA')

    # Extract all atom information first
    all_atoms = {
        'ids': u.atoms.ids,
        'types': u.atoms.types,
        'positions': u.atoms.positions,
        'masses': u.atoms.masses if hasattr(u.atoms, 'masses') else None
    }

    # Extract bond information
    bonds = []
    if hasattr(u, 'bonds') and len(u.bonds) > 0:
        for bond in u.bonds:
            bonds.append((bond.atoms[0].id, bond.atoms[1].id, bond.type))
    bonds = np.array(bonds) if bonds else np.array([]).reshape(0, 3)

    # Identify polymer atoms
    polymer_atoms, polymer_mask = identify_polymer_atoms(all_atoms, bonds)

    # Filter atoms to include only polymer atoms
    polymer_atoms = {
        'ids': all_atoms['ids'][polymer_mask],
        'types': all_atoms['types'][polymer_mask],
        'positions': all_atoms['positions'][polymer_mask],
        'masses': all_atoms['masses'][polymer_mask] if all_atoms['masses'] is not None else None
    }

    # Filter bonds to include only those between polymer atoms
    if len(bonds) > 0:
        polymer_atom_ids = set(polymer_atoms['ids'])
        polymer_bonds = []
        for bond in bonds:
            # Convert bond atom IDs to integers for comparison
            atom1_id = int(bond[0]) if isinstance(bond[0], str) else bond[0]
            atom2_id = int(bond[1]) if isinstance(bond[1], str) else bond[1]
            if atom1_id in polymer_atom_ids and atom2_id in polymer_atom_ids:
                polymer_bonds.append(bond)
        bonds = np.array(polymer_bonds)

    # Extract box information
    box = {
        'dimensions': u.dimensions,
        'lengths': u.dimensions[:3],
        'angles': u.dimensions[3:]
    }

    return polymer_atoms, all_atoms, bonds, box


def read_lammps_polymer_trajectory(dump_files, polymer_atom_ids=None):
    """
    Read LAMMPS dump files and extract trajectory information for polymer atoms only.

    Parameters:
    - dump_files: List of dump file paths or glob pattern
    - polymer_atom_ids: Set of polymer atom IDs to include (if None, include all)

    Returns:
    - trajectory: Dictionary with trajectory data
    """
    if isinstance(dump_files, str):
        # If it's a glob pattern, expand it
        if '*' in dump_files:
            dump_files = sorted(glob.glob(dump_files))
        else:
            dump_files = [dump_files]

    if len(dump_files) == 0:
        raise ValueError("No dump files found")

    # Read first file to get topology
    u = mda.Universe(dump_files[0], format='LAMMPSDUMP')

    # Determine which atoms to include
    if polymer_atom_ids is not None:
        # Create mask for polymer atoms
        atom_ids_array = np.array(u.atoms.ids)
        polymer_mask = np.isin(atom_ids_array, list(polymer_atom_ids))
        selected_atoms = u.atoms[polymer_mask]
    else:
        # Include all atoms
        selected_atoms = u.atoms
        polymer_mask = np.ones(len(u.atoms), dtype=bool)

    trajectory = {
        'polymer_atom_positions': [],
        'steps': [],
        'polymer_atom_ids': selected_atoms.ids,
        'polymer_atom_types': selected_atoms.types
    }

    # Read all frames
    step_counter = 0
    for dump_file in dump_files:
        u_temp = mda.Universe(dump_file, format='LAMMPSDUMP')
        if polymer_atom_ids is not None:
            # Apply the same mask to each frame
            atom_ids_array = np.array(u_temp.atoms.ids)
            polymer_mask = np.isin(atom_ids_array, list(polymer_atom_ids))
            selected_atoms_temp = u_temp.atoms[polymer_mask]
        else:
            selected_atoms_temp = u_temp.atoms

        for ts in u_temp.trajectory:
            trajectory['polymer_atom_positions'].append(selected_atoms_temp.positions.copy())
            trajectory['steps'].append(step_counter)
            step_counter += 1

    trajectory['polymer_atom_positions'] = np.array(trajectory['polymer_atom_positions'])
    trajectory['steps'] = np.array(trajectory['steps'])

    return trajectory


def read_simulation_data(data_file, dump_pattern):
    """
    Read both LAMMPS data file and trajectory dump files for polymer atoms only.

    Parameters:
    - data_file: Path to LAMMPS data file
    - dump_pattern: Glob pattern for dump files (e.g., "coord/dump.*.txt")

    Returns:
    - polymer_atoms: Polymer atom information
    - bonds: Polymer bond connectivity
    - trajectory: Polymer trajectory data
    - box: Box information
    """
    # Read topology and identify polymer atoms
    polymer_atoms, all_atoms, bonds, box = read_lammps_data(data_file)

    # Get polymer atom IDs for trajectory filtering
    polymer_atom_ids = set(polymer_atoms['ids'])

    # Read trajectory with polymer atom filtering
    polymer_trajectory = read_lammps_polymer_trajectory(dump_pattern, polymer_atom_ids)

    return polymer_atoms, all_atoms, bonds, polymer_trajectory, box
