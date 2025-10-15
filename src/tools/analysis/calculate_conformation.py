import numpy as np
from scipy.stats import linregress
from collections import defaultdict


def calculate_radius_of_gyration(positions):
    """
    Calculate the radius of gyration for given atom positions.

    Parameters:
    - positions: Array of atom positions (n_frames, n_atoms, 3) or (n_atoms, 3)

    Returns:
    - rg2: Radius of gyration value(s) squared
    """
    if positions.ndim == 2:
        # Single frame
        positions = positions[np.newaxis, :, :]

    rg2_values = []
    for frame_positions in positions:
        # Calculate center of mass
        com = np.mean(frame_positions, axis=0)
        # Calculate radius of gyration
        rg2_value = np.sum((frame_positions - com)**2) / len(frame_positions)
        rg2_values.append(rg2_value)

    return np.array(rg2_values)


def calculate_persistence_length(positions, bonds, atom_ids, polymer_type="linear", atom_types=None):
    """
    Calculate the persistence length from atom positions and bond connectivity.

    Parameters:
    - positions: Atom positions (n_atoms, 3)
    - bonds: Bond connectivity array (n_bonds, 3) with [atom1_id, atom2_id, bond_type]
    - atom_ids: Array of atom IDs corresponding to positions
    - polymer_type: Type of polymer ("linear", "ring", "brush", "star", "dendrimer")
    - atom_types: Array of atom types corresponding to atom_ids

    Returns:
    - lp: Persistence length (float for linear/ring, dict for others, None for dendrimer)
    """
    if polymer_type == "dendrimer":
        # Persistence length not well-defined for dendrimers due to branching
        return None

    if len(bonds) == 0:
        raise ValueError("No bonds found. Cannot calculate persistence length.")

    # Create mapping from atom ID to position index
    id_to_index = {atom_id: i for i, atom_id in enumerate(atom_ids)}

    if polymer_type in ["linear", "ring"]:
        # Original logic for single chain
        bond_vectors = []
        for bond in bonds:
            if len(bond) >= 2:
                atom1_id, atom2_id = int(bond[0]), int(bond[1])
                if atom1_id in id_to_index and atom2_id in id_to_index:
                    pos1 = positions[id_to_index[atom1_id]]
                    pos2 = positions[id_to_index[atom2_id]]
                    vec = pos2 - pos1
                    norm = np.linalg.norm(vec)
                    if norm > 0:
                        bond_vectors.append(vec / norm)

        if len(bond_vectors) < 2:
            raise ValueError(f"Not enough bonds to calculate persistence length. Found {len(bond_vectors)} bond vectors from {len(bonds)} bonds.")

        # Calculate cos(theta) for consecutive bonds
        cos_thetas = []
        for i in range(len(bond_vectors) - 1):
            cos_theta = np.dot(bond_vectors[i], bond_vectors[i+1])
            cos_thetas.append(cos_theta)

        # Clip individual cosines to [-1, 1] to avoid numerical issues
        cos_thetas = np.clip(cos_thetas, -1, 1)
        if np.any(np.isnan(cos_thetas)):
            lp = None
        else:
            avg_cos = np.mean(cos_thetas)
            if np.isnan(avg_cos) or avg_cos >= 1 or avg_cos <= -1:
                lp = np.inf
            else:
                lp = -1 / np.log(avg_cos)

        return lp
    else:
        # For brush and star polymers, calculate for each chain
        if atom_types is None:
            raise ValueError("atom_types must be provided for brush and star polymers")

        chains = get_chains(bonds, atom_ids, polymer_type, atom_types)
        lps = []
        for chain in chains:
            if len(chain) < 2:
                continue
            chain_positions = np.array([positions[id_to_index[aid]] for aid in chain])
            bond_vectors = []
            for i in range(len(chain)-1):
                pos1 = chain_positions[i]
                pos2 = chain_positions[i+1]
                vec = pos2 - pos1
                norm = np.linalg.norm(vec)
                if norm > 0:
                    bond_vectors.append(vec / norm)

            if len(bond_vectors) < 2:
                lp = np.inf if len(bond_vectors) == 1 else None
            else:
                cos_thetas = [np.dot(bond_vectors[i], bond_vectors[i+1]) for i in range(len(bond_vectors)-1)]
                # Clip individual cosines to [-1, 1] to avoid numerical issues
                cos_thetas = np.clip(cos_thetas, -1, 1)
                if np.any(np.isnan(cos_thetas)):
                    lp = None
                else:
                    avg_cos = np.mean(cos_thetas)
                    if np.isnan(avg_cos) or avg_cos >= 1 or avg_cos <= -1:
                        lp = np.inf
                    else:
                        lp = -1 / np.log(avg_cos)
            lps.append(lp)

        valid_lps = [lp for lp in lps if lp is not None and lp != np.inf]
        if valid_lps:
            overall_lp = np.mean(valid_lps)
        else:
            overall_lp = None

        return {"overall": overall_lp, "chains": lps}


def get_chains(bonds, atom_ids, polymer_type, atom_types):
    """
    Extract chains from bonds based on polymer type.

    Parameters:
    - bonds: Bond connectivity
    - atom_ids: List of atom IDs
    - polymer_type: "brush" or "star"
    - atom_types: List of atom types

    Returns:
    - chains: List of lists, each sublist is atom IDs in chain order
    """
    id_to_index = {atom_id: i for i, atom_id in enumerate(atom_ids)}
    id_to_type = {atom_id: atom_types[i] for i, atom_id in enumerate(atom_ids)}

    # Build adjacency
    adj = defaultdict(list)
    for bond in bonds:
        a1, a2 = int(bond[0]), int(bond[1])
        if a1 in id_to_index and a2 in id_to_index:
            adj[a1].append(a2)
            adj[a2].append(a1)

    chains = []
    if polymer_type == "brush":
        # Backbone: connected type 1 atoms
        backbone_atoms = [aid for aid in atom_ids if id_to_type[aid] == '1']
        if backbone_atoms:
            chain = sort_chain(backbone_atoms, adj)
            chains.append(chain)

    elif polymer_type == "star":
        # One arm: first connected component of type 2
        arm_atoms = [aid for aid in atom_ids if id_to_type[aid] == '2']
        if arm_atoms:
            visited = set()
            chain = dfs_chain(arm_atoms[0], adj, visited, id_to_type, '2')
            if len(chain) > 1:
                sorted_chain = sort_chain(chain, adj)
                chains.append(sorted_chain)

    else:
        # Fallback to single chain
        chain = sort_chain(list(atom_ids), adj)
        chains.append(chain)

    return chains


def dfs_chain(start, adj, visited, id_to_type, allowed_type):
    """DFS to get connected component of given type."""
    chain = []
    stack = [start]
    while stack:
        current = stack.pop()
        if current in visited:
            continue
        visited.add(current)
        if id_to_type[current] == allowed_type:
            chain.append(current)
            for neigh in adj[current]:
                if neigh not in visited and id_to_type[neigh] == allowed_type:
                    stack.append(neigh)
    return chain


def sort_chain(atoms, adj):
    """Sort atoms in chain order."""
    if not atoms:
        return []
    # Find end atoms (degree 1)
    degrees = {a: len([n for n in adj[a] if n in atoms]) for a in atoms}
    ends = [a for a in atoms if degrees[a] == 1]
    if not ends:
        start = atoms[0]
    else:
        start = ends[0]

    # Traverse
    chain = []
    visited = set()
    current = start
    prev = None
    while current not in visited:
        visited.add(current)
        chain.append(current)
        neighbors = [n for n in adj[current] if n != prev and n in atoms]
        if not neighbors:
            break
        prev = current
        current = neighbors[0]
    return chain


def calculate_end_to_end_distance(positions, atom_ids, bonds, polymer_type="linear", atom_types=None):
    """
    Calculate the end-to-end distance for given polymer type.

    Parameters:
    - positions: Atom positions (n_atoms, 3)
    - atom_ids: Array of atom IDs corresponding to positions
    - bonds: Bond connectivity array (n_bonds, 3) with [atom1_id, atom2_id, bond_type]
    - polymer_type: Type of polymer ("linear", "ring", "brush", "star", "dendrimer")
    - atom_types: Array of atom types corresponding to atom_ids (required for brush and star)

    Returns:
    - distance: End-to-end distance (float for linear, dict for brush/star, None for ring/dendrimer)
    """
    if polymer_type in ["ring", "dendrimer"]:
        return None

    if atom_types is None and polymer_type in ["brush", "star"]:
        raise ValueError("atom_types must be provided for brush and star polymers")

    chains = get_chains(bonds, atom_ids, polymer_type, atom_types)
    if not chains:
        return None

    if polymer_type == "linear":
        chain = chains[0]
        if len(chain) < 2:
            return None

        id_to_index = {atom_id: i for i, atom_id in enumerate(atom_ids)}
        pos1 = positions[id_to_index[chain[0]]]
        pos2 = positions[id_to_index[chain[-1]]]
        distance = np.linalg.norm(pos2 - pos1)
        return distance
    else:
        # for brush and star
        distances = []
        id_to_index = {atom_id: i for i, atom_id in enumerate(atom_ids)}
        for chain in chains:
            if len(chain) < 2:
                continue
            pos1 = positions[id_to_index[chain[0]]]
            pos2 = positions[id_to_index[chain[-1]]]
            dist = np.linalg.norm(pos2 - pos1)
            distances.append(dist)
        if not distances:
            return None
        overall = np.mean(distances)
        return {"overall": overall, "chains": distances}


def calculate_diffusion_coefficient(trajectory_positions, steps, timestep=0.01*1000):
    """
    Calculate the diffusion coefficient from trajectory positions.

    Parameters:
    - trajectory_positions: Array of positions over time (n_frames, n_atoms, 3)
    - steps: Array of step values (LAMMPS timesteps)
    - timestep: Time step between frames (in simulation units)

    Returns:
    - D: Diffusion coefficient
    """
    # Calculate center of mass trajectory
    com_trajectory = np.mean(trajectory_positions, axis=1)

    # Calculate mean square displacement
    msd = [0]  # for dts=0 case
    for dt in range(1, len(com_trajectory)):
        displacements = com_trajectory[dt:] - com_trajectory[:-dt]
        sq_disp = np.sum(displacements**2, axis=1)
        msd.append(np.mean(sq_disp))

    msd = np.array(msd)
    dts = (steps[1:] - steps[0]) * timestep  # Convert steps to time

    # Fit linear: MSD = 6 * D * t
    nfit = len(msd)//2
    if len(dts) > 1 and len(msd) > 1:
        slope, intercept, r_value, p_value, std_err = linregress(dts[:nfit], msd[:nfit])  # use the first 1/3 points for fitting
        D = slope / 6.0  # for 3D
    else:
        D = 0.0

    return D, msd, dts
