import numpy as np
import os
import json


#####################################################################
# some useful helper functions
####################################################################


def generate_gaussian_chain(chain_length, bond_length, angle_mean, angle_std):
    """
    Generate a Gaussian chain with specified bond length and angle distribution.

    Parameters:
    - chain_length: int, number of atoms in the chain
    - bond_length: float, distance between consecutive atoms
    - angle_mean: float, mean bond angle in radians
    - angle_std: float, standard deviation of bond angle in radians

    Returns:
    - list of np.array: positions of atoms in the chain
    """
    if chain_length < 2:
        return [np.array([0.0, 0.0, 0.0])]

    positions = [np.array([0.0, 0.0, 0.0])]  # Start at origin

    # First two atoms: along x-axis
    positions.append(np.array([bond_length, 0.0, 0.0]))

    # Generate remaining atoms
    for i in range(2, chain_length):
        # Sample bond angle from Gaussian distribution
        bond_angle = np.random.normal(angle_mean, angle_std)
        bond_angle = np.clip(bond_angle, 0.1, np.pi - 0.1)  # Keep within reasonable bounds

        # Get the previous two positions to determine the plane
        pos_prev2 = positions[i - 2]
        pos_prev1 = positions[i - 1]

        # Vector from prev2 to prev1 (previous bond)
        prev_bond = pos_prev1 - pos_prev2
        prev_bond = prev_bond / np.linalg.norm(prev_bond)

        # Generate a random torsion angle for 3D flexibility
        torsion_angle = np.random.uniform(0, 2 * np.pi)

        # Create rotation matrix to generate new bond direction
        # First, create a perpendicular vector to the previous bond
        if abs(prev_bond[0]) > 0.1:
            perp = np.array([-prev_bond[1] - prev_bond[2], prev_bond[0], prev_bond[0]])
        else:
            perp = np.array([prev_bond[1], -prev_bond[0] - prev_bond[2], prev_bond[1]])
        perp = perp / np.linalg.norm(perp)

        # Second perpendicular vector
        perp2 = np.cross(prev_bond, perp)

        # New bond direction: rotate around previous bond by bond_angle, then around torsion
        # Start with the previous bond direction rotated by bond_angle in the plane
        cos_angle = np.cos(bond_angle)
        sin_angle = np.sin(bond_angle)

        # Rotate previous bond by bond_angle around perp axis
        new_bond = cos_angle * prev_bond + sin_angle * np.cos(torsion_angle) * perp + sin_angle * np.sin(torsion_angle) * perp2

        new_bond = new_bond / np.linalg.norm(new_bond)  # New position
        new_pos = pos_prev1 + bond_length * new_bond
        positions.append(new_pos)

    return positions


def generate_perpendicular_direction(reference_vector):
    """
    Generate a random unit vector perpendicular to the reference vector.

    Parameters:
    - reference_vector: np.array, the vector to be perpendicular to

    Returns:
    - np.array: unit vector perpendicular to reference_vector
    """
    ref = reference_vector / np.linalg.norm(reference_vector)

    # Create a perpendicular vector
    if abs(ref[0]) > 0.1:
        perp = np.array([-ref[1] - ref[2], ref[0], ref[0]])
    else:
        perp = np.array([ref[1], -ref[0] - ref[2], ref[1]])

    perp = perp / np.linalg.norm(perp)

    # Rotate randomly around the reference vector
    angle = np.random.uniform(0, 2 * np.pi)
    perp2 = np.cross(ref, perp)

    result = perp * np.cos(angle) + perp2 * np.sin(angle)
    return result / np.linalg.norm(result)


def rotation_matrix_from_vectors(vec1, vec2):
    """
    Create a rotation matrix that rotates vec1 to align with vec2.

    Parameters:
    - vec1: np.array, source vector
    - vec2: np.array, target vector

    Returns:
    - np.array: 3x3 rotation matrix
    """
    vec1 = vec1 / np.linalg.norm(vec1)
    vec2 = vec2 / np.linalg.norm(vec2)

    # Cross product to find rotation axis
    cross = np.cross(vec1, vec2)
    cross_norm = np.linalg.norm(cross)

    # If vectors are already aligned (or anti-aligned), return identity or 180° rotation
    if cross_norm < 1e-10:
        if np.dot(vec1, vec2) > 0:
            return np.eye(3)  # Already aligned
        else:
            # 180° rotation around arbitrary perpendicular axis
            perp = np.array([1.0, 0.0, 0.0])
            if abs(np.dot(vec1, perp)) > 0.9:
                perp = np.array([0.0, 1.0, 0.0])
            perp = perp - np.dot(vec1, perp) * vec1
            perp = perp / np.linalg.norm(perp)
            return np.array(
                [
                    [2 * perp[0] ** 2 - 1, 2 * perp[0] * perp[1], 2 * perp[0] * perp[2]],
                    [2 * perp[0] * perp[1], 2 * perp[1] ** 2 - 1, 2 * perp[1] * perp[2]],
                    [2 * perp[0] * perp[2], 2 * perp[1] * perp[2], 2 * perp[2] ** 2 - 1],
                ]
            )

    # Rodrigues' rotation formula
    cos_angle = np.dot(vec1, vec2)
    sin_angle = cross_norm

    # Skew-symmetric matrix for cross product
    K = np.array([[0, -cross[2], cross[1]], [cross[2], 0, -cross[0]], [-cross[1], cross[0], 0]])

    # Rodrigues' formula: R = I + sinθ·K + (1-cosθ)·K²
    identity_matrix = np.eye(3)
    K2 = np.dot(K, K)

    rotation_matrix = identity_matrix + sin_angle * K + (1 - cos_angle) * K2
    return rotation_matrix


#####################################################################
# generate actual polymer configurations
####################################################################


def generate_linear_polymer_config(chain_length, box_size=50.0):
    """
    Generates LAMMPS datafile for a linear polymer (Gaussian chain).

    Parameters:
    - chain_length: int, number of beads in the chain
    - box_size: float, size of simulation box

    Returns:
    - str: path to the generated datafile
    """
    # Metadata
    metadata = {
        "type": "linear",
        "chain_length": chain_length,
        "box_size": box_size
    }

    # Initialize lists
    positions = []
    bonds = []
    angles = []
    atom_id = 1
    bond_id = 1
    angle_id = 1

    # Gaussian chain parameters
    bond_length = 1.0
    angle_mean = 0.0
    angle_std = np.pi * 60.0 / 180.0  # 60 degree standard deviation

    # Generate chain as Gaussian chain
    chain_positions = generate_gaussian_chain(chain_length, bond_length, angle_mean, angle_std)

    # Add atoms to positions
    for pos in chain_positions:
        positions.append((atom_id, 1, 1, 0, pos[0], pos[1], pos[2]))  # Type 1 for all atoms
        atom_id += 1

    # Bonds
    for i in range(chain_length - 1):
        bonds.append((bond_id, 1, i + 1, i + 2))
        bond_id += 1

    # Angles
    for i in range(chain_length - 2):
        angles.append((angle_id, 1, i + 1, i + 2, i + 3))
        angle_id += 1

    # Center the polymer at origin
    if positions:
        com = np.mean([[pos[4], pos[5], pos[6]] for pos in positions], axis=0)
        for i, pos in enumerate(positions):
            new_pos = (pos[0], pos[1], pos[2], pos[3], pos[4] - com[0], pos[5] - com[1], pos[6] - com[2])
            positions[i] = new_pos

    # Write datafile
    datafile_path = os.path.join(os.getcwd(), "polymer_linear.data")
    with open(datafile_path, "w") as f:
        f.write(f"# Metadata: {json.dumps(metadata)}\n")
        f.write("# Linear polymer datafile\n")
        f.write(f"{len(positions)} atoms\n")
        f.write(f"{len(bonds)} bonds\n")
        f.write(f"{len(angles)} angles\n")
        f.write("1 atom types\n")
        f.write("1 bond types\n")
        f.write("1 angle types\n")
        f.write(f"{-box_size/2} {box_size/2} xlo xhi\n")
        f.write(f"{-box_size/2} {box_size/2} ylo yhi\n")
        f.write(f"{-box_size/2} {box_size/2} zlo zhi\n")
        f.write("\nMasses\n\n")
        f.write("1 1.00\n")
        f.write("\n")
        f.write("Atoms #full\n\n")
        for atom in positions:
            f.write(f"{atom[0]} {atom[1]} {atom[2]} {atom[3]} {atom[4]} {atom[5]} {atom[6]}\n")
        f.write("\n")
        f.write("Bonds\n\n")
        for bond in bonds:
            f.write(f"{bond[0]} {bond[1]} {bond[2]} {bond[3]}\n")
        f.write("\n")
        f.write("Angles\n\n")
        for angle in angles:
            f.write(f"{angle[0]} {angle[1]} {angle[2]} {angle[3]} {angle[4]}\n")
        f.write("\n")

    return datafile_path


def generate_ring_polymer_config(ring_length, box_size=50.0):
    """
    Generates LAMMPS datafile for a ring polymer (cyclic polymer).

    Parameters:
    - ring_length: int, number of beads in the ring
    - box_size: float, size of simulation box

    Returns:
    - str: path to the generated datafile
    """
    if ring_length < 3:
        raise ValueError("Ring length must be at least 3 for a valid ring polymer")

    # Metadata
    metadata = {
        "type": "ring",
        "ring_length": ring_length,
        "box_size": box_size
    }

    # Initialize lists
    positions = []
    bonds = []
    angles = []
    atom_id = 1
    bond_id = 1
    angle_id = 1

    # Gaussian chain parameters
    bond_length = 1.0
    angle_mean = 2 * np.pi / ring_length  # Ideal angle for regular polygon
    angle_std = np.pi * 30.0 / 180.0  # 30 degree standard deviation

    # Generate initial chain as Gaussian chain
    chain_positions = generate_gaussian_chain(ring_length, bond_length, angle_mean, angle_std)

    # For ring closure, we need to adjust the positions to form a closed loop
    # Calculate the vector needed to close the ring
    closure_vector = chain_positions[0] - chain_positions[-1]
    closure_distance = np.linalg.norm(closure_vector)

    if closure_distance > bond_length * 1.5:  # If gap is too large, redistribute positions
        # Generate positions in a regular polygon for better ring closure
        radius = (ring_length * bond_length) / (2 * np.pi)  # Circumradius
        for i in range(ring_length):
            angle = 2 * np.pi * i / ring_length
            x = radius * np.cos(angle)
            y = radius * np.sin(angle)
            z = np.random.normal(0, 0.1)  # Small random z-variation
            chain_positions[i] = np.array([x, y, z])

    # Add atoms to positions
    for pos in chain_positions:
        positions.append((atom_id, 1, 1, 0, pos[0], pos[1], pos[2]))  # Type 1 for all atoms
        atom_id += 1

    # Bonds - including the closing bond
    for i in range(ring_length):
        next_i = (i + 1) % ring_length
        bonds.append((bond_id, 1, i + 1, next_i + 1))
        bond_id += 1

    # Angles - all possible angles in the ring
    for i in range(ring_length):
        prev_i = (i - 1) % ring_length
        next_i = (i + 1) % ring_length
        angles.append((angle_id, 1, prev_i + 1, i + 1, next_i + 1))
        angle_id += 1

    # Center the polymer at origin
    if positions:
        com = np.mean([[pos[4], pos[5], pos[6]] for pos in positions], axis=0)
        for i, pos in enumerate(positions):
            new_pos = (pos[0], pos[1], pos[2], pos[3], pos[4] - com[0], pos[5] - com[1], pos[6] - com[2])
            positions[i] = new_pos

    # Write datafile
    datafile_path = os.path.join(os.getcwd(), "polymer_ring.data")
    with open(datafile_path, "w") as f:
        f.write(f"# Metadata: {json.dumps(metadata)}\n")
        f.write("# Ring polymer datafile\n")
        f.write(f"{len(positions)} atoms\n")
        f.write(f"{len(bonds)} bonds\n")
        f.write(f"{len(angles)} angles\n")
        f.write("1 atom types\n")
        f.write("1 bond types\n")
        f.write("1 angle types\n")
        f.write(f"{-box_size/2} {box_size/2} xlo xhi\n")
        f.write(f"{-box_size/2} {box_size/2} ylo yhi\n")
        f.write(f"{-box_size/2} {box_size/2} zlo zhi\n")
        f.write("\nMasses\n\n")
        f.write("1 1.00\n")
        f.write("\n")
        f.write("Atoms #full\n\n")
        for atom in positions:
            f.write(f"{atom[0]} {atom[1]} {atom[2]} {atom[3]} {atom[4]} {atom[5]} {atom[6]}\n")
        f.write("\n")
        f.write("Bonds\n\n")
        for bond in bonds:
            f.write(f"{bond[0]} {bond[1]} {bond[2]} {bond[3]}\n")
        f.write("\n")
        f.write("Angles\n\n")
        for angle in angles:
            f.write(f"{angle[0]} {angle[1]} {angle[2]} {angle[3]} {angle[4]}\n")
        f.write("\n")

    return datafile_path


def generate_brush_polymer_config(backbone_length, grafting_density, side_chain_length, box_size=50.0):
    """
    Generates LAMMPS datafile for a brush polymer.

    Parameters:
    - backbone_length: int, number of beads in backbone
    - grafting_density: float, probability (0-1) of grafting a side chain at each backbone bead
    - side_chain_length: int, number of beads in each side chain
    - box_size: float, size of simulation box

    Returns:
    - str: path to the generated datafile
    """
    # Metadata
    metadata = {
        "type": "brush",
        "backbone_length": backbone_length,
        "grafting_density": grafting_density,
        "side_chain_length": side_chain_length,
        "box_size": box_size
    }

    # Initialize lists
    positions = []
    bonds = []
    angles = []
    atom_id = 1
    bond_id = 1
    angle_id = 1

    # Gaussian chain parameters
    bond_length = 1.0
    backbone_angle_mean = 0.0
    backbone_angle_std = np.pi * 60.0 / 180.0  # 10 degree standard deviation
    side_chain_angle_mean = 0.0
    side_chain_angle_std = np.pi * 60.0 / 180.0  # More flexible side chains

    # Generate backbone as Gaussian chain
    backbone_positions = generate_gaussian_chain(backbone_length, bond_length, backbone_angle_mean, backbone_angle_std)

    # Add backbone atoms to positions
    for pos in backbone_positions:
        positions.append((atom_id, 1, 1, 0, pos[0], pos[1], pos[2]))  # Type 1 for backbone
        atom_id += 1

    # Bonds for backbone
    for i in range(backbone_length - 1):
        bonds.append((bond_id, 1, i + 1, i + 2))
        bond_id += 1

    # Angles for backbone
    for i in range(backbone_length - 2):
        angles.append((angle_id, 1, i + 1, i + 2, i + 3))
        angle_id += 1

    # Side chains - uniformly distributed, Gaussian chains
    num_side_chains = int(backbone_length * grafting_density)
    if num_side_chains > 0:
        # Calculate spacing between side chains for uniform distribution
        spacing = max(1, backbone_length // num_side_chains)
        side_chain_positions = []

        # Distribute side chains uniformly along the backbone
        for i in range(0, backbone_length, spacing):
            if len(side_chain_positions) < num_side_chains:
                side_chain_positions.append(i)

        # If we have fewer side chains than desired, add more at the end
        while len(side_chain_positions) < num_side_chains and len(side_chain_positions) < backbone_length:
            next_pos = (side_chain_positions[-1] + spacing) % backbone_length
            if next_pos not in side_chain_positions:
                side_chain_positions.append(next_pos)

        # Generate side chains at the calculated positions
        for backbone_idx in side_chain_positions:
            back_pos = backbone_positions[backbone_idx]

            # Generate side chain as Gaussian chain
            side_chain_atoms = generate_gaussian_chain(side_chain_length, bond_length, side_chain_angle_mean, side_chain_angle_std)

            # Position the side chain relative to backbone
            # First, determine the direction perpendicular to the backbone at this point
            if backbone_idx == 0:
                # For first backbone atom, use a random perpendicular direction
                backbone_dir = backbone_positions[1] - backbone_positions[0]
            elif backbone_idx == backbone_length - 1:
                # For last backbone atom, use direction from previous
                backbone_dir = backbone_positions[backbone_idx] - backbone_positions[backbone_idx - 1]
            else:
                # For middle atoms, use average of adjacent bonds
                dir1 = backbone_positions[backbone_idx] - backbone_positions[backbone_idx - 1]
                dir2 = backbone_positions[backbone_idx + 1] - backbone_positions[backbone_idx]
                backbone_dir = (dir1 + dir2) / 2

            backbone_dir = backbone_dir / np.linalg.norm(backbone_dir)

            # Generate perpendicular direction for side chain
            perp_dir = generate_perpendicular_direction(backbone_dir)

            # Rotate side chain to align with perp_dir
            if side_chain_length >= 2:
                initial_dir = np.array([1.0, 0.0, 0.0])  # Direction from first to second atom
                rot_matrix = rotation_matrix_from_vectors(initial_dir, perp_dir)
                # Apply rotation to all side chain positions
                side_chain_atoms = [np.dot(rot_matrix, pos) for pos in side_chain_atoms]

            # Position side chain atoms
            for j, chain_pos in enumerate(side_chain_atoms):
                # Offset from backbone and rotate to perpendicular direction
                if j == 0:
                    # First side chain atom is attached to backbone
                    absolute_pos = back_pos + perp_dir * bond_length
                else:
                    # Subsequent atoms follow the chain
                    absolute_pos = back_pos + perp_dir * bond_length + chain_pos

                positions.append((atom_id, 1, 2, 0, absolute_pos[0], absolute_pos[1], absolute_pos[2]))  # Type 2 for side chains
                atom_id += 1
                if j > 0:
                    bonds.append((bond_id, 1, atom_id - 2, atom_id - 1))
                    bond_id += 1

            # Bond from backbone to first side chain bead
            bonds.append((bond_id, 1, backbone_idx + 1, atom_id - side_chain_length))
            bond_id += 1

            # Angles for side chain
            start = atom_id - side_chain_length
            if side_chain_length >= 2:
                # Angle between backbone, first side, second side
                angles.append((angle_id, 1, backbone_idx + 1, start, start + 1))
                angle_id += 1
            # Angles along the side chain
            for j in range(side_chain_length - 2):
                angles.append((angle_id, 1, start + j, start + j + 1, start + j + 2))
                angle_id += 1  # Center the polymer at origin (0,0,0)
    if positions:
        # Calculate center of mass
        com = np.mean([[pos[4], pos[5], pos[6]] for pos in positions], axis=0)

        # Translate all positions to center at origin
        for i, pos in enumerate(positions):
            new_pos = (pos[0], pos[1], pos[2], pos[3], pos[4] - com[0], pos[5] - com[1], pos[6] - com[2])
            positions[i] = new_pos

    # Write datafile
    datafile_path = os.path.join(os.getcwd(), "polymer_brush.data")
    with open(datafile_path, "w") as f:
        f.write(f"# Metadata: {json.dumps(metadata)}\n")
        f.write("# Brush polymer datafile\n")
        f.write(f"{len(positions)} atoms\n")
        f.write(f"{len(bonds)} bonds\n")
        f.write(f"{len(angles)} angles\n")
        f.write("2 atom types\n")
        f.write("1 bond types\n")
        f.write("1 angle types\n")
        f.write(f"{-box_size/2} {box_size/2} xlo xhi\n")
        f.write(f"{-box_size/2} {box_size/2} ylo yhi\n")
        f.write(f"{-box_size/2} {box_size/2} zlo zhi\n")
        f.write("\nMasses\n\n")
        f.write("1 1.00\n")  # Backbone atoms
        f.write("2 1.00\n")  # Side chain atoms
        f.write("\n")
        f.write("Atoms #full\n\n")
        for atom in positions:
            f.write(f"{atom[0]} {atom[1]} {atom[2]} {atom[3]} {atom[4]} {atom[5]} {atom[6]}\n")
        f.write("\n")
        f.write("Bonds\n\n")
        for bond in bonds:
            f.write(f"{bond[0]} {bond[1]} {bond[2]} {bond[3]}\n")
        f.write("\n")
        f.write("Angles\n\n")
        for angle in angles:
            f.write(f"{angle[0]} {angle[1]} {angle[2]} {angle[3]} {angle[4]}\n")
        f.write("\n")

    return datafile_path


def generate_star_polymer_config(arm_length, num_arms, box_size=50.0):
    """
    Generates LAMMPS datafile for a star polymer.

    Parameters:
    - arm_length: int, number of beads in each arm (excluding core)
    - num_arms: int, number of arms radiating from the core
    - box_size: float, size of simulation box

    Returns:
    - str: path to the generated datafile
    """
    if num_arms < 3:
        raise ValueError("Number of arms must be at least 3 for a valid star polymer")
    if arm_length < 1:
        raise ValueError("Arm length must be at least 1")

    # Metadata
    metadata = {
        "type": "star",
        "arm_length": arm_length,
        "num_arms": num_arms,
        "box_size": box_size
    }

    # Initialize lists
    positions = []
    bonds = []
    angles = []
    atom_id = 1
    bond_id = 1
    angle_id = 1

    # Gaussian chain parameters
    bond_length = 1.0
    angle_mean = np.pi * 0.0
    angle_std = np.pi * 60.0 / 180.0  # 30 degree standard deviation

    # Core atom at origin
    core_pos = np.array([0.0, 0.0, 0.0])
    positions.append((atom_id, 1, 1, 0, core_pos[0], core_pos[1], core_pos[2]))  # Type 1 for core
    core_id = atom_id
    atom_id += 1

    # Generate arms radiating from core
    for arm_idx in range(num_arms):
        # Calculate direction for this arm (symmetric arrangement)
        if num_arms == 3:
            # For 3 arms, use 120-degree spacing in xy-plane
            angle_xy = 2 * np.pi * arm_idx / num_arms
            direction = np.array([np.cos(angle_xy), np.sin(angle_xy), 0.0])
        elif num_arms == 4:
            # For 4 arms, use tetrahedral arrangement
            if arm_idx == 0:
                direction = np.array([1.0, 1.0, 1.0])
            elif arm_idx == 1:
                direction = np.array([1.0, -1.0, -1.0])
            elif arm_idx == 2:
                direction = np.array([-1.0, 1.0, -1.0])
            else:  # arm_idx == 3
                direction = np.array([-1.0, -1.0, 1.0])
        else:
            # For more arms, distribute evenly on sphere
            phi = np.arccos(1 - 2 * (arm_idx + 0.5) / num_arms)  # Golden angle spiral
            theta = np.pi * (1 + np.sqrt(5)) * arm_idx
            direction = np.array([np.sin(phi) * np.cos(theta), np.sin(phi) * np.sin(theta), np.cos(phi)])

        direction = direction / np.linalg.norm(direction)

        # Generate arm as Gaussian chain
        arm_positions = generate_gaussian_chain(arm_length, bond_length, angle_mean, angle_std)

        # Compute the direction of the first bond in the Gaussian chain
        if len(arm_positions) > 1:
            # Use the vector from first to second atom as the chain direction
            first_bond_direction = arm_positions[1] - arm_positions[0]
            first_bond_direction = first_bond_direction / np.linalg.norm(first_bond_direction)
        else:
            first_bond_direction = np.array([1.0, 0.0, 0.0])  # Default direction for single atom

        # Create rotation matrix to align first_bond_direction with target direction
        rotation_matrix = rotation_matrix_from_vectors(first_bond_direction, direction)

        # Rotate and position arm atoms
        for i, arm_pos in enumerate(arm_positions):
            # Rotate the atom position
            rotated_pos = np.dot(rotation_matrix, arm_pos)

            # Position relative to core (first atom should be at bond_length from core)
            if i == 0:
                absolute_pos = core_pos + direction * bond_length
                # Store the position of the first atom after rotation for shifting
                first_atom_offset = rotated_pos
            else:
                # Shift all atoms so that the first atom aligns with the target position
                shifted_pos = rotated_pos - first_atom_offset + direction * bond_length
                absolute_pos = core_pos + shifted_pos

            positions.append((atom_id, 1, 2, 0, absolute_pos[0], absolute_pos[1], absolute_pos[2]))  # Type 2 for arm atoms
            atom_id += 1

            # Bonds within arm
            if i > 0:
                bonds.append((bond_id, 1, atom_id - 2, atom_id - 1))
                bond_id += 1

        # Bond from core to first arm atom
        first_arm_atom_id = atom_id - arm_length
        bonds.append((bond_id, 1, core_id, first_arm_atom_id))
        bond_id += 1

        # Angles for arm
        if arm_length >= 2:
            # Angle between core, first arm atom, second arm atom
            angles.append((angle_id, 1, core_id, first_arm_atom_id, first_arm_atom_id + 1))
            angle_id += 1

        # Angles along the arm
        for i in range(arm_length - 2):
            angles.append((angle_id, 1, first_arm_atom_id + i, first_arm_atom_id + i + 1, first_arm_atom_id + i + 2))
            angle_id += 1

    # Center the polymer at origin (already at origin, but ensure numerical precision)
    if positions:
        com = np.mean([[pos[4], pos[5], pos[6]] for pos in positions], axis=0)
        for i, pos in enumerate(positions):
            new_pos = (pos[0], pos[1], pos[2], pos[3], pos[4] - com[0], pos[5] - com[1], pos[6] - com[2])
            positions[i] = new_pos

    # Write datafile
    datafile_path = os.path.join(os.getcwd(), "polymer_star.data")
    with open(datafile_path, "w") as f:
        f.write(f"# Metadata: {json.dumps(metadata)}\n")
        f.write("# Star polymer datafile\n")
        f.write(f"{len(positions)} atoms\n")
        f.write(f"{len(bonds)} bonds\n")
        f.write(f"{len(angles)} angles\n")
        f.write("2 atom types\n")
        f.write("1 bond types\n")
        f.write("1 angle types\n")
        f.write(f"{-box_size/2} {box_size/2} xlo xhi\n")
        f.write(f"{-box_size/2} {box_size/2} ylo yhi\n")
        f.write(f"{-box_size/2} {box_size/2} zlo zhi\n")
        f.write("\nMasses\n\n")
        f.write("1 1.00\n")  # Core atoms
        f.write("2 1.00\n")  # Arm atoms
        f.write("\n")
        f.write("Atoms #full\n\n")
        for atom in positions:
            f.write(f"{atom[0]} {atom[1]} {atom[2]} {atom[3]} {atom[4]} {atom[5]} {atom[6]}\n")
        f.write("\n")
        f.write("Bonds\n\n")
        for bond in bonds:
            f.write(f"{bond[0]} {bond[1]} {bond[2]} {bond[3]}\n")
        f.write("\n")
        f.write("Angles\n\n")
        for angle in angles:
            f.write(f"{angle[0]} {angle[1]} {angle[2]} {angle[3]} {angle[4]}\n")
        f.write("\n")

    return datafile_path



def generate_dendrimer_config(generations, branching_factor, spacer=5, box_size=50.0):
    """
    Generates LAMMPS datafile for a dendrimer.
    Parameters:
    - generations: int, number of generations (layers) in the dendrimer
    - branching_factor: int, number of branches per node (typically 2 or 3)
    - spacer: int, number of atoms per branch segment (1 = current behavior, >1 adds more atoms per line)
    - box_size: float, size of simulation box
    Returns:
    - str: path to the generated datafile
    """
    if generations < 1:
        raise ValueError("Generations must be at least 1")
    if branching_factor < 2:
        raise ValueError("Branching factor must be at least 2")
    if spacer < 1:
        raise ValueError("Spacer must be at least 1")

    # Metadata
    metadata = {
        "type": "dendrimer",
        "generations": generations,
        "branching_factor": branching_factor,
        "spacer": spacer,
        "box_size": box_size
    }

    # Initialize lists
    positions = []
    bonds = []
    atom_id = 1
    bond_id = 1
    # Angles will be added later

    # Bond length (scale slightly with generation for better spreading, optional)
    base_bond_length = 1.0

    # Core atom (generation 0)
    core_pos = np.array([0.0, 0.0, 0.0])
    positions.append((atom_id, 1, 1, 0, core_pos[0], core_pos[1], core_pos[2]))  # Molecule 1, type 1, charge 0
    atom_id += 1

    # Keep track of branch points (parents) in each generation
    # Each entry: (atom_id, position, incoming_direction)
    gen_branch_points = [[(1, core_pos, np.array([0.0, 0.0, 1.0]))]]  # Generation 0: core has arbitrary incoming direction

    # Build dendrimer generation by generation
    for gen in range(1, generations + 1):
        gen_branch_points.append([])
        parent_atoms = gen_branch_points[gen - 1]
        for parent_id, parent_pos, incoming_dir in parent_atoms:
            # Number of new branches from this parent
            if gen == 1:
                num_new = branching_factor  # From core
            else:
                num_new = branching_factor - 1  # From branch points

            # Generate directions distributed around the incoming direction
            directions = []
            if gen == 1:
                # For the first generation from core, distribute uniformly in 3D space
                if num_new == 1:
                    # Single branch: random direction
                    dir_vec = np.random.normal(0, 1, 3)
                    dir_vec /= np.linalg.norm(dir_vec)
                    directions.append(dir_vec)
                else:
                    # Multiple branches: use golden spiral for uniform distribution on sphere
                    golden_ratio = (1 + 5**0.5) / 2
                    for i in range(num_new):
                        theta = 2 * np.pi * i / golden_ratio
                        phi = np.arccos(1 - 2 * (i + 0.5) / num_new)
                        dir_vec = np.array([
                            np.sin(phi) * np.cos(theta),
                            np.sin(phi) * np.sin(theta),
                            np.cos(phi)
                        ])
                        directions.append(dir_vec)
            else:
                # For subsequent generations, distribute around the incoming direction
                if num_new == 1:
                    # Single branch: random direction perpendicular to incoming
                    perp_dir = generate_perpendicular_direction(incoming_dir)
                    directions.append(perp_dir)
                else:
                    # Multiple branches: distribute around incoming direction using spherical coordinates
                    for i in range(num_new):
                        # Use golden angle for even distribution
                        golden_angle = 2 * np.pi * i / ((1 + np.sqrt(5)) / 2)
                        # Elevation angle from the incoming direction plane
                        elevation = np.arccos(1 - 2 * (i + 0.5) / num_new) - np.pi/2  # Center around equator

                        # Create rotation matrix to align with incoming direction
                        # First, create a reference direction perpendicular to incoming
                        ref_dir = generate_perpendicular_direction(incoming_dir)
                        ref_dir2 = np.cross(incoming_dir, ref_dir)

                        # Rotate around incoming direction by golden_angle, then tilt by elevation
                        cos_golden = np.cos(golden_angle)
                        sin_golden = np.sin(golden_angle)
                        cos_elev = np.cos(elevation)
                        sin_elev = np.sin(elevation)

                        # Direction in the plane perpendicular to incoming
                        perp_component = cos_golden * ref_dir + sin_golden * ref_dir2
                        # Combine with elevation
                        dir_vec = cos_elev * perp_component + sin_elev * incoming_dir
                        dir_vec = dir_vec / np.linalg.norm(dir_vec)
                        directions.append(dir_vec)

            # For each new branch: add a chain of 'spacer' atoms along the direction (straight line)
            for dir_vec in directions:
                current_id = parent_id
                current_pos = parent_pos
                branch_direction = dir_vec  # Store the outgoing direction for the branch point
                for s in range(spacer):
                    # Scale bond length mildly with gen to fan out (helps with overlaps)
                    bond_length = base_bond_length * (1 + 0.02 * gen)  # Optional: adjust 0.1 for more/less spreading
                    new_pos = current_pos + branch_direction * bond_length
                    positions.append((atom_id, 1, 2, 0, new_pos[0], new_pos[1], new_pos[2]))  # Type 2
                    # Bond to previous
                    bonds.append((bond_id, 1, current_id, atom_id))
                    bond_id += 1
                    # Update current
                    current_id = atom_id
                    current_pos = new_pos
                    atom_id += 1
                # The last atom in the chain is the new branch point
                gen_branch_points[gen].append((current_id, current_pos, branch_direction))

    # Now add angles correctly using adjacency list
    angles = []
    angle_id = 1
    # Build neighbors dict
    neighbors = {i: [] for i in range(1, atom_id)}
    for _, _, a1, a2 in bonds:
        neighbors[a1].append(a2)
        neighbors[a2].append(a1)
    # For each atom (central), add angles for pairs of neighbors
    for central in neighbors:
        neigh_list = sorted(neighbors[central])  # Sort to avoid duplicates
        for i in range(len(neigh_list)):
            for j in range(i + 1, len(neigh_list)):
                angles.append((angle_id, 1, neigh_list[i], central, neigh_list[j]))
                angle_id += 1

    # Center the dendrimer at origin
    if positions:
        com = np.mean([[p[4], p[5], p[6]] for p in positions], axis=0)
        positions = [(p[0], p[1], p[2], p[3], p[4] - com[0], p[5] - com[1], p[6] - com[2]) for p in positions]

    # Write datafile
    datafile_path = os.path.join(os.getcwd(), "polymer_dendrimer.data")
    with open(datafile_path, "w") as f:
        f.write(f"# Metadata: {json.dumps(metadata)}\n")
        f.write("# Dendrimer datafile\n")
        f.write(f"{len(positions)} atoms\n")
        f.write(f"{len(bonds)} bonds\n")
        f.write(f"{len(angles)} angles\n")
        f.write("2 atom types\n")
        f.write("1 bond types\n")
        f.write("1 angle types\n")
        f.write(f"{-box_size/2} {box_size/2} xlo xhi\n")
        f.write(f"{-box_size/2} {box_size/2} ylo yhi\n")
        f.write(f"{-box_size/2} {box_size/2} zlo zhi\n")
        f.write("\nMasses\n\n")
        f.write("1 1.00\n")  # Core
        f.write("2 1.00\n")  # Others
        f.write("Atoms #full\n\n")
        for atom in positions:
            f.write(f"{atom[0]} {atom[1]} {atom[2]} {atom[3]} {atom[4]:.6f} {atom[5]:.6f} {atom[6]:.6f}\n")
        f.write("\nBonds\n\n")
        for bond in bonds:
            f.write(f"{bond[0]} {bond[1]} {bond[2]} {bond[3]}\n")
        f.write("\nAngles\n\n")
        for angle in angles:
            f.write(f"{angle[0]} {angle[1]} {angle[2]} {angle[3]} {angle[4]}\n")

    return datafile_path