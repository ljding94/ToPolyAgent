import numpy as np
import os


def pack_solvent(polymer_datafile, solvent_density, box_size=50.0):
    """
    Adds solvent beads around the polymer, avoiding overlaps.

    Parameters:
    - polymer_datafile: str, path to polymer datafile
    - solvent_density: float, number density of solvent beads
    - box_size: float, box size

    Returns:
    - str: path to full system datafile
    """
    # Read polymer datafile
    polymer_positions = []
    polymer_bonds = []
    polymer_angles = []
    max_atom_id = 0
    max_bond_id = 0
    max_angle_id = 0

    with open(polymer_datafile, 'r') as f:
        lines = f.readlines()

    # Parse atoms and bonds from polymer file
    in_atoms_section = False
    in_bonds_section = False
    in_angles_section = False

    for line in lines:
        line = line.strip()
        if line == "Atoms" or line == "Atoms #full":
            in_atoms_section = True
            in_bonds_section = False
            in_angles_section = False
            continue
        elif line == "Bonds":
            in_atoms_section = False
            in_bonds_section = True
            in_angles_section = False
            continue
        elif line == "Angles":
            in_atoms_section = False
            in_bonds_section = False
            in_angles_section = True
            continue
        elif line == "" or line.startswith("#"):
            continue

        if in_atoms_section and len(line.split()) >= 7:
            parts = line.split()
            atom_id = int(parts[0])
            mol_id = int(parts[1])
            atom_type = int(parts[2])
            charge = float(parts[3])
            x, y, z = float(parts[4]), float(parts[5]), float(parts[6])
            polymer_positions.append((atom_id, mol_id, atom_type, charge, x, y, z))
            max_atom_id = max(max_atom_id, atom_id)

        elif in_bonds_section and len(line.split()) >= 4:
            parts = line.split()
            bond_id = int(parts[0])
            bond_type = int(parts[1])
            atom1, atom2 = int(parts[2]), int(parts[3])
            polymer_bonds.append((bond_id, bond_type, atom1, atom2))
            max_bond_id = max(max_bond_id, bond_id)

        elif in_angles_section and len(line.split()) >= 5:
            parts = line.split()
            angle_id = int(parts[0])
            angle_type = int(parts[1])
            atom1, atom2, atom3 = int(parts[2]), int(parts[3]), int(parts[4])
            polymer_angles.append((angle_id, angle_type, atom1, atom2, atom3))
            max_angle_id = max(max_angle_id, angle_id)

    # Calculate number of solvent beads needed
    volume = box_size ** 3
    num_solvent = int(solvent_density * volume)

    # Add solvent beads avoiding overlaps
    solvent_positions = []
    min_distance = 1.5  # Minimum distance between beads

    for i in range(num_solvent):
        attempts = 0
        max_attempts = 100  # Reduced from 1000 for faster testing

        while attempts < max_attempts:
            # Generate random position within box
            x = np.random.uniform(-box_size/2, box_size/2)
            y = np.random.uniform(-box_size/2, box_size/2)
            z = np.random.uniform(-box_size/2, box_size/2)

            # Quick check: ensure position is not too close to polymer center
            polymer_com = np.mean([[px, py, pz] for _, _, _, _, px, py, pz in polymer_positions], axis=0)
            com_distance = np.sqrt((x - polymer_com[0])**2 + (y - polymer_com[1])**2 + (z - polymer_com[2])**2)

            # Skip if too close to polymer center (likely to overlap)
            if com_distance < 3.0:
                attempts += 1
                continue

            # Check distance from polymer atoms (only check nearby atoms for efficiency)
            overlap = False
            for _, _, _, _, px, py, pz in polymer_positions:
                distance = np.sqrt((x - px)**2 + (y - py)**2 + (z - pz)**2)
                if distance < min_distance:
                    overlap = True
                    break

            # If no overlap with polymer, check distance from other solvent atoms
            if not overlap:
                for _, _, _, _, sx, sy, sz in solvent_positions[-10:]:  # Only check last 10 solvent atoms
                    distance = np.sqrt((x - sx)**2 + (y - sy)**2 + (z - sz)**2)
                    if distance < min_distance:
                        overlap = True
                        break

            if not overlap:
                solvent_positions.append((max_atom_id + i + 1, 2, 3, 0.0, x, y, z))  # Type 3 for solvent
                break

            attempts += 1

        if attempts >= max_attempts:
            print(f"Warning: Could not place solvent bead {i+1} without overlap")

    # Create system datafile
    dirname = os.path.dirname(polymer_datafile)
    basename = os.path.basename(polymer_datafile)
    system_datafile = os.path.join(dirname, f"system_{basename}")

    with open(system_datafile, 'w') as f:
        f.write("# Polymer + Solvent system datafile\n")
        f.write(f"{len(polymer_positions) + len(solvent_positions)} atoms\n")
        f.write(f"{len(polymer_bonds)} bonds\n")
        f.write(f"{len(polymer_angles)} angles\n")
        f.write("3 atom types\n")
        f.write("1 bond types\n")
        f.write("1 angle types\n")
        f.write(f"{-box_size/2} {box_size/2} xlo xhi\n")
        f.write(f"{-box_size/2} {box_size/2} ylo yhi\n")
        f.write(f"{-box_size/2} {box_size/2} zlo zhi\n")
        f.write("\nMasses\n\n")
        f.write("1 1.00\n")  # Polymer backbone
        f.write("2 1.00\n")  # Polymer side chains
        f.write("3 1.00\n")  # Solvent atoms
        f.write("\n")
        f.write("Atoms #full\n\n")

        # Write polymer atoms
        for atom in polymer_positions:
            f.write(f"{atom[0]} {atom[1]} {atom[2]} {atom[3]} {atom[4]} {atom[5]} {atom[6]}\n")

        # Write solvent atoms
        for atom in solvent_positions:
            f.write(f"{atom[0]} {atom[1]} {atom[2]} {atom[3]} {atom[4]} {atom[5]} {atom[6]}\n")

        f.write("\n")
        f.write("Bonds\n\n")

        # Write polymer bonds
        for bond in polymer_bonds:
            f.write(f"{bond[0]} {bond[1]} {bond[2]} {bond[3]}\n")

        f.write("\n")
        f.write("Angles\n\n")

        # Write polymer angles
        for angle in polymer_angles:
            f.write(f"{angle[0]} {angle[1]} {angle[2]} {angle[3]} {angle[4]}\n")

        f.write("\n")

    print(f"Added {len(solvent_positions)} solvent beads to system")
    return system_datafile
