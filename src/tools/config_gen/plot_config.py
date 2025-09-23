import numpy as np
import MDAnalysis as mda
import pyvista as pv
import os


def plot_configuration(data_file_path, solvent_size_factor=0.5, interactive=False, plot_box=True, plot_solvent=True):
    """
    Load a LAMMPS .data file and plot the configuration using PyVista.
    Shows atoms as spheres and bonds as lines for better polymer visualization.

    Parameters:
    - data_file_path: Path to the .data file
    - solvent_size_factor: Factor to scale down solvent particle sizes (default: 0.5)
    - interactive: If True, display interactive plot window; if False, save screenshot (default: False)
    - plot_box: If True, plot the simulation box boundaries (default: True, but disabled if only polymer atoms)
    - plot_solvent: If True, plot solvent atoms; if False, only plot polymer atoms (default: True)
    """

    # Generate output image path
    base_path = os.path.splitext(data_file_path)[0]
    if plot_solvent:
        output_image_path = base_path + '.png'
    else:
        output_image_path = base_path + '_nosolvent.png'

    # Load the data using MDAnalysis
    u = mda.Universe(data_file_path, format='DATA')

    # Get atom positions and types
    positions = u.atoms.positions
    types = u.atoms.types.astype(int)

    # Filter atoms based on plot_solvent option
    if not plot_solvent:
        atom_mask = types <= 2  # Only polymer atoms
        positions = positions[atom_mask]
        types = types[atom_mask]

    # Determine if to plot box: disable if only polymer atoms
    if not np.any(types > 2):
        plot_box = False

    # Create PyVista plotter
    plotter = pv.Plotter(off_screen=not interactive, window_size=(900, 900))

    # Enable better lighting
    #plotter.enable_eye_dome_lighting()
    #light = pv.Light()
    #light.set_direction_angle(30, -30)
    #plotter.add_light(light)

    # Get unique atom types and assign colors
    unique_types = np.unique(types)
    #colors = ['blue', 'red', 'yellow', 'green', 'orange', 'purple', 'brown', 'pink', 'gray', 'olive', 'cyan']
    colors = ["black","royalblue", "tomato", "gold", "limegreen"]

    # Plot atoms by type with physical sphere sizes
    for i, atom_type in enumerate(unique_types):
        mask = types == atom_type
        # Set radius: polymer types (1-2) have radius 0.75, solvent 0.375
        radius = 0.5 if atom_type <= 2 else 0.5 * solvent_size_factor
        points = pv.PolyData(positions[mask])
        spheres = points.glyph(geom=pv.Sphere(radius=radius), scale=False)
        plotter.add_mesh(spheres, color=colors[atom_type], label=f'Atom Type {atom_type}', show_edges=False)

    # Plot bonds for polymer atoms (types <= 2) as cylinders
    full_types = u.atoms.types.astype(int)  # Use full types for bond filtering
    polymer_mask = full_types <= 2
    if hasattr(u, 'bonds') and len(u.bonds) > 0:
        # Read box bounds for PBC handling
        bounds = []
        with open(data_file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if 'xlo xhi' in line:
                    parts = line.split()
                    bounds.extend([float(parts[0]), float(parts[1])])
                elif 'ylo yhi' in line:
                    parts = line.split()
                    bounds.extend([float(parts[0]), float(parts[1])])
                elif 'zlo zhi' in line:
                    parts = line.split()
                    bounds.extend([float(parts[0]), float(parts[1])])
                    break
        if len(bounds) == 6:
            Lx = bounds[1] - bounds[0]
            Ly = bounds[3] - bounds[2]
            Lz = bounds[5] - bounds[4]
        else:
            Lx = Ly = Lz = 0  # Fallback if bounds not found

        # Get polymer atom IDs for filtering bonds
        polymer_ids = set(u.atoms.ids[polymer_mask])

        for bond in u.bonds:
            # Check if both atoms in the bond are polymer atoms
            atom1_id, atom2_id = bond.atoms[0].id, bond.atoms[1].id
            if atom1_id in polymer_ids and atom2_id in polymer_ids:
                pos1 = bond.atoms[0].position
                pos2 = bond.atoms[1].position
                vec = pos2 - pos1
                original_distance = np.linalg.norm(vec)

                if original_distance > 2:  # Threshold for PBC crossing
                    # Calculate shift for minimum image
                    shift = np.zeros(3)
                    Ls = [Lx, Ly, Lz]
                    for i in range(3):
                        L = Ls[i]
                        if vec[i] > L/2:
                            shift[i] = -L
                        elif vec[i] < -L/2:
                            shift[i] = L

                    # Mirror positions
                    mirror_pos2 = pos2 + shift
                    mirror_pos1 = pos1 - shift

                    # Plot first bond: pos1 to mirror_pos2
                    vec1 = mirror_pos2 - pos1
                    height1 = np.linalg.norm(vec1)
                    if height1 > 0:
                        direction1 = vec1 / height1
                        center1 = (pos1 + mirror_pos2) / 2
                        cylinder1 = pv.Cylinder(center=center1, direction=direction1, radius=0.2, height=height1)
                        plotter.add_mesh(cylinder1, color='gray')

                    # Plot second bond: mirror_pos1 to pos2
                    vec2 = pos2 - mirror_pos1
                    height2 = np.linalg.norm(vec2)
                    if height2 > 0:
                        direction2 = vec2 / height2
                        center2 = (mirror_pos1 + pos2) / 2
                        cylinder2 = pv.Cylinder(center=center2, direction=direction2, radius=0.2, height=height2)
                        plotter.add_mesh(cylinder2, color='gray')
                else:
                    # Normal bond plotting
                    height = original_distance
                    if height > 0:
                        direction = vec / height
                        center = (pos1 + pos2) / 2
                        cylinder = pv.Cylinder(center=center, direction=direction, radius=0.2, height=height)
                        plotter.add_mesh(cylinder, color='gray')

    # Plot simulation box if enabled
    if plot_box:
        # Read box bounds directly from the .data file
        bounds = []
        with open(data_file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if 'xlo xhi' in line:
                    parts = line.split()
                    bounds.extend([float(parts[0]), float(parts[1])])
                elif 'ylo yhi' in line:
                    parts = line.split()
                    bounds.extend([float(parts[0]), float(parts[1])])
                elif 'zlo zhi' in line:
                    parts = line.split()
                    bounds.extend([float(parts[0]), float(parts[1])])
                    break  # Assuming zlo zhi is the last box line
        if len(bounds) == 6:
            box = pv.Box(bounds=bounds)
            plotter.add_mesh(box, style='wireframe', color='black', line_width=1)

    # Set background color
    # plotter.set_background('white')

    # Tight layout
    plotter.reset_camera()

    # Set camera position for better view
    plotter.camera_position = 'iso'

    if interactive:
        # Show the plot interactively
        plotter.show()
        print("Configuration plot displayed interactively")
    else:
        # Save screenshot
        plotter.screenshot(output_image_path, scale=3)
        # Close plotter
        plotter.close()
        print(f"Configuration plot saved to {output_image_path}")



