import os
from tools.simulation.run_lammps import run_lammps





if __name__ == "__main__":
    # linear polymer
    datafile_path = "/Users/ldq/Work/ToPolyAgent/data/test/system_linear_polymer_30.data"

    # ring polymer
    datafile_path = "/Users/ldq/Work/ToPolyAgent/data/test/system_ring_polymer_20.data"

    # brush polymer
    #datafile_path = "/Users/ldq/Work/ToPolyAgent/data/test/system_brush_polymer_50_0.5_10.data"

    # start polymer
    datafile_path = "/Users/ldq/Work/ToPolyAgent/data/test/system_star_polymer_8_4.data"

    dump_path = datafile_path.replace(".data", "")
    run_lammps(dump_path, datafile_path)
