import subprocess
import os


def run_lammps(dump_path, datafile_path):
    """
    Runs LAMMPS simulation with the static script.

    Parameters:
    - datafile_path: str, path to datafile
    - temp: float, temperature
    - steps: int, number of steps

    Returns:
    - dict: paths to dump and log files
    """

    # lmp_serial -in simulation.lammps -v data_file "/Users/ldq/Work/ToPolyAgent/data/test/system_linear_polymer_30.data"
    script_path = os.path.join(os.path.dirname(__file__), "simulation.lammps")
    cmd = ["lmp_serial",
           "-in", script_path,
           "-var", "dump_path", dump_path,
           "-var", "datafile_path", datafile_path]
    subprocess.run(cmd)


if __name__ == "__main__":

    datafile_path = "/Users/ldq/Work/ToPolyAgent/data/test/system_linear_polymer_30.data"

    datafile_path = "/Users/ldq/Work/ToPolyAgent/data/test/system_brush_polymer_50_0.5_10.data"

    dump_path = datafile_path.replace(".data", "")
    run_lammps(dump_path, datafile_path)
