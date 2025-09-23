import subprocess
import os


def run_lammps(dump_path, datafile_path, thermostat="langevin", interaction_params={"pp": 0.3, "ss": 0.3, "sp": 1.0}, run_steps=10000):
    """
    Runs LAMMPS simulation with the static script and configurable parameters.

    Parameters:
    - dump_path: str, path for output dumps/logs
    - datafile_path: str, path to input datafile
    - thermostat: str, "langevin" or "nose-hoover"
    - ensemble: str, "nvt" or "npt"
    - interaction_params: dict, {"pp": float, "ss": float, "sp": float} for LJ epsilons
    - run_steps: int, number of production steps

    Returns:
    - dict: {"dump_files": str (pattern), "final_config": str, "log": str}
    """
    script_path = os.path.join(os.path.dirname(__file__), "simulation.lammps")
    cmd = [
        "lmp_serial",
        "-in", script_path,
        "-var", "dump_path", dump_path,
        "-var", "datafile_path", datafile_path,
        "-var", "thermostat", thermostat,
        "-var", "eps_pp", str(interaction_params["pp"]),
        "-var", "eps_ss", str(interaction_params["ss"]),
        "-var", "eps_sp", str(interaction_params["sp"]),
        "-var", "prun", str(run_steps),
    ]
    subprocess.run(cmd, check=True)  # Raise error on failure

    return {"dump_files": f"{dump_path}/coord/dump.*.txt", "final_config": f"{dump_path}/final_state.data", "log": f"{dump_path}/log.lammps"}


if __name__ == "__main__":

    datafile_path = "/Users/ldq/Work/ToPolyAgent/data/test/system_linear.data"

    #datafile_path = "/Users/ldq/Work/ToPolyAgent/data/test/system_ring.data"

    #datafile_path = "/Users/ldq/Work/ToPolyAgent/data/test/system_brush.data"

    #datafile_path = "/Users/ldq/Work/ToPolyAgent/data/test/system_star.data"

    #datafile_path = "/Users/ldq/Work/ToPolyAgent/data/test/system_dendrimer.data"

    dump_path = datafile_path.replace(".data", "")
    results = run_lammps(dump_path, datafile_path, thermostat="nosehoover", interaction_params={"pp": 0.3, "ss": 0.3, "sp": 1.0}, run_steps=20000)
    print(results)
