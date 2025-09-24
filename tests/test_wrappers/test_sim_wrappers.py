import unittest
from unittest.mock import patch
from src.wrappers.sim_wrappers import RunLammpsTool


class TestRunLammpsTool(unittest.TestCase):
    @patch('src.wrappers.sim_wrappers.run_lammps')
    def test_run_simulation(self, mock_run):
        mock_run.return_value = {"dump_files": "dump.*", "final_config": "final.data", "log": "log.lammps"}
        tool = RunLammpsTool()
        result = tool._run("/path/to/data.data")
        self.assertEqual(result, {"dump_files": "dump.*", "final_config": "final.data", "log": "log.lammps"})
        mock_run.assert_called_once_with("/path/to/data", "/path/to/data.data", "langevin", {"pp": 0.5, "ss": 0.5, "sp": 0.5}, 100000)

    @patch('src.wrappers.sim_wrappers.run_lammps')
    def test_run_simulation_error(self, mock_run):
        mock_run.side_effect = Exception("Sim failed")
        tool = RunLammpsTool()
        result = tool._run("/path/to/data.data")
        self.assertEqual(result, {"error": "Sim failed"})


if __name__ == '__main__':
    unittest.main()
