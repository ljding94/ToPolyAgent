import os
import sys
import unittest
from unittest.mock import patch, MagicMock

# Add project root to path for imports
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(os.path.dirname(current_dir))
sys.path.insert(0, project_root)

from src.wrappers.workflow_wrappers import FullWorkflowTool, run_full_workflow


class TestFullWorkflowTool(unittest.TestCase):

    def setUp(self):
        """Set up test fixtures."""
        self.config = {
            "output_dir": "/tmp/test_output",
            "polymer_type": "linear",
            "box_size": 20.0,
            "solvent_density": 0.3,
            "run_steps": 1000,
            "thermostat": "langevin",
            "interaction_params": {"pp": 0.3, "ss": 0.3, "sp": 1.5},
            "polymer_params": {"chain_length": 20}
        }

    @patch('src.wrappers.workflow_wrappers.GenerateLinearPolymerTool')
    @patch('src.wrappers.workflow_wrappers.RunLammpsTool')
    @patch('src.wrappers.workflow_wrappers.ComprehensiveAnalysisTool')
    @patch('src.wrappers.workflow_wrappers.PlotAnalysisTool')
    @patch('src.wrappers.workflow_wrappers.os.makedirs')
    @patch('src.wrappers.workflow_wrappers.os.chdir')
    @patch('builtins.open', new_callable=unittest.mock.mock_open)
    @patch('json.dump')
    def test_run_full_workflow_linear(self, mock_json_dump, mock_open, mock_chdir, mock_makedirs, mock_plot_analysis,
                                      mock_analysis, mock_sim, mock_gen):
        """Test the full workflow with a linear polymer."""
        # Mock the tool instances
        mock_gen_instance = MagicMock()
        mock_gen.return_value = mock_gen_instance
        mock_gen_instance._run.return_value = "/path/to/system.data"

        mock_sim_instance = MagicMock()
        mock_sim.return_value = mock_sim_instance
        mock_sim_instance._run.return_value = {
            "dump_files": "/path/to/dumps/*.txt",
            "final_config": "/path/to/final.data",
            "final_polymer_config": "/path/to/final_polymer.data",
            "log": "/path/to/log.lammps"
        }

        mock_analysis_instance = MagicMock()
        mock_analysis.return_value = mock_analysis_instance
        mock_analysis_instance._run.return_value = "/path/to/analysis_results.json"

        mock_plot_analysis_instance = MagicMock()
        mock_plot_analysis.return_value = mock_plot_analysis_instance
        mock_plot_analysis_instance._run.return_value = "/path/to/analysis_plot.png"

        # Run the workflow
        result = run_full_workflow(self.config)

        # Verify the result structure
        self.assertEqual(result["status"], "completed")
        self.assertIn("output_paths", result)
        self.assertIn("steps", result)
        self.assertIn("system", result["output_paths"])
        self.assertIn("analysis_results", result["output_paths"])

        # Verify all steps were called
        self.assertEqual(result["steps"]["configuration_generation"], "completed")
        self.assertEqual(result["steps"]["simulation"], "completed")
        self.assertEqual(result["steps"]["analysis"], "completed")

    def test_run_full_workflow_invalid_polymer_type(self):
        """Test workflow with invalid polymer type."""
        invalid_config = self.config.copy()
        invalid_config["polymer_type"] = "invalid"

        with self.assertRaises(ValueError):
            run_full_workflow(invalid_config)

    def test_run_full_workflow_invalid_thermostat(self):
        """Test workflow with invalid thermostat."""
        invalid_config = self.config.copy()
        invalid_config["thermostat"] = "invalid"

        with self.assertRaises(ValueError) as context:
            run_full_workflow(invalid_config)

        self.assertIn("Invalid thermostat", str(context.exception))
        self.assertIn("nosehoover", str(context.exception))
        self.assertIn("langevin", str(context.exception))

    def test_run_full_workflow_missing_output_dir(self):
        """Test that output_dir is required."""
        # This should work since we have default output dir creation
        result = run_full_workflow(self.config)
        self.assertIn("output_dir", result)

    @patch('src.wrappers.workflow_wrappers.run_full_workflow')
    def test_full_workflow_tool(self, mock_workflow):
        """Test the FullWorkflowTool class."""
        mock_workflow.return_value = {"status": "completed", "output_paths": {}}

        tool = FullWorkflowTool()
        result = tool._run(
            output_dir=self.config["output_dir"],
            polymer_type=self.config["polymer_type"],
            polymer_params_json='{"chain_length": 20}',
            box_size=self.config["box_size"],
            solvent_density=self.config["solvent_density"],
            run_steps=self.config["run_steps"],
            thermostat=self.config["thermostat"],
            interaction_params_json='{"pp": 0.3, "ss": 0.3, "sp": 1.5}'
        )

        # Verify that run_full_workflow was called with the correct config
        expected_config = {
            'output_dir': self.config["output_dir"],
            'polymer_type': self.config["polymer_type"],
            'polymer_params': {"chain_length": 20},
            'box_size': self.config["box_size"],
            'solvent_density': self.config["solvent_density"],
            'run_steps': self.config["run_steps"],
            'thermostat': self.config["thermostat"],
            'interaction_params': {"pp": 0.3, "ss": 0.3, "sp": 1.5}
        }
        mock_workflow.assert_called_once_with(expected_config)
        self.assertEqual(result["status"], "completed")

    def test_full_workflow_tool_invalid_thermostat(self):
        """Test the FullWorkflowTool class with invalid thermostat."""
        tool = FullWorkflowTool()

        with self.assertRaises(ValueError) as context:
            tool._run(
                output_dir=self.config["output_dir"],
                polymer_type=self.config["polymer_type"],
                polymer_params_json='{"chain_length": 20}',
                thermostat="invalid_thermostat"
            )

        self.assertIn("Invalid thermostat", str(context.exception))
        self.assertIn("nosehoover", str(context.exception))
        self.assertIn("langevin", str(context.exception))

    def test_config_validation(self):
        """Test that config requires polymer_type."""
        invalid_config = self.config.copy()
        del invalid_config["polymer_type"]

        with self.assertRaises(KeyError):
            run_full_workflow(invalid_config)


if __name__ == '__main__':
    unittest.main()
