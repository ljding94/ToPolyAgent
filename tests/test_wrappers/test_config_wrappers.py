import os
import sys

# Add project root to path for imports
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(os.path.dirname(current_dir))
sys.path.insert(0, project_root)

import unittest
from unittest.mock import patch
from src.wrappers.config_wrappers import (
    GenerateLinearPolymerTool,
    GenerateRingPolymerTool,
    GenerateBrushPolymerTool,
    GenerateStarPolymerTool,
    GenerateDendrimerTool,
    PackSolventTool
)


class TestGenerateLinearPolymerTool(unittest.TestCase):
    @patch('src.wrappers.config_wrappers.generate_linear_polymer_config')
    def test_generate_linear(self, mock_gen):
        mock_gen.return_value = "/path/to/linear.data"
        tool = GenerateLinearPolymerTool()
        result = tool._run(chain_length=30, output_dir="/tmp")
        self.assertEqual(result, "/path/to/linear.data")
        mock_gen.assert_called_once_with(30, 30.0, "/tmp")

    def test_invalid_chain_length(self):
        tool = GenerateLinearPolymerTool()
        with self.assertRaises(ValueError):
            tool._run(chain_length=1, output_dir="/tmp")

    def test_missing_output_dir(self):
        tool = GenerateLinearPolymerTool()
        with self.assertRaises(ValueError):
            tool._run(chain_length=30)


class TestGenerateRingPolymerTool(unittest.TestCase):
    @patch('src.wrappers.config_wrappers.generate_ring_polymer_config')
    def test_generate_ring(self, mock_gen):
        mock_gen.return_value = "/path/to/ring.data"
        tool = GenerateRingPolymerTool()
        result = tool._run(chain_length=30, output_dir="/tmp")
        self.assertEqual(result, "/path/to/ring.data")
        mock_gen.assert_called_once_with(30, 30.0, "/tmp")

    def test_invalid_chain_length(self):
        tool = GenerateRingPolymerTool()
        with self.assertRaises(ValueError):
            tool._run(chain_length=2, output_dir="/tmp")


class TestGenerateBrushPolymerTool(unittest.TestCase):
    @patch('src.wrappers.config_wrappers.generate_brush_polymer_config')
    def test_generate_brush(self, mock_gen):
        mock_gen.return_value = "/path/to/brush.data"
        tool = GenerateBrushPolymerTool()
        result = tool._run(backbone_length=50, grafting_density=0.3, side_chain_length=10, output_dir="/tmp")
        self.assertEqual(result, "/path/to/brush.data")
        mock_gen.assert_called_once_with(50, 0.3, 10, 30.0, "/tmp")

    def test_invalid_grafting_density(self):
        tool = GenerateBrushPolymerTool()
        with self.assertRaises(ValueError):
            tool._run(backbone_length=50, grafting_density=1.5, side_chain_length=10, output_dir="/tmp")


class TestGenerateStarPolymerTool(unittest.TestCase):
    @patch('src.wrappers.config_wrappers.generate_star_polymer_config')
    def test_generate_star(self, mock_gen):
        mock_gen.return_value = "/path/to/star.data"
        tool = GenerateStarPolymerTool()
        result = tool._run(arm_length=10, num_arms=4, output_dir="/tmp")
        self.assertEqual(result, "/path/to/star.data")
        mock_gen.assert_called_once_with(10, 4, 30.0, "/tmp")

    def test_invalid_num_arms(self):
        tool = GenerateStarPolymerTool()
        with self.assertRaises(ValueError):
            tool._run(arm_length=10, num_arms=2, output_dir="/tmp")


class TestGenerateDendrimerTool(unittest.TestCase):
    @patch('src.wrappers.config_wrappers.generate_dendrimer_config')
    def test_generate_dendrimer(self, mock_gen):
        mock_gen.return_value = "/path/to/dendrimer.data"
        tool = GenerateDendrimerTool()
        result = tool._run(generation=3, branching_factor=3, output_dir="/tmp")
        self.assertEqual(result, "/path/to/dendrimer.data")
        mock_gen.assert_called_once_with(3, 3, 5, 30.0, "/tmp")

    def test_invalid_generation(self):
        tool = GenerateDendrimerTool()
        with self.assertRaises(ValueError):
            tool._run(generation=0, branching_factor=3, output_dir="/tmp")


class TestPackSolventTool(unittest.TestCase):
    @patch('src.wrappers.config_wrappers.pack_solvent')
    def test_pack_solvent(self, mock_pack):
        mock_pack.return_value = "/path/to/system.data"
        tool = PackSolventTool()
        result = tool._run("/path/to/polymer.data", 0.2)
        self.assertEqual(result, "/path/to/system.data")
        mock_pack.assert_called_once_with("/path/to/polymer.data", 0.2, 30.0, None)


if __name__ == '__main__':
    unittest.main()
