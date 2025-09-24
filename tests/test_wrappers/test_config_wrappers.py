import unittest
from unittest.mock import patch
from src.wrappers.config_wrappers import PolymerGeneratorTool, PackSolventTool, PlotConfigTool


class TestPolymerGeneratorTool(unittest.TestCase):
    @patch('src.wrappers.config_wrappers.generate_linear_polymer_config')
    def test_generate_linear(self, mock_gen):
        mock_gen.return_value = "/path/to/linear.data"
        tool = PolymerGeneratorTool()
        result = tool._run("linear", {"chain_length": 30})
        self.assertEqual(result, "/path/to/linear.data")
        mock_gen.assert_called_once_with(30, 50.0)

    def test_invalid_type(self):
        tool = PolymerGeneratorTool()
        with self.assertRaises(ValueError):
            tool._run("invalid", {})


class TestPackSolventTool(unittest.TestCase):
    @patch('src.wrappers.config_wrappers.pack_solvent')
    def test_pack_solvent(self, mock_pack):
        mock_pack.return_value = "/path/to/system.data"
        tool = PackSolventTool()
        result = tool._run("/path/to/polymer.data", 0.2)
        self.assertEqual(result, "/path/to/system.data")
        mock_pack.assert_called_once_with("/path/to/polymer.data", 0.2, 50.0)


class TestPlotConfigTool(unittest.TestCase):
    @patch('src.wrappers.config_wrappers.plot_configuration')
    def test_plot_config(self, mock_plot):
        mock_plot.return_value = "/path/to/plot.png"
        tool = PlotConfigTool()
        result = tool._run("/path/to/data.data")
        self.assertEqual(result, "/path/to/plot.png")
        mock_plot.assert_called_once_with("/path/to/data.data", plot_solvent=True)


if __name__ == '__main__':
    unittest.main()
