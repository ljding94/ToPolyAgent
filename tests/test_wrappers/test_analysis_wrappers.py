import unittest
from unittest.mock import patch
from src.wrappers.analysis_wrappers import ComprehensiveAnalysisTool, PlotAnalysisTool


class TestComprehensiveAnalysisTool(unittest.TestCase):
    @patch('src.wrappers.analysis_wrappers.save_analysis_results')
    @patch('src.wrappers.analysis_wrappers.run_complete_analysis')
    @patch('src.wrappers.analysis_wrappers.os.makedirs')
    def test_run_analysis(self, mock_makedirs, mock_analysis, mock_save):
        mock_analysis.return_value = {"metadata": {"success": True}}
        tool = ComprehensiveAnalysisTool()
        result = tool._run("/tmp/test/system_linear.data", "dump.*")
        expected_path = "/tmp/test/system_linear/analysis_results.json"
        self.assertEqual(result, expected_path)
        mock_analysis.assert_called_once_with("/tmp/test/system_linear.data", "dump.*")
        mock_makedirs.assert_called_once_with("/tmp/test/system_linear", exist_ok=True)
        mock_save.assert_called_once_with({"metadata": {"success": True}}, expected_path)


class TestPlotAnalysisTool(unittest.TestCase):
    @patch('src.wrappers.analysis_wrappers.plot_analysis_results')
    def test_plot_analysis(self, mock_plot):
        tool = PlotAnalysisTool()
        result = tool._run("/path/to/analysis_results.json")
        self.assertEqual(result, "/path/to/conformation_analysis.png")
        mock_plot.assert_called_once_with("/path/to", "/path/to/analysis_results.json")


if __name__ == '__main__':
    unittest.main()
