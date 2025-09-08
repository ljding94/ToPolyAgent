import unittest
import os
from tools.config_gen.generate_polymer import generate_brush_polymer_config, generate_linear_polymer_config, generate_ring_polymer_config, generate_star_polymer_config


class TestGeneratePolymer(unittest.TestCase):

    def setUp(self):
        # Create test directory
        self.test_dir = "/Users/ldq/Work/ToPolyAgent/data/test"
        os.makedirs(self.test_dir, exist_ok=True)
        self.original_cwd = os.getcwd()

    def tearDown(self):
        # Clean up test files
        os.chdir(self.original_cwd)
        # Temporarily disabled cleanup for inspection
        # for filename in os.listdir(self.test_dir):
        #     if filename.startswith('brush_polymer_'):
        #         os.remove(os.path.join(self.test_dir, filename))

    def test_generate_linear_polymer_config(self):
        # Test linear polymer generation
        chain_length = 50
        box_size = 50.0

        os.chdir(self.test_dir)

        datafile_path = generate_linear_polymer_config(chain_length, box_size)

        # Basic checks
        self.assertTrue(os.path.exists(datafile_path))
        self.assertTrue(os.path.getsize(datafile_path) > 0)

    def test_generate_ring_polymer_config(self):
        # Test ring polymer generation
        ring_length = 20
        box_size = 50.0

        os.chdir(self.test_dir)

        datafile_path = generate_ring_polymer_config(ring_length, box_size)

        # Basic checks
        self.assertTrue(os.path.exists(datafile_path))
        self.assertTrue(os.path.getsize(datafile_path) > 0)

    def test_generate_brush_polymer_config(self):
        # Test basic polymer generation
        backbone_length = 100
        grafting_density = 0.5
        side_chain_length = 20
        box_size = 100.0

        os.chdir(self.test_dir)

        datafile_path = generate_brush_polymer_config(backbone_length, grafting_density, side_chain_length, box_size)

        # Basic checks
        self.assertTrue(os.path.exists(datafile_path))
        self.assertTrue(os.path.getsize(datafile_path) > 0)

    def test_generate_brush_polymer_config_no_grafting(self):
        # Test with no side chains
        backbone_length = 100
        grafting_density = 0.0
        side_chain_length = 5
        box_size = 100.0

        os.chdir(self.test_dir)

        datafile_path = generate_brush_polymer_config(backbone_length, grafting_density, side_chain_length, box_size)

        # Basic checks
        self.assertTrue(os.path.exists(datafile_path))
        self.assertTrue(os.path.getsize(datafile_path) > 0)

    def test_generate_star_polymer_config(self):
        # Test star polymer generation
        arm_length = 10
        num_arms = 6
        box_size = 40.0

        os.chdir(self.test_dir)

        datafile_path = generate_star_polymer_config(arm_length, num_arms, box_size)

        # Basic checks
        self.assertTrue(os.path.exists(datafile_path))
        self.assertTrue(os.path.getsize(datafile_path) > 0)


if __name__ == "__main__":
    unittest.main()
