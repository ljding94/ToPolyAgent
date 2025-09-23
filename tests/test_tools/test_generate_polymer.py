import unittest
import os
import json
from tools.config_gen.generate_polymer import generate_brush_polymer_config, generate_linear_polymer_config, generate_ring_polymer_config, generate_star_polymer_config, generate_dendrimer_config


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
        chain_length = 30
        box_size = 30.0

        os.chdir(self.test_dir)

        datafile_path = generate_linear_polymer_config(chain_length, box_size)

        # Basic checks
        self.assertTrue(os.path.exists(datafile_path))
        self.assertTrue(os.path.getsize(datafile_path) > 0)

        # Check filename
        expected_filename = "polymer_linear.data"
        self.assertEqual(os.path.basename(datafile_path), expected_filename)

        # Check metadata
        with open(datafile_path, 'r') as f:
            first_line = f.readline().strip()
            self.assertTrue(first_line.startswith("# Metadata:"))
            metadata_str = first_line[len("# Metadata:"):].strip()
            metadata = json.loads(metadata_str)
            self.assertEqual(metadata["type"], "linear")
            self.assertEqual(metadata["chain_length"], chain_length)
            self.assertEqual(metadata["box_size"], box_size)

    def test_generate_ring_polymer_config(self):
        # Test ring polymer generation
        ring_length = 40
        box_size = 30.0

        os.chdir(self.test_dir)

        datafile_path = generate_ring_polymer_config(ring_length, box_size)

        # Basic checks
        self.assertTrue(os.path.exists(datafile_path))
        self.assertTrue(os.path.getsize(datafile_path) > 0)

        # Check filename
        expected_filename = "polymer_ring.data"
        self.assertEqual(os.path.basename(datafile_path), expected_filename)

        # Check metadata
        with open(datafile_path, 'r') as f:
            first_line = f.readline().strip()
            self.assertTrue(first_line.startswith("# Metadata:"))
            metadata_str = first_line[len("# Metadata:"):].strip()
            metadata = json.loads(metadata_str)
            self.assertEqual(metadata["type"], "ring")
            self.assertEqual(metadata["ring_length"], ring_length)
            self.assertEqual(metadata["box_size"], box_size)

    def test_generate_brush_polymer_config(self):
        # Test basic polymer generation
        backbone_length = 40
        grafting_density = 0.1
        side_chain_length = 10
        box_size = 30.0

        os.chdir(self.test_dir)

        datafile_path = generate_brush_polymer_config(backbone_length, grafting_density, side_chain_length, box_size)

        # Basic checks
        self.assertTrue(os.path.exists(datafile_path))
        self.assertTrue(os.path.getsize(datafile_path) > 0)

        # Check filename
        expected_filename = "polymer_brush.data"
        self.assertEqual(os.path.basename(datafile_path), expected_filename)

        # Check metadata
        with open(datafile_path, 'r') as f:
            first_line = f.readline().strip()
            self.assertTrue(first_line.startswith("# Metadata:"))
            metadata_str = first_line[len("# Metadata:"):].strip()
            metadata = json.loads(metadata_str)
            self.assertEqual(metadata["type"], "brush")
            self.assertEqual(metadata["backbone_length"], backbone_length)
            self.assertEqual(metadata["grafting_density"], grafting_density)
            self.assertEqual(metadata["side_chain_length"], side_chain_length)
            self.assertEqual(metadata["box_size"], box_size)


    def test_generate_star_polymer_config(self):
        # Test star polymer generation
        arm_length = 10
        num_arms = 6
        box_size = 30.0

        os.chdir(self.test_dir)

        datafile_path = generate_star_polymer_config(arm_length, num_arms, box_size)

        # Basic checks
        self.assertTrue(os.path.exists(datafile_path))
        self.assertTrue(os.path.getsize(datafile_path) > 0)

        # Check filename
        expected_filename = "polymer_star.data"
        self.assertEqual(os.path.basename(datafile_path), expected_filename)

        # Check metadata
        with open(datafile_path, 'r') as f:
            first_line = f.readline().strip()
            self.assertTrue(first_line.startswith("# Metadata:"))
            metadata_str = first_line[len("# Metadata:"):].strip()
            metadata = json.loads(metadata_str)
            self.assertEqual(metadata["type"], "star")
            self.assertEqual(metadata["arm_length"], arm_length)
            self.assertEqual(metadata["num_arms"], num_arms)
            self.assertEqual(metadata["box_size"], box_size)

    def test_generate_dendrimer_config(self):
        # Test dendrimer generation
        generations = 3
        branching_factor = 3
        spacer = 3
        box_size = 30.0

        os.chdir(self.test_dir)

        datafile_path = generate_dendrimer_config(generations, branching_factor, spacer, box_size)

        # Basic checks
        self.assertTrue(os.path.exists(datafile_path))
        self.assertTrue(os.path.getsize(datafile_path) > 0)

        # Check filename
        expected_filename = "polymer_dendrimer.data"
        self.assertEqual(os.path.basename(datafile_path), expected_filename)

        # Check metadata
        with open(datafile_path, 'r') as f:
            first_line = f.readline().strip()
            self.assertTrue(first_line.startswith("# Metadata:"))
            metadata_str = first_line[len("# Metadata:"):].strip()
            metadata = json.loads(metadata_str)
            self.assertEqual(metadata["type"], "dendrimer")
            self.assertEqual(metadata["generations"], generations)
            self.assertEqual(metadata["branching_factor"], branching_factor)
            self.assertEqual(metadata["spacer"], spacer)
            self.assertEqual(metadata["box_size"], box_size)


if __name__ == "__main__":
    unittest.main()
