import unittest
import os
import json
from tools.config_gen.generate_polymer import generate_brush_polymer_config, generate_linear_polymer_config, generate_ring_polymer_config, generate_star_polymer_config, generate_dendrimer_config
from tools.config_gen.pack_solvent import pack_solvent


class TestPackSolvent(unittest.TestCase):

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
        #     if filename.startswith(('brush_polymer_', 'system_')):
        #         os.remove(os.path.join(self.test_dir, filename))

    def test_pack_solvent_linear_polymer(self):
        # Test packing solvent around a linear polymer
        chain_length = 30
        box_size = 20.0
        solvent_density = 0.4

        os.chdir(self.test_dir)

        # Generate linear polymer
        polymer_file = generate_linear_polymer_config(chain_length, box_size)

        # Pack solvent
        system_file = pack_solvent(polymer_file, solvent_density, box_size)

        # Basic checks
        self.assertTrue(os.path.exists(system_file))
        self.assertTrue(os.path.getsize(system_file) > 0)

        # Check filename
        expected_filename = "system_linear.data"
        self.assertEqual(os.path.basename(system_file), expected_filename)

        # Check metadata
        with open(system_file, 'r') as f:
            first_line = f.readline().strip()
            self.assertTrue(first_line.startswith("# Metadata:"))
            metadata_str = first_line[len("# Metadata:"):].strip()
            metadata = json.loads(metadata_str)
            self.assertEqual(metadata["type"], "linear")
            self.assertEqual(metadata["chain_length"], chain_length)
            self.assertEqual(metadata["solvent_density"], solvent_density)
            self.assertIn("num_solvent", metadata)

    def test_pack_solvent_ring_polymer(self):
        # Test packing solvent around a ring polymer
        ring_length = 60
        box_size = 20.0
        solvent_density = 0.4

        os.chdir(self.test_dir)

        # Generate ring polymer
        polymer_file = generate_ring_polymer_config(ring_length, box_size)

        # Pack solvent
        system_file = pack_solvent(polymer_file, solvent_density, box_size)

        # Basic checks
        self.assertTrue(os.path.exists(system_file))
        self.assertTrue(os.path.getsize(system_file) > 0)

        # Check filename
        expected_filename = "system_ring.data"
        self.assertEqual(os.path.basename(system_file), expected_filename)

        # Check metadata
        with open(system_file, 'r') as f:
            first_line = f.readline().strip()
            self.assertTrue(first_line.startswith("# Metadata:"))
            metadata_str = first_line[len("# Metadata:"):].strip()
            metadata = json.loads(metadata_str)
            self.assertEqual(metadata["type"], "ring")
            self.assertEqual(metadata["ring_length"], ring_length)
            self.assertEqual(metadata["solvent_density"], solvent_density)
            self.assertIn("num_solvent", metadata)

    def test_pack_solvent_brush_polymer_config(self):
        # Test packing solvent around a brush polymer
        backbone_length = 30
        grafting_density = 0.3
        side_chain_length = 5
        box_size = 20.0
        solvent_density = 0.4

        os.chdir(self.test_dir)

        # Generate polymer with side chains
        polymer_file = generate_brush_polymer_config(backbone_length, grafting_density, side_chain_length, box_size)

        # Pack solvent
        system_file = pack_solvent(polymer_file, solvent_density, box_size)

        # Basic checks
        self.assertTrue(os.path.exists(system_file))
        self.assertTrue(os.path.getsize(system_file) > 0)

        # Check filename
        expected_filename = "system_brush.data"
        self.assertEqual(os.path.basename(system_file), expected_filename)

        # Check metadata
        with open(system_file, 'r') as f:
            first_line = f.readline().strip()
            self.assertTrue(first_line.startswith("# Metadata:"))
            metadata_str = first_line[len("# Metadata:"):].strip()
            metadata = json.loads(metadata_str)
            self.assertEqual(metadata["type"], "brush")
            self.assertEqual(metadata["backbone_length"], backbone_length)
            self.assertEqual(metadata["grafting_density"], grafting_density)
            self.assertEqual(metadata["side_chain_length"], side_chain_length)
            self.assertEqual(metadata["solvent_density"], solvent_density)
            self.assertIn("num_solvent", metadata)

    def test_pack_solvent_star_polymer(self):
        # Test packing solvent around a star polymer
        arm_length = 10
        num_arms = 5
        box_size = 20.0
        solvent_density = 0.4

        os.chdir(self.test_dir)

        # Generate star polymer
        polymer_file = generate_star_polymer_config(arm_length, num_arms, box_size)

        # Pack solvent
        system_file = pack_solvent(polymer_file, solvent_density, box_size)

        # Basic checks
        self.assertTrue(os.path.exists(system_file))
        self.assertTrue(os.path.getsize(system_file) > 0)

        # Check filename
        expected_filename = "system_star.data"
        self.assertEqual(os.path.basename(system_file), expected_filename)

        # Check metadata
        with open(system_file, 'r') as f:
            first_line = f.readline().strip()
            self.assertTrue(first_line.startswith("# Metadata:"))
            metadata_str = first_line[len("# Metadata:"):].strip()
            metadata = json.loads(metadata_str)
            self.assertEqual(metadata["type"], "star")
            self.assertEqual(metadata["arm_length"], arm_length)
            self.assertEqual(metadata["num_arms"], num_arms)
            self.assertEqual(metadata["solvent_density"], solvent_density)
            self.assertIn("num_solvent", metadata)

    def test_pack_solvent_dendrimer_config(self):
        # Test packing solvent around a dendrimer
        generations = 3
        branching_factor = 3
        spacer = 3
        box_size = 20.0
        solvent_density = 0.4

        os.chdir(self.test_dir)

        # Generate dendrimer
        polymer_file = generate_dendrimer_config(generations, branching_factor, spacer, box_size)

        # Pack solvent
        system_file = pack_solvent(polymer_file, solvent_density, box_size)

        # Basic checks
        self.assertTrue(os.path.exists(system_file))
        self.assertTrue(os.path.getsize(system_file) > 0)

        # Check filename
        expected_filename = "system_dendrimer.data"
        self.assertEqual(os.path.basename(system_file), expected_filename)

        # Check metadata
        with open(system_file, 'r') as f:
            first_line = f.readline().strip()
            self.assertTrue(first_line.startswith("# Metadata:"))
            metadata_str = first_line[len("# Metadata:"):].strip()
            metadata = json.loads(metadata_str)
            self.assertEqual(metadata["type"], "dendrimer")
            self.assertEqual(metadata["generations"], generations)
            self.assertEqual(metadata["branching_factor"], branching_factor)
            self.assertEqual(metadata["spacer"], spacer)
            self.assertEqual(metadata["solvent_density"], solvent_density)
            self.assertIn("num_solvent", metadata)


if __name__ == "__main__":
    unittest.main()
