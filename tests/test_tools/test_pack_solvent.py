import unittest
import os
from tools.config_gen.generate_polymer import generate_brush_polymer_config, generate_linear_polymer_config, generate_ring_polymer_config, generate_star_polymer_config
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
        box_size = 30.0
        solvent_density = 0.1

        os.chdir(self.test_dir)

        # Generate linear polymer
        polymer_file = generate_linear_polymer_config(chain_length, box_size)

        # Pack solvent
        system_file = pack_solvent(polymer_file, solvent_density, box_size)

        # Basic checks
        self.assertTrue(os.path.exists(system_file))
        self.assertTrue(os.path.getsize(system_file) > 0)

    def test_pack_solvent_ring_polymer(self):
        # Test packing solvent around a ring polymer
        ring_length = 20
        box_size = 40.0
        solvent_density = 0.08

        os.chdir(self.test_dir)

        # Generate ring polymer
        polymer_file = generate_ring_polymer_config(ring_length, box_size)

        # Pack solvent
        system_file = pack_solvent(polymer_file, solvent_density, box_size)

        # Basic checks
        self.assertTrue(os.path.exists(system_file))
        self.assertTrue(os.path.getsize(system_file) > 0)

    def test_pack_solvent_brush_polymer_config(self):
        # Test packing solvent around a brush polymer
        backbone_length = 50
        grafting_density = 0.3
        side_chain_length = 10
        box_size = 50.0
        solvent_density = 0.05

        os.chdir(self.test_dir)

        # Generate polymer with side chains
        polymer_file = generate_brush_polymer_config(backbone_length, grafting_density, side_chain_length, box_size)

        # Pack solvent
        system_file = pack_solvent(polymer_file, solvent_density, box_size)

        # Basic checks
        self.assertTrue(os.path.exists(system_file))
        self.assertTrue(os.path.getsize(system_file) > 0)

    def test_pack_solvent_brush_polymer_config_no_grafting(self):
        # Test packing solvent around a brush polymer with no side chains
        backbone_length = 50
        grafting_density = 0.0
        side_chain_length = 3
        box_size = 20.0
        solvent_density = 0.05

        os.chdir(self.test_dir)

        # Generate polymer
        polymer_file = generate_brush_polymer_config(backbone_length, grafting_density, side_chain_length, box_size)

        # Pack solvent
        system_file = pack_solvent(polymer_file, solvent_density, box_size)

        # Basic checks
        self.assertTrue(os.path.exists(system_file))
        self.assertTrue(os.path.getsize(system_file) > 0)

    def test_pack_solvent_star_polymer(self):
        # Test packing solvent around a star polymer
        arm_length = 10
        num_arms = 6
        box_size = 40.0
        solvent_density = 0.06

        os.chdir(self.test_dir)

        # Generate star polymer
        polymer_file = generate_star_polymer_config(arm_length, num_arms, box_size)

        # Pack solvent
        system_file = pack_solvent(polymer_file, solvent_density, box_size)

        # Basic checks
        self.assertTrue(os.path.exists(system_file))
        self.assertTrue(os.path.getsize(system_file) > 0)


if __name__ == "__main__":
    unittest.main()
