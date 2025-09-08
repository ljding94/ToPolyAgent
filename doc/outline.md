### Revised Project Development Plan

The core idea remains: a phased approach to build a multi-agent system for coarse-grained polymer simulations, with tools organized by function (config generation, simulation, analysis) and AI-friendly wrappers for CrewAI integration. The main revision is in the simulation tools: instead of generating a LAMMPS script dynamically, we'll use a fixed `simulation.lammps` script that reads the agent-generated datafile as input, with `run_lammps.py` handling execution. This reduces complexity in the simulation step and lets agents focus on passing the correct datafile path.

#### Revised Project Structure
The folder structure is updated to reflect the static LAMMPS script. Here's the new layout:

```
polymer_sim_project/
├── src/
│   ├── tools/                  # Core tools (pure Python functions/classes)
│   │   ├── config_gen/         # Tools for generating initial structures/configs
│   │   │   ├── __init__.py
│   │   │   ├── generate_polymer.py  # Functions for star, brush, etc.
│   │   │   └── pack_solvent.py      # Wraps Packmol
│   │   ├── simulation/         # Tools for running simulations
│   │   │   ├── __init__.py
│   │   │   ├── simulation.lammps    # Standard LAMMPS script with variable datafile
│   │   │   └── run_lammps.bash        # Executes LAMMPS with script
│   │   └── analysis/           # Tools for post-processing
│   │       ├── __init__.py
│   │       ├── calculate_rg.py
│   │       ├── calculate_pq.py
│   │       ├── calculate_kuhn.py
│   │       ├── calculate_diffusion.py
│   │       └── batch_analyze.py
│   ├── wrappers/               # AI-friendly wrappers (CrewAI Tool classes)
│   │   ├── __init__.py
│   │   ├── config_wrappers.py
│   │   ├── sim_wrappers.py
│   │   └── analysis_wrappers.py
│   ├── agents/                 # Agent definitions and crew setups
│   │   ├── __init__.py
│   │   ├── structure_agent.py  # build the polymer structure and place solvent, output datafile as system initial configuration
│   │   ├── simulation_agent.py # runs the lammps simulation using the static simulation.lammps script and the datafile from structure_agent
│   │   ├── analysis_agent.py   # run analysis tools on the dump file from simulation_agent, will only run analysis requested by user
│   │   ├── reporting_agent.py  # compile results into a report (Markdown)
│   │   └── main_crew.py
│   └── utils/                  # Shared utilities
│       ├── __init__.py
│       ├── constants.py        # Force field params, defaults
│       ├── file_handlers.py    # File I/O helpers
│       └── logging.py         # Logging setup
├── tests/                      # Unit/integration tests
│   ├── test_tools/
│   ├── test_wrappers/
│   └── test_agents/
├── data/                       # Sample datafiles, dumps (git-ignored)
├── docs/                       # Documentation
├── requirements.txt            # Dependencies
├── setup.py                    # Optional packaging
└── README.md                   # Project overview
```

**Key Change**:
- Replaced `generate_lammps_script.py` with `simulation.lammps`, a static LAMMPS script that uses a variable (e.g., `variable datafile string "path/to/datafile"`) to read the agent-generated datafile. This script will define standard settings (e.g., LJ pair style, FENE bonds, Langevin/Nose-Hoover thermostat, NVT/NPT ensemble) with customizable params via variables or command-line arguments.

#### Phase 1: Develop Core Tools
**Focus**: Build and test standalone tools, including the static LAMMPS script. Validate by running a manual simulation (e.g., brush polymer) using the tools.

**Estimated Time**: 1-2 weeks (10-20 hours/week).

**Sub-Tasks**:
1. **Setup**:
   - Install deps: `pip install crewai mdanalysis numpy scipy`. Install LAMMPS/Packmol binaries via conda or manual builds.
   - Write `README.md` with setup instructions.

2. **Config Generation Tools (src/tools/config_gen/)**:
   - `generate_polymer.py`: Generate configs for each topology (e.g., `generate_brush_polymer_config(backbone_length, grafting_density, side_chain_length, box_size)` → LAMMPS datafile).
     - Use NumPy for random walks, ensure no overlaps.
     - Test: Check datafile format, bead positions.
   - `pack_solvent.py`: Wrap Packmol to add solvent beads around polymer.
     - Inputs: Polymer datafile path, solvent density, box size.
     - Outputs: Full system datafile path.
     - Test: Validate solvent packing with Packmol's output.

3. **Simulation Tools (src/tools/simulation/)**:
   - `simulation.lammps`: Write a standard LAMMPS script with:
     - Variable for datafile: `variable datafile string "path/to/system.data"`.
     - Standard settings: LJ pair style, FENE bonds, thermostat (Langevin or Nose-Hoover), ensemble (NVT or NPT), dump for trajectories, thermo outputs (pressure, energy).
     - Example snippet:
       ```lammps
       # simulation.lammps
       units lj
       atom_style bond
       variable datafile string "system.data"
       read_data ${datafile}
       pair_style lj/cut 2.5
       pair_coeff 1 1 1.0 1.0
       bond_style fene
       bond_coeff 1 30.0 1.5 1.0 1.0
       fix 1 all nvt temp 1.0 1.0 0.1  # or npt for pressure control
       dump 1 all custom 1000 traj.dump id type x y z
       thermo 100
       run 100000
       ```
     - Make params like thermostat, ensemble, or run steps configurable via variables or command-line args (e.g., `lmp -var datafile system.data -var temp 1.0`).
     - Test: Run with a sample datafile, check dump/log outputs.
   - `run_lammps.py`: Execute LAMMPS with `simulation.lammps`.
     - Inputs: Datafile path, optional params (e.g., MPI cores, temp).
     - Outputs: Paths to dump file, log file, status.
     - Use subprocess: `subprocess.run(['lmp', '-in', 'simulation.lammps', '-var', 'datafile', datafile_path])`.
     - Test: Verify execution, file outputs.

4. **Analysis Tools (src/tools/analysis/)**:
   - `calculate_rg.py`: Compute radius of gyration (use MDAnalysis).
   - `calculate_pq.py`: Compute structure factor P(q) .
   - `calculate_kuhn.py`: Estimate Kuhn length from bond correlations.
   - `calculate_diffusion.py`: Compute diffusion constant from MSD.
   - `batch_analyze.py`: Run multiple analyses in one call.
   - Test: Use sample dump files to verify results.

5. **Utils (src/utils/)**:
   - `constants.py`: Define defaults (LJ params, box size, etc.).
   - `file_handlers.py`: Helpers for reading/writing datafiles.
   - `logging.py`: Setup logging for debugging.

**Deliverables**: Working tools, a `manual_sim.py` script to run a full simulation (e.g., brush polymer), and unit tests in `tests/test_tools/`.

#### Phase 2: Build Agent-Friendly Tools (Wrappers)
**Focus**: Wrap core tools into CrewAI-compatible Tool classes for LLM use.

**Estimated Time**: 1 week.

**Sub-Tasks**:
1. **Config Wrappers (src/wrappers/config_wrappers.py)**:
   - Wrap each config tool (e.g., BrushPolymerGeneratorTool).
   - Example:
     ```python
     from crewai import Tool
     from src.tools.config_gen.generate_polymer import generate_brush_polymer_config

     class BrushPolymerGeneratorTool(Tool):
         name = "generate_brush_polymer_config"
         description = """Generates LAMMPS datafile for a brush polymer.
         Inputs: backbone_length (int), grafting_density (float, 0-1), side_chain_length (int), box_size (float, default 50.0).
         Outputs: Path to datafile."""
         def _run(self, backbone_length: int, grafting_density: float, side_chain_length: int, box_size: float = 50.0) -> str:
             return generate_brush_polymer_config(backbone_length, grafting_density, side_chain_length, box_size)
     ```

2. **Simulation Wrappers (src/wrappers/sim_wrappers.py)**:
   - Wrap `run_lammps.py` (no wrapper needed for `simulation.lammps` since it's static).
   - Example:
     ```python
     class RunLammpsTool(Tool):
         name = "run_lammps"
         description = """Runs LAMMPS simulation with standard script.
         Inputs: datafile_path (str), optional: temp (float), steps (int).
         Outputs: Paths to dump and log files."""
         def _run(self, datafile_path: str, temp: float = 1.0, steps: int = 100000) -> dict:
             return run_lammps(datafile_path, temp, steps)
     ```

3. **Analysis Wrappers (src/wrappers/analysis_wrappers.py)**:
   - Wrap each analysis tool + batch analyzer.
   - Add a report generator tool (takes analysis results, outputs Markdown/JSON).

4. **Testing**:
   - Unit tests: Mock LLM inputs to wrappers.
   - Integration: Test a partial crew (e.g., config → sim).

**Deliverables**: Wrapped tools, testable via a simple CrewAI script.

#### Phase 3: Wrap Up Entire System
**Focus**: Build agents, integrate with CrewAI, test end-to-end pipeline.

**Estimated Time**: 1-2 weeks.

**Sub-Tasks**:
1. **Define Agents (src/agents/)**:
   - **StructureAgent**: Uses config wrappers to parse prompt, generate polymer, pack solvent.
   - **SimulationAgent**: Uses sim wrapper to run LAMMPS with `simulation.lammps`.
   - **AnalysisAgent**: Uses analysis wrappers to compute properties.
   - **ReportingAgent**: Compiles results into a report.
   - **MainCrew**: Orchestrates agents, passes file paths via shared memory.

2. **Integration**:
   - Create `main.py`: Entry point to take prompt, run crew, output report.
   - Ensure agents pass paths correctly (e.g., datafile from StructureAgent to SimulationAgent).

3. **Testing**:
   - Unit/Integration tests: Cover all tools, wrappers, agents.
   - End-to-End: Test with multiple polymers (brush, star, etc.).
   - Handle errors: Retry on Packmol failure, validate inputs.

4. **Polish**:
   - Add logging for agent actions.
   - Document: Update README with usage examples.
   - Optional: Dockerize for portability.

**Deliverables**: Fully automated system, demo script, comprehensive docs.

#### Notes on `simulation.lammps`
- **Why Static?**: A fixed script reduces agent complexity (no need to generate scripts dynamically). It can still be flexible with variables (e.g., `variable temp index 1.0`, `variable steps index 100000`).
- **Customization**: Allow agents to pass params via command-line vars (e.g., `lmp -var datafile system.data -var temp 1.2`). Hardcode defaults for LJ/FENE params, thermostat, etc., but make them overrideable if needed.
- **Future Flexibility**: If you later need dynamic scripts (e.g., for different ensembles), we can add a `generate_lammps_script.py` tool without major changes.

#### Next Steps
This revised plan keeps the simulation step lean while preserving modularity. To kick off Phase 1, I suggest starting with the brush polymer config tool (`generate_brush_polymer_config`). If you share specs (e.g., bead-spring model details, LJ params, or box size), I can propose a code outline. Alternatively, we could draft the `simulation.lammps` script first to define the simulation baseline. What do you think—where should we start?