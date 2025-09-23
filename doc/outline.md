# ToPolyAgent Project Development Plan

This document outlines the development plan for ToPolyAgent, a multi-agent system using CrewAI to automate coarse-grained polymer simulations with LAMMPS. The system supports topological polymers (linear, ring, brush, star) and automates configuration generation, solvent packing, simulation, and analysis. The pipeline is driven by user prompts specifying polymer type, parameters, solvent volume fraction, thermostat, interaction parameters, and simulation duration. A static LAMMPS script (`simulation.lammps`) ensures simplicity, with variables for flexibility.

## Project Structure

```
ToPolyAgent/
├── README.md                   # Project overview and setup instructions
├── requirements.txt            # Python dependencies
├── setup.py                    # Package setup script
├── TODOS.md                    # Task list
├── data/
│   └── test/
├── doc/
│   ├── note.md                # Project notes and ideas
│   └── outline.md             # This development plan
├── src/
│   ├── __init__.py            # Package init
│   ├── agents/
│   │   └── __init__.py
│   ├── tools/
│   │   ├── analysis/
│   │   │   ├── __init__.py
│   │   │   ├── calculate_conformation.py  # Rg, persistence length, diffusion
│   │   │   ├── calculate_pq.py           # Structure factor P(q)
│   │   │   ├── read_config.py            # Reads LAMMPS data/trajectories
│   │   │   └── run_analysis.py           # Comprehensive analysis
│   │   ├── config_gen/
│   │   │   ├── __init__.py
│   │   │   ├── generate_polymer.py       # Generates polymer configs
│   │   │   ├── pack_solvent.py           # Adds solvent beads
│   │   │   ├── plot_config.py            # Plots configurations
│   │   │   ├── polymer_brush.data        # Sample polymer data
│   │   │   ├── run_config_generation.py  # Runs config generation
│   │   │   └── system_brush.data         # Sample system data
│   │   └── simulation/
│   │       ├── __init__.py
│   │       ├── log.lammps                # LAMMPS log template
│   │       ├── run_lammps.py            # Executes LAMMPS
│   │       └── simulation.lammps         # Static LAMMPS script
│   ├── topolyagent.egg-info/
│   │   ├── dependency_links.txt
│   │   ├── PKG-INFO
│   │   ├── requires.txt
│   │   ├── SOURCES.txt
│   │   └── top_level.txt
│   └── wrappers/
│       ├── __init__.py
│       └── analysis_wrappers.py          # Placeholder for CrewAI wrappers
├── tests/
│   ├── __init__.py
│   ├── test_agents/
│   │   └── __init__.py
│   ├── test_tools/
│   │   ├── __init__.py
│   │   ├── log.lammps                    # Test log template
│   │   ├── test_analysis_tools.py        # Tests analysis tools
│   │   ├── test_full_workflow.py         # Tests full workflow
│   │   ├── test_generate_polymer.py      # Tests polymer generation
│   │   ├── test_pack_solvent.py          # Tests solvent packing
│   │   ├── test_run_lammps.py           # Tests LAMMPS execution
│   │   └── brush_workflow_50_0.5_10_20250916_163418/
│   │       ├── analysis_results.json
│   │       ├── conformation_analysis.png
│   │       ├── log.lammps
│   │       ├── polymer_brush.data
│   │       ├── polymer_brush.png
│   │       ├── system_brush.data
│   │       ├── system_brush.png
│   │       └── ...                       # Additional output files
│   └── test_wrappers/
│       └── __init__.py
```

## Pipeline Steps
The system follows a three-step pipeline, orchestrated by CrewAI agents, with inputs from user prompts and outputs as files or results.

1. **Generate Polymer Structure and Add Solvent**
   - **Inputs**:
     - Polymer type: linear, ring, brush, star.
     - Type-specific parameters:
       - Linear: chain_length (int, e.g., 30).
       - Ring: ring_length (int, e.g., 20).
       - Brush: backbone_length (int, e.g., 50), grafting_density (float 0-1, e.g., 0.3), side_chain_length (int, e.g., 10).
       - Star: arm_length (int, e.g., 8), num_arms (int, e.g., 4).
     - solvent_volume_fraction (float 0-1, e.g., 0.1).
     - box_size (float, default 50.0).
   - **Process**: Generate polymer using `generate_{type}_polymer_config` (Gaussian chains, bond length 1.0, angle constraints). Add solvent with `pack_solvent`, converting volume fraction to density (density = volume_fraction, assuming unit mass/volume in LJ units; num_solvent = int(volume_fraction * box_size^3)).
   - **Outputs**: LAMMPS system datafile (e.g., system_brush_polymer_50_0.3_10.data).
   - **Agent**: StructureAgent.
   - **Tools**: LinearPolymerGeneratorTool, RingPolymerGeneratorTool, BrushPolymerGeneratorTool, StarPolymerGeneratorTool, PackSolventTool.

2. **Run Simulation Using LAMMPS**
   - **Inputs**:
     - System datafile path (from Step 1).
     - Thermostat: langevin, nose-hoover.
     - Ensemble: nvt, npt.
     - Interaction parameters (LJ epsilon for pair_coeff):
       - pp: polymer-polymer (float, e.g., 1.0).
       - ss: solvent-solvent (float, e.g., 0.5).
       - sp: solvent-polymer (float, e.g., 0.8).
     - run_steps (int, e.g., 100000).
   - **Process**: Execute `simulation.lammps` via `run_lammps.py`, passing variables (e.g., -var eps_pp ${pp}, -var thermostat langevin).
   - **Outputs**: Dump files (coord/dump.*.txt), final configuration (final_state.data).
   - **Agent**: SimulationAgent.
   - **Tools**: RunLammpsTool.

3. **Analyze the Data**
   - **Inputs**:
     - System datafile (for bonds, atom types).
     - Dump files pattern (e.g., coord/dump.*.txt).
   - **Process**: Use `read_simulation_data` to load data/trajectories, `run_analysis.py` to compute radius of gyration (Rg), persistence length (lp), diffusion coefficient (D), structure factor (P(q)).
   - **Outputs**: JSON/dict with results (e.g., {"mean_rg": float, "mean_lp": float, "diffusion_coeff": float, "scattering_pq": array}), optional plots (conformation_analysis.png).
   - **Agent**: AnalysisAgent.
   - **Tools**: RadiusOfGyrationTool, PersistenceLengthTool, DiffusionCoefficientTool, PairCorrelationTool, ComprehensiveAnalysisTool.

## Development Phases

### Phase 1: Develop Core Tools
**Objective**: Build and test standalone tools for configuration, simulation, and analysis.

**Status**: ~90% complete (see doc/progress.md).

**Sub-Tasks**:
1. **Setup**:
   - Dependencies: crewai, mdanalysis, numpy, scipy (requirements.txt). LAMMPS/Packmol via conda.
   - Update README.md with setup instructions, usage examples.

2. **Config Generation Tools (src/tools/config_gen/)**:
   - `generate_polymer.py`: Implements linear, ring, brush, star polymers. Inputs: type-specific params, box_size. Outputs: LAMMPS datafile.
   - `pack_solvent.py`: Adds solvent beads, avoiding overlaps (min distance 1.5). Update to accept solvent_volume_fraction (density = volume_fraction).
   - Tests: tests/test_tools/test_generate_polymer.py, test_pack_solvent.py. Verify atom counts, file formats.

3. **Simulation Tools (src/tools/simulation/)**:
   - `simulation.lammps`: Static script with variables (datafile_path, dump_path, pdump=1000, prun=50000, pequi=5000, eps_pp, eps_ss, eps_sp, thermostat, ensemble, time_step=0.01). Uses LJ pair style (cutoff 2.5), FENE bonds, cosine angles. Example:
     ```lammps
     units lj
     atom_style full
     read_data ${datafile_path}
     pair_style lj/cut 2.5
     pair_coeff 1 1 ${eps_pp} 1.0 2.5  # Polymer-polymer
     pair_coeff 3 3 ${eps_ss} 1.0 2.5  # Solvent-solvent
     pair_coeff 1 3 ${eps_sp} 1.0 2.5  # Solvent-polymer
     if "${thermostat} == langevin" then "fix 1 all langevin 1.0 1.0 1.0 ${random}" else "fix 1 all npt temp 1.0 1.0 0.1 ..."
     dump 2 all custom ${pdump} ${dump_path}/coord/dump.*.txt id type x y z xu yu zu
     run ${prun}
     write_data ${dump_path}/final_state.data
     ```
   - `run_lammps.py`: Executes LAMMPS with variables for thermostat, ensemble, interactions, steps.
   - Tests: tests/test_tools/test_run_lammps.py. Verify dump/log outputs.

4. **Analysis Tools (src/tools/analysis/)**:
   - `calculate_conformation.py`: Computes Rg, lp, diffusion coefficient.
   - `calculate_pq.py`: Computes P(q) with q based on Rg.
   - `read_config.py`: Reads datafiles/trajectories, filters polymer atoms.
   - `run_analysis.py`: Runs all analyses, saves JSON, generates plots.
   - Tests: tests/test_tools/test_analysis_tools.py. Validate metrics for sample systems.

5. **Utils (src/utils/)**:
   - Planned: constants.py (LJ/FENE params), file_handlers.py (I/O), logging.py (debug logs).

**Deliverables**: Working tools, manual demo script (e.g., demo_pipeline.py for brush polymer), tests for all steps.

### Phase 2: Build Agent-Friendly Wrappers
**Objective**: Wrap tools into CrewAI Tool classes for LLM integration.

**Status**: Not started.

**Sub-Tasks**:
1. **Config Wrappers (src/wrappers/config_wrappers.py)**:
   - Tools: PolymerConfigTool (routes to generate_{type}_polymer_config), PackSolventTool.
   - Example:
     ```python
     from crewai import Tool
     from src.tools.config_gen.generate_polymer import generate_brush_polymer_config

     class PolymerConfigTool(Tool):
         name = "GeneratePolymerConfig"
         description = """Generates LAMMPS datafile for a polymer. Inputs: type (str: linear/ring/brush/star), params (dict: e.g., {'backbone_length': int, 'grafting_density': float, 'side_chain_length': int} for brush), box_size (float, default 50.0). Outputs: Datafile path."""
         def _run(self, type: str, params: dict, box_size: float = 50.0) -> str:
             if type == "brush":
                 return generate_brush_polymer_config(params["backbone_length"], params["grafting_density"], params["side_chain_length"], box_size)
             # Add other types
     ```

2. **Simulation Wrappers (src/wrappers/sim_wrappers.py)**:
   - Tool: RunLammpsTool. Inputs: datafile_path, dump_path, thermostat, ensemble, interaction_params, run_steps.
   - Example:
     ```python
     class RunLammpsTool(Tool):
         name = "RunLammps"
         description = """Runs LAMMPS simulation. Inputs: datafile_path (str), dump_path (str), thermostat (str: langevin/nose-hoover), ensemble (str: nvt/npt), interaction_params (dict: {'pp': float, 'ss': float, 'sp': float}), run_steps (int). Outputs: Dict with dump and final config paths."""
         def _run(self, datafile_path: str, dump_path: str, thermostat: str = "langevin", ensemble: str = "nvt", interaction_params: dict = {"pp": 1.0, "ss": 1.0, "sp": 1.0}, run_steps: int = 50000) -> dict:
             return run_lammps(dump_path, datafile_path, thermostat, ensemble, interaction_params, run_steps)
     ```

3. **Analysis Wrappers (src/wrappers/analysis_wrappers.py)**:
   - Tools: RadiusOfGyrationTool, PersistenceLengthTool, DiffusionCoefficientTool, PairCorrelationTool, ComprehensiveAnalysisTool.
   - ComprehensiveAnalysisTool runs run_analysis.py, returns JSON and plot paths.

4. **Testing**:
   - Unit tests: Mock LLM inputs (tests/test_wrappers/).
   - Integration: Test config → sim → analysis with sample prompt.

**Deliverables**: Wrapped tools, testable via CrewAI script.

### Phase 3: Wrap Up Entire System
**Objective**: Build agents, integrate with CrewAI, test end-to-end pipeline.

**Status**: Not started.

**Sub-Tasks**:
1. **Define Agents (src/agents/)**:
   - **StructureAgent**: Parses prompt for type/params, generates polymer + solvent.
   - **SimulationAgent**: Runs LAMMPS with user-specified settings.
   - **AnalysisAgent**: Computes requested metrics (Rg, lp, D, Pq).
   - **ReportingAgent**: Compiles results into Markdown with tables/plots.
   - **MainCrew**: Sequential tasks, passes file paths via shared memory.

2. **Integration**:
   - `main.py`: Parses prompt (e.g., "Simulate brush polymer, backbone 50, grafting 0.3, side 10, solvent 0.1, langevin, pp=1.0 ss=0.5 sp=0.8, steps=100k, analyze all"), runs crew, outputs report.
   - Prompt parsing: Use regex/LLM to extract type, params, settings.

3. **Testing**:
   - Unit tests: Each agent (tests/test_agents/).
   - End-to-End: Test prompts for multiple topologies, validate outputs.
   - Error handling: Retry on failures (e.g., solvent packing), validate inputs (e.g., grafting_density ≤ 1).

4. **Polish**:
   - Logging: Add src/utils/logging.py for agent/tool debugging.
   - Docs: Update README with example prompts.
   - Optional: Dockerize for portability.

**Deliverables**: Automated system, demo script, comprehensive docs.

## Notes
- **Configurability**: simulation.lammps uses variables for thermostat (langevin/nose-hoover), ensemble (nvt/npt), interactions (eps_pp, eps_ss, eps_sp), steps (prun). Future: Dynamic script generation if needed.
- **Extensibility**: Per doc/note.md, consider SearchAgent (web search for params) or HumanTool for clarifications.
- **Performance**: Optimize pack_solvent.py for high volume fractions (e.g., KD-tree). Test large systems (backbone=1000).
- **Validation**: Compare Rg (~sqrt(N/6) for ideal chains), diffusion to theoretical values. Validate P(q) with known systems.
- **ACS 2026**: Demo by Sept 29, 2025, emphasizing AI-driven topology exploration.

## Next Steps
- Update pack_solvent.py for volume_fraction, simulation.lammps/run_lammps.py for thermostat/ensemble/interactions.
- Implement Phase 2 wrappers, starting with config_wrappers.py.
- Test pipeline manually with demo_pipeline.py (config → sim → analysis).
- Begin Phase 3 with StructureAgent, test with simple prompt.