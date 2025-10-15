# Revised ToPolyAgent Project Development Plan (Version 2.0 - September 23, 2025)

This revised plan updates the original `outline.md` to reflect the current state post-Phase 1 completion. Key changes include:
- **Actual Directory Structure**: Incorporated the implemented structure from the project (e.g., added `plot_config.py`, `run_config_generation.py`, and specific test artifacts like workflow directories). Removed unimplemented items (e.g., `setup.py`, `TODOS.md`, `data/test/`) and adjusted for what's built (e.g., no `utils/` yet, wrappers are placeholders).
- **Phase 1 Status**: Marked as complete, with summaries of achievements, including support for 5 topologies (linear, ring, brush, star, dendrimer), solvent packing, static LAMMPS simulation, and analysis tools. Added notes on tested features (e.g., thermostats: Langevin/Nose-Hoover; ensembles: NVT/NPT implied via fixes).
- **Refinements for Phases 2-3**: Refreshed details based on current tools (e.g., wrappers will directly reference existing functions like `generate_brush_polymer_config`). Emphasized integration with CrewAI, prompt parsing, and extensibility (e.g., from `note.md`: agents for varying topologies/solvents/thermostats, potential human/search tools).
- **General Updates**: Static `simulation.lammps` confirmed as effective; added error handling and logging recommendations. Ensured plan aligns with bead-spring model (LJ/FENE, single-bead solvent). Prioritized modularity for non-trivial features (e.g., multi-agent optimization of params).

The system automates coarse-grained polymer simulations via CrewAI agents, driven by user prompts (e.g., "Simulate star polymer with 4 arms of length 10, solvent fraction 0.2, Langevin thermostat, analyze Rg and P(q)"). Pipeline: config generation → simulation → analysis → report.

## Project Structure

```
ljding94-topolyagent/
├── README.md                   # Project overview, setup, usage
├── requirements.txt            # Dependencies (crewai, mdanalysis, numpy, etc.)
├── doc/
│   ├── note.md                 # Ideas, references
│   └── outline.md              # This plan
├── src/
│   ├── __init__.py
│   ├── agents/                 # (Phase 3: Agent defs, to be added)
│   │   └── __init__.py
│   ├── tools/
│   │   ├── analysis/
│   │   │   ├── __init__.py
│   │   │   ├── calculate_conformation.py  # Rg, persistence length, diffusion
│   │   │   ├── calculate_pq.py           # Structure factor P(q)
│   │   │   ├── read_config.py            # Reads data/trajectories
│   │   │   └── run_analysis.py           # Batch analysis + JSON/plots
│   │   ├── config_gen/
│   │   │   ├── __init__.py
│   │   │   ├── generate_polymer.py       # Topology-specific generators
│   │   │   ├── pack_solvent.py           # Adds solvent via Packmol
│   │   │   ├── plot_config.py            # Plots configs (PyVista/Matplotlib)
│   │   │   └── run_config_generation.py  # Orchestrates config gen (optional)
│   │   └── simulation/
│   │       ├── __init__.py
│   │       ├── log.lammps                # Template/sample log
│   │       ├── run_lammps.py             # Executes LAMMPS
│   │       └── simulation.lammps         # Static script with variables
│   └── wrappers/                # (Phase 2: CrewAI Tool wrappers)
│       ├── __init__.py
│       └── analysis_wrappers.py  # Placeholder; expand for all tools
└── tests/
    ├── __init__.py
    ├── test_agents/            # (Phase 3: To be added)
    │   └── __init__.py
    ├── test_tools/
    │   ├── __init__.py
    │   ├── log.lammps          # Test log
    │   ├── test_analysis_tools.py
    │   ├── test_full_workflow.py  # End-to-end for all topologies
    │   ├── test_generate_polymer.py
    │   ├── test_pack_solvent.py
    │   └── test_run_lammps.py
    └── test_wrappers/          # (Phase 2: To be added)
        └── __init__.py
```

**Notes on Structure**:
- Data outputs (e.g., `data/test/`, workflow dirs like `brush_workflow_50_0.5_10_20250916_163418/`) are git-ignored.
- No `utils/` yet; add in Phase 2 if needed (e.g., `logging.py`, `constants.py` for defaults like bond_length=1.0, LJ_cutoff=2.5).
- Future: Add `main.py` for entry point, `demo_pipeline.py` for manual runs.

## Pipeline Steps
Unchanged from original, but now validated:
1. **Generate Polymer + Solvent**: Inputs: Type (str), params (dict), solvent_density (float 0-1), box_size (float=50.0). Outputs: System datafile (e.g., `system_brush.data`) with embedded JSON metadata.
2. **Run Simulation**: Inputs: Datafile path, thermostat (str: "langevin"/"nosehoover"), interaction_params (dict: {"pp": float, "ss": float, "sp": float}), run_steps (int). Uses static `simulation.lammps` (LJ pair, FENE bonds, dumps every 1000 steps). Outputs: Dump files (`coord/dump.*.txt`), log, final_state.data.
3. **Analyze Data**: Inputs: Datafile, dump pattern. Computes Rg, persistence length (lp), diffusion (D), P(q). Outputs: `analysis_results.json`, plots (e.g., `conformation_analysis.png`).

## Development Phases

### Phase 1: Develop Core Tools (Completed)
**Achievements**:
- Setup: Dependencies installed; LAMMPS/Packmol via conda.
- Config Gen: All 5 topologies implemented; solvent packing handles density→num_solvent calc (e.g., int(density * box_size**3)); plots with/without solvent.
- Simulation: Static script flexible via vars (e.g., `-var eps_pp 0.5 -var thermostat langevin`); supports NVT/NPT via fixes; tested with 10k-20k steps.
- Analysis: Full metrics; handles polymer-specific calcs (e.g., chain-wise lp for brush/star); JSON outputs with means/errors.
- Utils: Implicit in tools (e.g., file I/O in `read_config.py`).
- Testing: Unit tests pass; `test_full_workflow.py` runs for all topologies, producing dirs with datafiles, dumps, JSON, plots.
- Validation: Matches theory (e.g., good solvent: sp>0.5 collapses less); no major bugs.

**Deliverables**: All tools functional; manual workflow via `test_full_workflow.py`.

### Phase 2: Build Agent-Friendly Wrappers
**Objective**: Wrap core tools as CrewAI `Tool` classes for LLM integration. Focus on type-safe inputs/outputs; handle errors gracefully. Add utils for shared needs.


**Sub-Tasks**:
1. **Setup Utils (src/utils/)**:
   - Add `__init__.py`, `constants.py` (defaults: box_size=50.0, bond_length=1.0, interaction_params={"pp":0.5, "ss":0.5, "sp":0.5}), `logging.py` (setup logger for tool debug), `file_handlers.py` (metadata parsers, path validators).

2. **Config Wrappers (src/wrappers/config_wrappers.py)**:
   - Create one Tool per major function, routing via params.
   - Example for Polymer Generation:
     ```python
     from crewai import Tool
     from src.tools.config_gen.generate_polymer import generate_linear_polymer_config, generate_ring_polymer_config, ...  # Import all topology funcs

     class PolymerGeneratorTool(Tool):
         name = "GeneratePolymerConfig"
         description = """Generates LAMMPS datafile for specified polymer topology.
         Inputs: polymer_type (str: 'linear'/'ring'/'brush'/'star'/'dendrimer'), params (dict: topology-specific, e.g., {'chain_length':30} for linear), box_size (float=50.0).
         Outputs: Path to polymer datafile (str)."""
         def _run(self, polymer_type: str, params: dict, box_size: float = 50.0) -> str:
             if polymer_type == "linear":
                 return generate_linear_polymer_config(params.get("chain_length", 30), box_size)
             elif polymer_type == "brush":
                 return generate_brush_polymer_config(params.get("backbone_length", 50), params.get("grafting_density", 0.3), params.get("side_chain_length", 10), box_size)
             # Add cases for ring, star, dendrimer; raise ValueError on invalid type
             logging.info(f"Generated {polymer_type} config at {path}")
             return path
     ```
   - `PackSolventTool`: Inputs: polymer_file (str), solvent_density (float), box_size (float). Outputs: System datafile path.
   - `PlotConfigTool`: Optional, for visuals.

3. **Simulation Wrappers (src/wrappers/sim_wrappers.py)**:
   - `RunLammpsTool`:
     ```python
     class RunLammpsTool(Tool):
         name = "RunLammpsSimulation"
         description = """Executes LAMMPS using static script.
         Inputs: datafile_path (str), thermostat (str='langevin'/'nosehoover'), interaction_params (dict={'pp':0.5, 'ss':0.5, 'sp':0.5}), run_steps (int=100000).
         Outputs: Dict with 'dump_files' (str), 'final_config' (str), 'log' (str)."""
         def _run(self, datafile_path: str, thermostat: str = "langevin", interaction_params: dict = None, run_steps: int = 100000) -> dict:
             if interaction_params is None:
                 interaction_params = constants.DEFAULT_INTERACTIONS
             dump_path = datafile_path.replace(".data", "")
             try:
                 results = run_lammps(dump_path, datafile_path, thermostat, interaction_params, run_steps)
                 return results
             except Exception as e:
                 logging.error(f"Simulation failed: {e}")
                 return {"error": str(e)}
     ```

4. **Analysis Wrappers (src/wrappers/analysis_wrappers.py)**:
   - Expand existing placeholder.
   - `ComprehensiveAnalysisTool`: Inputs: datafile_path (str), dump_pattern (str). Outputs: Path to analysis_results.json.
   - Individual tools: e.g., `RadiusOfGyrationTool` for targeted analysis.

5. **Testing (tests/test_wrappers/)**:
   - Add tests: Mock _run calls, assert paths exist, handle errors (e.g., invalid thermostat → log warning).

**Deliverables**: Wrapped tools; simple script to test via CrewAI (e.g., create agent, assign tools, run task like "Generate brush polymer").

### Phase 3: Wrap Up Entire System
**Objective**: Define agents, orchestrate via CrewAI, enable prompt-driven runs. Add reporting and extensibility.

**Sub-Tasks**:
1. **Define Agents (src/agents/)**:
   - **StructureAgent**: Goal: "Parse prompt for polymer type/params/solvent, generate config." Tools: PolymerGeneratorTool, PackSolventTool. Outputs: System datafile path.
   - **SimulationAgent**: Goal: "Run simulation with specified thermostat/params/steps." Tools: RunLammpsTool. Outputs: Dump/log paths.
   - **AnalysisAgent**: Goal: "Analyze dumps for requested metrics (e.g., all or specific like Rg)." Tools: Analysis wrappers. Dynamically select based on prompt (e.g., if "analyze Rg and P(q)").
   - **ReportingAgent**: Goal: "Compile results into Markdown report with tables/plots." Tools: Custom ReportTool (reads JSON, embeds images).

2. **Integration (src/agents/main_crew.py)**:
   - Use CrewAI's Crew: Sequential tasks, shared memory for paths (e.g., datafile from Structure to Simulation).
   - Prompt Parsing: In StructureAgent, use LLM or regex to extract (e.g., "brush polymer, backbone 50, grafting 0.3, side 10, solvent 0.1, langevin, steps=100k, analyze all").
   - Entry Point: Add `main.py` – Inputs: User prompt. Runs crew, outputs report.

3. **Extensibility**:
   - Per `note.md`: Add optional agents (e.g., ParameterOptimizerAgent: Varies solvent density/thermostat; uses web_search tool for refs). HumanTool for clarifications.
   - Error Handling: Agents retry on failures (e.g., Packmol overlap → reduce density).

4. **Testing (tests/test_agents/)**:
   - Unit: Mock tools for each agent.
   - End-to-End: Test prompts for topologies; validate reports.

5. **Polish**:
   - Logging: Integrate across agents.
   - Optional: Docker for portability.


## Next Steps
- Start Phase 2: Implement `config_wrappers.py` first (focus on brush for quick win).
- By Sept 25: Phase 2 done; test partial crew.

This fresh outline positions us for rapid Phase 2/3 implementation. Thoughts on priorities?