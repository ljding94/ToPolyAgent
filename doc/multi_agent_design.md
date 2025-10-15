
# Multi-Agent Design for ToPolyAgent

## Outline

Our multi-agent system should have two modes: 1) autonomous mode and 2) interactive modes.

The 1) autonomous mode will only have crew of two agents: a workflow agent that uses the workflow tool to complete the entire config generation, simulation, analysis pipeline and a writer_agent that writes a report.

The 2) interactive mode will have 3 agents: a) config generation agent that generates config, b) simulation agent that runs the simulation and analyzes the data, c) writer agent that writes the report, and each agent will ask for user feedback and confirmation before passing the job to the next agent.

---

## Overview

ToPolyAgent implements a multi-agent system for polymer simulation workflows with two operational modes: **Autonomous Mode** and **Interactive Mode**. The system leverages CrewAI to coordinate specialized agents that handle different aspects of polymer simulation, from configuration generation to analysis and reporting.

## System Architecture

### Core Components

1. **Wrapper Tools**: Specialized tools that encapsulate simulation functionality
   - `FullWorkflowTool`: Complete end-to-end simulation pipeline
   - Individual tools for config generation, simulation, and analysis

2. **Agents**: Specialized AI agents with domain expertise
   - Workflow Agent: Orchestrates complete simulation workflows
   - Config Generation Agent: Handles polymer configuration creation
   - Simulation Agent: Manages molecular dynamics simulations
   - Writer Agent: Generates comprehensive reports

3. **LLM Integration**: OpenRouter API with configurable models for natural language processing

## Mode 1: Autonomous Mode

### Architecture
```
Workflow Agent → Writer Agent
```

### Agent Roles

#### Workflow Agent
- **Purpose**: Executes the complete simulation pipeline using the FullWorkflowTool, handling natural language input directly
- **Input**: User prompt (e.g., "Simulate a linear polymer with 50 monomers")
- **Output**: Complete simulation results including all output paths
- **Tools**: FullWorkflowTool
- **Responsibilities**:
  - Parse natural language prompt into simulation parameters
  - Polymer configuration generation
  - Solvent packing
  - Molecular dynamics simulation
  - Trajectory analysis
  - Result visualization
- **Output**: Complete simulation results including all output paths
- **Tools**: FullWorkflowTool
- **Responsibilities**:
  - Polymer configuration generation
  - Solvent packing
  - Molecular dynamics simulation
  - Trajectory analysis
  - Result visualization

#### Writer Agent
- **Purpose**: Generates comprehensive Markdown reports from simulation results
- **Input**: Simulation results and analysis data
- **Output**: Formatted report with figures, metrics, and conclusions
- **Tools**: Report generation utilities

### Workflow Execution
1. User provides natural language description
2. Workflow Agent parses prompt and runs complete pipeline
3. Writer Agent generates final report
4. All results saved to timestamped directory

## Mode 2: Interactive Mode

### Architecture
```
Config Agent → [User Feedback] → Simulation Agent → [User Feedback] → Writer Agent
```

### Agent Roles

#### Config Generation Agent
- **Purpose**: Creates polymer configurations with user validation, handling natural language input directly
- **Input**: User prompt with polymer specifications
- **Output**: Polymer and system configuration files with plots
- **Tools**: Individual config generation tools (GenerateLinearPolymerTool, etc.)
- **Responsibilities**:
  - Parse natural language prompt into polymer parameters
  - Generate polymer structure
  - Pack solvent molecules
  - Create initial visualizations
  - Present results for user approval
- **Output**: Polymer and system configuration files with plots
- **Tools**: Individual config generation tools (GenerateLinearPolymerTool, etc.)
- **Responsibilities**:
  - Generate polymer structure
  - Pack solvent molecules
  - Create initial visualizations
  - Present results for user approval

#### Simulation Agent
- **Purpose**: Runs molecular dynamics simulations with user oversight
- **Input**: System configuration from Config Agent
- **Output**: Simulation trajectories and final configurations
- **Tools**: RunLammpsTool for simulation execution
- **Responsibilities**:
  - Execute LAMMPS simulations
  - Monitor convergence
  - Generate trajectory data
  - Present results for user feedback

#### Writer Agent
- **Same as Autonomous Mode**: Generates comprehensive reports

### Interactive Workflow Execution
1. User provides natural language description
2. **Config Generation Phase**:
   - Config Agent parses prompt and generates polymer and system
   - User reviews results and provides feedback
   - Parameters modified if requested, step repeated
3. **Simulation Phase**:
   - Simulation Agent runs MD simulation
   - User reviews trajectory and final state
   - Parameters modified if requested, step repeated
4. **Analysis & Reporting Phase**:
   - Analysis performed automatically
   - Writer Agent generates report
   - User can request report regeneration

## Parameter Structures

### Config Dictionary Format
```python
config = {
    "polymer_type": "linear",  # "linear", "ring", "brush", "star", "dendrimer"
    "polymer_params": {
        # Type-specific parameters
        "chain_length": 50  # for linear/ring
        # OR
        "backbone_length": 20,      # for brush
        "grafting_density": 0.5,
        "side_chain_length": 10
        # OR
        "arm_length": 25,    # for star
        "num_arms": 4
        # OR
        "generation": 3,     # for dendrimer
        "branching_factor": 2
    },
    "box_size": 20.0,              # nm
    "solvent_density": 0.3,        # relative
    "run_steps": 10000,            # simulation steps
    "thermostat": "langevin",      # thermostat type
    "interaction_params": {        # polymer-polymer, solvent-solvent, polymer-solvent
        "pp": 0.5,
        "ss": 0.5,
        "sp": 0.5
    }
}
```

## User Feedback System

### Feedback Interpretation
The interactive mode uses LLM-powered feedback interpretation to understand user requests:

- **Proceed**: Continue to next step ("yes", "looks good", "continue")
- **Redo**: Repeat current step ("redo", "try again", "run it again")
- **Modify**: Change parameters and redo ("make it longer", "increase temperature")
- **Stop**: End workflow ("stop", "quit", "cancel")

### Parameter Modification
When users request changes, the system:
1. Uses LLM to interpret the request
2. Maps natural language to parameter changes
3. Applies modifications with proper nesting
4. Re-executes the current step

## Output Structure

### Directory Organization
```
data/agent_run/{polymer_type}_workflow_{timestamp}/
├── polymer.xyz              # Initial polymer configuration
├── system.xyz               # Polymer + solvent system
├── polymer_plot.png         # Initial polymer visualization
├── system_plot.png          # System with solvent visualization
├── simulation/
│   ├── log.lammps          # LAMMPS log file
│   ├── dump_*.lammpstrj    # Trajectory files
│   └── final_config.xyz     # Final configuration
├── analysis/
│   ├── results.json        # Analysis metrics
│   └── analysis_plot.png   # Analysis visualizations
├── final_plots/
│   ├── final_no_solvent.png
│   └── final_with_solvent.png
└── report.md               # Comprehensive report
```

## Advantages of Multi-Agent Design

### Autonomous Mode Benefits
- **Simplicity**: Direct natural language to results pipeline
- **Efficiency**: Complete pipeline execution without user intervention
- **Consistency**: Single agent handles parameter parsing and execution
- **Scalability**: Easy to batch process multiple simulations

### Interactive Mode Benefits
- **User Control**: Step-by-step validation and parameter adjustment
- **Flexibility**: Each agent handles natural language context independently
- **Learning**: Users can understand each phase of the simulation
- **Quality Assurance**: Human oversight prevents erroneous simulations

## Future Extensions

### Planned Enhancements
1. **Advanced Parameter Optimization**: Agent that suggests optimal parameters based on polymer type
2. **Comparative Analysis**: Multi-simulation comparison agent
3. **Real-time Monitoring**: Live simulation progress tracking
4. **Collaborative Mode**: Multiple users contributing to single workflow

### Integration Points
- **Database Storage**: Persistent storage of simulation results
- **Web Interface**: Browser-based interactive mode
- **API Endpoints**: RESTful API for external integrations
- **Cloud Execution**: Distributed simulation on HPC clusters



