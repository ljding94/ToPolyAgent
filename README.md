# ToPolyAgent: AI Agents for Coarse-Grained Topological Polymer Simulations

🤖 **Multi-Agent AI Framework for Polymer Simulations**

[![arXiv](https://img.shields.io/badge/arXiv-2510.12091-b31b1b.svg)](https://arxiv.org/abs/2510.12091)

ToPolyAgent is a cutting-edge multi-agent AI system that enables natural language-driven coarse-grained molecular dynamics simulations of topological polymers. By integrating Large Language Models (LLMs) with specialized computational tools, it democratizes complex polymer simulations for researchers across disciplines.

## 📖 Understanding the Repository

To understand this repository, and let AI help you navigate and setup the project,
use [gitingest](https://gitingest.com/github.com/ljding94/ToPolyAgent) to create a LLM-friendly version of this repo.
Then your coding agent should be able to help you with any questions about the codebase!

## ✨ Key Features

- **Natural Language Interface**: Describe polymer simulations in plain English
- **Multi-Agent Architecture**: Specialized agents for configuration, simulation, analysis, and reporting
- **Dual Operation Modes**: Interactive (with user feedback) and autonomous workflows
- **Comprehensive Polymer Support**: Linear, ring, brush, star, and dendrimer architectures
- **Full Analysis Pipeline**: Conformation analysis, scattering functions, and visualization
- **Extensible Design**: Modular tools and wrappers for easy customization

## 📋 Table of Contents

- [Understanding the Repository](#understanding-the-repository)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Usage](#usage)
- [Architecture](#architecture)
- [Examples](#examples)
- [Testing](#testing)
- [Contributing](#contributing)
- [Citation](#citation)
- [License](#license)

## 🚀 Installation

### Prerequisites

- Python 3.8+
- LAMMPS (molecular dynamics simulator)
- OpenRouter API key (for LLM integration)

### Install LAMMPS

ToPolyAgent uses the `lmp_serial` command to run LAMMPS simulations. Ensure `lmp_serial` is installed and available in your PATH.

```bash
# Using conda (recommended)
conda install -c conda-forge lammps

# Alternative: Using Homebrew (macOS)
brew install lammps

# Verify installation
which lmp_serial

# If lmp_serial is not found, but you have another LAMMPS executable (e.g., lmp_mpi), create an alias:
# alias lmp_serial='lmp_mpi'
# Or add to your shell profile (e.g., ~/.bashrc or ~/.zshrc):
# echo "alias lmp_serial='lmp_mpi'" >> ~/.zshrc
# source ~/.zshrc
```

### Install ToPolyAgent

```bash
# Clone the repository
git clone https://github.com/ljding94/ToPolyAgent.git
cd ToPolyAgent

# Install Python dependencies
pip install -e .
```

### LLM Integration (Optional but Recommended)

1. Get an API key from [OpenRouter](https://openrouter.ai/)
2. Set the environment variable:
   ```bash
   export OPENROUTER_API_KEY="your-api-key-here"
   ```

## 🎯 Quick Start

### Autonomous Mode (Recommended)

Run complete simulations with natural language prompts:

```python
from src.agents.main_crew import run_autonomous

# Simple linear polymer simulation
report_path = run_autonomous("Simulate a linear polymer with 50 beads in good solvent")

# Complex brush polymer with custom parameters
report_path = run_autonomous(
    "Create a brush polymer with backbone of 100 beads, grafting density 0.3, "
    "side chains of 20 beads each, in theta solvent conditions"
)
```

### Interactive Mode

For step-by-step control with user feedback:

```bash
python ToPolyAgent_cli.py --mode interactive
```

### Command Line Interface

```bash
# Autonomous mode
python ToPolyAgent_cli.py --mode autonomous

# Interactive mode
python ToPolyAgent_cli.py --mode interactive
```

## 📖 Usage

### Supported Polymer Architectures

ToPolyAgent supports five major topological polymer types:

- **Linear**: Simple chain polymers
- **Ring**: Cyclic polymers
- **Star**: Polymers with arms radiating from a core
- **Brush**: Polymers with side chains grafted to a backbone
- **Dendrimer**: Branched, tree-like polymers

### Example Prompts

```python
# Linear polymers
"Simulate a linear polymer with 100 monomers in good solvent"

# Star polymers
"Create a star polymer with 6 arms, each 15 beads long, in theta conditions"

# Brush polymers
"Run a brush polymer simulation with backbone of 80 beads, grafting density 0.25, side chains of 12 beads"

# Ring polymers
"Simulate a ring polymer with 60 beads in poor solvent conditions"

# Dendrimers
"Create a dendrimer with 4 generations using good solvent parameters"
```

### Advanced Parameters

Control simulation conditions, thermostats, and analysis:

```python
# Custom thermostat and solvent conditions
report_path = run_autonomous(
    "Simulate a star polymer with 4 arms of 20 beads each using Nose-Hoover thermostat "
    "at temperature 1.0, with solvent quality epsilon=0.8"
)

# Extended simulation time
report_path = run_autonomous(
    "Run a long simulation of linear polymer with 50 beads for 1 million steps "
    "to study diffusion properties"
)
```

## 🏗️ Architecture

ToPolyAgent employs a sophisticated multi-agent system built on CrewAI:

### Core Agents

- **Parser Agent**: Converts natural language to simulation parameters using LLM (primary) or regex (fallback)
- **Structure Agent**: Generates polymer configurations and packs solvent molecules
- **Simulation Agent**: Executes LAMMPS molecular dynamics simulations
- **Analysis Agent**: Performs conformational analysis (Rg, P(q), diffusion, persistence length)
- **Reporting Agent**: Generates comprehensive Markdown reports with visualizations

### Tool Ecosystem

- **Configuration Tools**: Polymer generation, solvent packing, visualization
- **Simulation Tools**: LAMMPS execution with customizable parameters
- **Analysis Tools**: Trajectory analysis, scattering functions, statistical mechanics
- **Wrapper Tools**: CrewAI-compatible interfaces for all functionality

### Data Flow

```
Natural Language Prompt → Parser Agent → Structure Agent → Simulation Agent → Analysis Agent → Report Agent
                                      ↓              ↓              ↓              ↓              ↓
                                 Config Files → LAMMPS Input → Trajectory → Analysis Results → Markdown Report
```

## 📚 Examples

Explore the `examples/` directory for comprehensive demonstrations:

- `autonomous_workflow_example.py`: End-to-end autonomous simulations
- `interactive_workflow_example.py`: Interactive mode with user feedback
- `llm_integration_example.py`: Advanced LLM parsing capabilities
- `generate_report.py`: Report generation utilities

Run any example:

```bash
python examples/autonomous_workflow_example.py
```

## 🧪 Testing

ToPolyAgent includes comprehensive test suites:

```bash
# Run all tests
pytest

# Run specific test categories
pytest tests/test_tools/      # Tool functionality tests
pytest tests/test_agents/     # Agent integration tests
pytest tests/test_utils/      # Utility function tests
pytest tests/test_wrappers/   # Wrapper tool tests
```

## 🤝 Contributing

We welcome contributions! Please see our [Contributing Guide](CONTRIBUTING.md) for details.

### Development Setup

```bash
# Fork and clone
git clone https://github.com/your-username/ToPolyAgent.git
cd ToPolyAgent

# Create virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install in development mode
pip install -e .

# Run tests
pytest
```

## 📄 Citation

If you use ToPolyAgent in your research, please cite our paper:

```bibtex
@article{ding2025topolyagent,
      title={{ToPolyAgent}: {AI Agents} for Coarse-Grained Topological Polymer Simulations},
      author={Ding, Lijie and Carrillo, Jan-Michael and Do, Changwoo},
      journal={arXiv preprint arXiv:2510.12091},
      year={2025},
}
```

## 📄 License

MIT License.

---

**ToPolyAgent**: Bridging natural language with computational polymer science through AI agents.
