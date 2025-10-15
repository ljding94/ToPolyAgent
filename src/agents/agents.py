# src/agents/agents.py
from crewai import Agent
from src.wrappers.workflow_wrappers import FullWorkflowTool
from src.wrappers.config_wrappers import (
    GenerateLinearPolymerTool,
    GenerateRingPolymerTool,
    GenerateBrushPolymerTool,
    GenerateStarPolymerTool,
    GenerateDendrimerTool,
    PackSolventTool
)
from src.wrappers.sim_wrappers import SimulationTool
from src.wrappers.custom_tools import ReportTool, HumanFeedbackTool, ReadFileTool
from src.utils.llm_utils import get_default_llm


class WorkflowAgent(Agent):
    """Agent for autonomous mode - handles natural language parsing and complete workflow execution."""

    def __init__(self, llm=None):
        llm = llm or get_default_llm()
        super().__init__(
            role="Workflow Orchestrator",
            goal="Parse natural language prompts and execute complete polymer simulation workflows",
            backstory="You are an expert at understanding polymer simulation requirements from natural language and orchestrating the complete simulation pipeline from configuration to analysis.",
            verbose=True,
            tools=[FullWorkflowTool()],
            allow_delegation=False,
            llm=llm
        )


class ConfigAgent(Agent):
    """Agent for interactive mode - handles config generation with natural language input."""

    def __init__(self, llm=None, workflow_dir=None):
        llm = llm or get_default_llm()
        super().__init__(
            role="Configuration Specialist",
            goal=f"Parse natural language requirements to generate polymer configurations, which includes packing with solvent and plotting, saving all outputs to the designated directory: {workflow_dir}",
            backstory="You are an expert at creating polymer configurations. You use tools that generate the polymer structure, pack it with solvent, and plot the results in a single step. You MUST use the provided output directory for all generated files.",
            verbose=True,
            tools=[
                GenerateLinearPolymerTool(),
                GenerateRingPolymerTool(),
                GenerateBrushPolymerTool(),
                GenerateStarPolymerTool(),
                GenerateDendrimerTool(),
                HumanFeedbackTool()
            ],
            allow_delegation=False,
            llm=llm
        )
        # Store workflow_dir as a regular attribute (not a Pydantic field)
        object.__setattr__(self, 'workflow_dir', workflow_dir)


class SimulationAgent(Agent):
    """Agent for interactive mode - handles simulation execution and analysis."""

    def __init__(self, llm=None):
        llm = llm or get_default_llm()
        super().__init__(
            role="Simulation & Analysis Specialist",
            goal="Execute the complete simulation pipeline: run LAMMPS simulation, plot final configurations, analyze trajectory, and plot analysis results in a single streamlined process.",
            backstory="You are an expert at running comprehensive polymer simulations. You execute the entire simulation pipeline from molecular dynamics to final analysis visualization in one efficient workflow.",
            verbose=True,
            tools=[
                SimulationTool(),
                HumanFeedbackTool()
            ],
            allow_delegation=False,
            llm=llm
        )


class WriterAgent(Agent):
    """Agent for generating comprehensive reports from simulation results."""

    def __init__(self, llm=None):
        llm = llm or get_default_llm()
        super().__init__(
            role="Report Generator",
            goal="Compile simulation results into comprehensive, well-formatted reports",
            backstory="You are an expert at synthesizing complex simulation data into clear, actionable reports that highlight key findings and insights.",
            verbose=True,
            tools=[ReportTool(), ReadFileTool()],
            allow_delegation=False,
            llm=llm
        )

    def generate_report(self, analysis_json_path: str) -> str:
        """Generate a Markdown report from analysis results."""
        return self.tools[0]._run(analysis_json_path)
