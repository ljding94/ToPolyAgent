# src/agents/autonomous_crew.py
from crewai import Crew, Task
from src.agents.agents import WorkflowAgent, WriterAgent
from src.utils.llm_utils import setup_llm
import os
from datetime import datetime


def run_autonomous_workflow(prompt: str, api_key: str = None, model: str = None) -> str:
    """
    Run the autonomous workflow: parse prompt, execute full simulation, generate report.

    Args:
        prompt (str): Natural language description of the simulation.
        api_key (str): OpenRouter API key (optional, uses env var if not provided).
        model (str): LLM model to use (optional, defaults to gpt-4o-mini).

    Returns:
        str: Path to the generated report.
    """
    # Setup LLM if API key provided
    llm = None
    if api_key or os.getenv('OPENROUTER_API_KEY'):
        llm = setup_llm(api_key, model)

    # Initialize agents
    workflow_agent = WorkflowAgent(llm)
    writer_agent = WriterAgent(llm)

    # Define project root and create a single, unique output directory for this run
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
    agent_run_dir = os.path.join(project_root, "data", "agent_run")
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_output_dir = os.path.join(agent_run_dir, f"autonomous_workflow_{timestamp}")
    os.makedirs(run_output_dir, exist_ok=True)

    # Task 1: Run complete workflow
    workflow_task = Task(
        description=f"""1. Parse this natural language prompt to extract polymer simulation parameters: '{prompt}'.
        2. Call the 'FullWorkflowTool' with the following parameters:
            - output_dir: '{run_output_dir}'
            - polymer_type: the type of polymer (e.g., 'linear', 'ring', 'brush', 'star', 'dendrimer')
            - polymer_params_json: JSON string with polymer-specific parameters
            - box_size, solvent_density, run_steps, thermostat, interaction_params_json: optional parameters with defaults
        The final output MUST be the dictionary returned by the tool.""",
        agent=workflow_agent,
        expected_output="A dictionary containing the full results of the workflow, including the 'output_dir' key which must match the provided path."
    )

    # Task 2: Generate report
    report_task = Task(
        description="Use the ReportTool to generate a comprehensive Markdown report from the analysis results. Call ReportTool with the analysis_results.json path from the workflow output. The report should include all configurations and analysis plots with relative paths.",
        agent=writer_agent,
        expected_output="A string containing the complete Markdown report.",
        context=[workflow_task]
    )

    # Create and run crew
    crew = Crew(
        agents=[workflow_agent, writer_agent],
        tasks=[workflow_task, report_task],
        verbose=True
    )

    # Execute workflow
    result = crew.kickoff()

    # The result from the crew should be the report string from the last task
    report_content = str(result)

    # Save the final report in the directory created at the start
    report_path = os.path.join(run_output_dir, "report.md")
    with open(report_path, 'w') as f:
        f.write(report_content)

    return report_path
