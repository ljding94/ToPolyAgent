# src/agents/interactive_crew.py
import os
from datetime import datetime

from crewai import Crew, Task

from .agents import ConfigAgent, SimulationAgent, WriterAgent
from ..utils.llm_utils import setup_llm
from ..wrappers.custom_tools import HumanFeedbackTool


def run_interactive_workflow(prompt: str, api_key: str = None, model: str = None, auto_mode: bool = False, auto_responses: list = None) -> str:
    """
    Run an interactive polymer simulation workflow driven by agent-native feedback.

    Args:
        prompt (str): Natural language description of the simulation.
        api_key (str): OpenRouter API key (optional, uses env var if not provided).
        model (str): LLM model to use (optional, defaults to gpt-4o-mini).
        auto_mode (bool): If True, use predefined responses for demo purposes.
        auto_responses (list): List of predefined responses for auto mode (optional, uses defaults if not provided).

    Returns:
        str: Path to the generated report.
    """
    llm = None
    if api_key or os.getenv('OPENROUTER_API_KEY'):
        llm = setup_llm(api_key, model)

    # Create workflow directory
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    workflow_name = f"interactive_workflow_{timestamp}"
    workflow_dir = os.path.join("data", "agent_run", workflow_name)
    os.makedirs(workflow_dir, exist_ok=True)
    print(f"📂 Workflow directory created at: {workflow_dir}")

    # Set up auto-responses if in auto mode
    if auto_mode:
        print("\n🤖 Running in AUTO-FEEDBACK MODE")
        print("   Predefined responses will be used for demonstration.\n")
        if auto_responses is None:
            auto_responses = [
                "make the backbone 2x long",  # First feedback: modify config (Trial 1->2)
                "make the solvent density 0.7",  # Second feedback: modify config (Trial 2->3)
                "proceed",  # Third feedback: accept modified config
                "let's use pp interaction 0.3, ss 0.4 and sp 1.5, use langevin thermostat, run for 10000 steps",  # Fourth feedback: modify simulation
                "looks good!",  # Fifth feedback: accept simulation results
            ]
        HumanFeedbackTool.set_auto_responses(auto_responses)
    else:
        HumanFeedbackTool.clear_auto_responses()

    # Create interaction log file
    interaction_log_path = os.path.join(workflow_dir, "interaction_log.txt")
    with open(interaction_log_path, 'w') as f:
        f.write(f"ToPolyAgent Interactive Workflow - {timestamp}\n")
        f.write("=" * 60 + "\n")
        f.write(f"Initial Prompt: {prompt}\n")
        f.write("=" * 60 + "\n\n")

    # Set the log path for HumanFeedbackTool
    HumanFeedbackTool.set_log_path(interaction_log_path)

    print(f"📝 Interaction log: {interaction_log_path}\n")

    # Initialize agents
    config_agent = ConfigAgent(llm, workflow_dir)
    simulation_agent = SimulationAgent(llm)
    writer_agent = WriterAgent(llm)

    # Define the tasks with human feedback integrated
    config_task = Task(
        description=f"""
        You are generating polymer configurations iteratively with human feedback.

        1. Parse the user's initial prompt: '{prompt}'.
        2. OUTPUT DIRECTORY: All files MUST be saved to: '{workflow_dir}'
           You MUST pass 'output_dir="{workflow_dir}"' to ALL tools.

        3. TRIAL TRACKING - Keep a trial counter starting from 1:
           For EACH trial (iteration), you must REGENERATE ALL FILES:

           a) Generate polymer configuration using the appropriate tool.
              The tool will create the polymer, pack it with solvent, and generate two plots:
              - polymer_<type>_nosolvent.png (polymer only)
              - system_<type>.png (polymer + solvent)
              It returns the path to the final system data file (e.g., system_linear.data).

        4. Present results to user with HumanFeedbackTool:
           - Show current trial number.
           - Show all generated file paths (system data and both plots).
           - Show parameters used.
           - Ask: "Are you satisfied with trial N, or would you like modifications?"

        5. Handle user feedback:
           - If user says "yes", "proceed", "looks good", etc. → Accept and finish.
           - If user requests changes → Increment trial counter, update parameters, go back to step 3.
           - Keep iterating until user approves.

        6. Final output must be a JSON string with:
           {{"trial_number": N, "polymer_file": "...", "system_file": "...",
             "polymer_plot": "...", "system_plot": "...", "params": {{...}}}}
        """,
        agent=config_agent,
        expected_output="A JSON string with final trial number, all file paths (including both polymer-only and system plots), and parameters."
    )

    simulation_task = Task(
        description=f"""
        1. Take the approved configuration results from the previous step (config_task context).
        2. Extract the 'system_file' path from the previous task's output - this is the datafile_path for simulation.

        3. OUTPUT DIRECTORY: All files MUST be saved to: '{workflow_dir}'
           You MUST pass 'output_dir="{workflow_dir}"' to the simulation tool.

        4. FIRST: Ask user for simulation parameters using HumanFeedbackTool:
           - Ask: "What simulation parameters would you like to use? Please specify interaction type (good/ideal/poor solvent), thermostat (langevin/nose-hoover), and number of run steps."
           - Parse the response to extract:
             * interaction: 'good', 'ideal', or 'poor' → interaction_params accordingly
             * OR custom parameters: look for specific pp/ss/sp values mentioned
             * thermostat: 'langevin' or 'nose-hoover'
             * run_steps: integer number

           For custom interaction parameters, extract numerical values:
           - Look for "polymer-polymer interaction X" or "pp X" → set pp = X
           - Look for "solvent-solvent interaction X" or "ss X" → set ss = X
           - Look for "polymer-solvent interaction X" or "ps X" or "sp X" → set sp = X
           - If no custom values provided, use defaults based on interaction type

        5. SIMULATION EXECUTION:
           Use the CompleteSimulationPipeline tool to run the entire simulation workflow in one step:
           - datafile_path: use the system_file from config results
           - thermostat: from user input (default: 'langevin')
           - interaction_params: based on user input or defaults
           - run_steps: from user input (default: 20000)
           - output_dir: '{workflow_dir}'

           The tool will automatically:
           * Run LAMMPS simulation
           * Plot final configurations (with and without solvent)
           * Run trajectory analysis
           * Plot analysis results

        6. Present results to user with HumanFeedbackTool:
           - Show simulation parameters used (thermostat, run_steps, interaction_params)
           - Show all generated file paths from the simulation pipeline
           - Ask: "Are you satisfied with the simulation results, or would you like modifications?"

        7. Handle user feedback:
           - If user says "yes", "proceed", "looks good", etc. → Accept and finish
           - If user requests changes (e.g., "change to theta solvent", "run longer") → Update parameters and re-run
           - Keep iterating until user approves

        8. Final output must be a JSON string with:
           {{"simulation_results": "...", "final_config_plots": "...",
             "analysis_results_path": "...", "analysis_plot_path": "..."}}
        """,
        agent=simulation_agent,
        context=[config_task],
        expected_output="A JSON string with simulation results and all generated file paths."
    )

    writer_task = Task(
        description=f"""
        Generate a comprehensive report based on the simulation and analysis results.

        STEP 1: Read the interaction log for workflow context
           - Use ReadFileTool(file_path='{interaction_log_path}')
           - This contains the full history of user interactions, polymer type, and modifications

        STEP 2: Extract analysis results path
           - Get 'analysis_json_path' from the previous task output

        STEP 3: Generate comprehensive report
           - Use ReportTool(analysis_path=<path>, interaction_context=<summary from log>)
           - Include in interaction_context:
             * Initial user request
             * Polymer type generated
             * Any modifications made through feedback
             * Final parameters used

        The ReportTool will automatically:
           - Extract polymer type from analysis data
           - Include system metrics (Rg, Ree, etc.)
           - Include initial and final configuration images
           - Include conformation analysis plot
           - Reference the interaction log file

        Return the complete Markdown report content.
        """,
        agent=writer_agent,
        context=[simulation_task],
        expected_output="Complete Markdown report with workflow context and analysis results."
    )

    # Create memory directory if it doesn't exist
    memory_dir = os.path.join("data", "memory")
    os.makedirs(memory_dir, exist_ok=True)

    # Create and run the main crew
    interactive_crew = Crew(
        agents=[config_agent, simulation_agent, writer_agent],
        tasks=[config_task, simulation_task, writer_task],
        verbose=auto_mode,  # Only verbose in auto mode for debugging
        memory_config={
            "provider": "chromadb",
            "config": {"persist_directory": memory_dir}
        }
    )

    print("\n🚀 Kicking off interactive workflow...")
    print("=" * 40)
    final_result = interactive_crew.kickoff()

    # Add a workflow summary to the interaction log for the report agent
    with open(interaction_log_path, 'a') as f:
        f.write("\n" + "=" * 60 + "\n")
        f.write("WORKFLOW SUMMARY\n")
        f.write("=" * 60 + "\n")
        f.write(f"Initial Prompt: {prompt}\n")
        f.write(f"Workflow Directory: {workflow_dir}\n")
        f.write(f"Timestamp: {timestamp}\n")
        f.write("\nThis log contains the complete interaction history.\n")
        f.write("The report agent should reference this file for workflow context.\n")
        f.write("=" * 60 + "\n")

    # The result from the crew should be the report string from the last task
    report_content = str(final_result)

    # Save the final report
    report_path = os.path.join(workflow_dir, "report.md")
    with open(report_path, 'w') as f:
        f.write(report_content)

    # Create a workflow summary
    summary_path = os.path.join(workflow_dir, "workflow_summary.txt")
    with open(summary_path, 'w') as f:
        f.write("ToPolyAgent Interactive Workflow Summary\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"Timestamp: {timestamp}\n")
        f.write(f"Initial Prompt: {prompt}\n")
        f.write(f"Mode: {'Auto-feedback (Demo)' if auto_mode else 'Interactive'}\n\n")
        f.write(f"Output Directory: {workflow_dir}\n\n")
        f.write("Generated Files:\n")
        f.write("-" * 60 + "\n")

        # List all files in the workflow directory
        try:
            for item in sorted(os.listdir(workflow_dir)):
                item_path = os.path.join(workflow_dir, item)
                if os.path.isfile(item_path):
                    size = os.path.getsize(item_path)
                    f.write(f"  - {item} ({size} bytes)\n")
        except Exception as e:
            f.write(f"  Error listing files: {e}\n")

        f.write("\n" + "=" * 60 + "\n")
        f.write("\nFor detailed interaction history, see: interaction_log.txt\n")
        f.write("For full report, see: report.md\n")

    print("\n🎉 Interactive workflow completed!")
    print(f"📄 Final report saved to: {report_path}")
    print(f"📋 Workflow summary saved to: {summary_path}")
    print(f"📝 Interaction log saved to: {interaction_log_path}")

    return report_path
