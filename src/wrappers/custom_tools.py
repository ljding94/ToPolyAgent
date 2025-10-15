# src/wrappers/custom_tools.py
from crewai.tools.base_tool import BaseTool
from typing import Optional
import json
import os
from src.utils.logging import setup_logger

# Rich UI imports
from rich.console import Console
from rich.panel import Panel
from rich.prompt import Prompt

logger = setup_logger()


class ReportTool(BaseTool):
    name: str = "GenerateReport"
    description: str = """Compiles analysis JSON and plots into Markdown report.
    Inputs: analysis_path (str: path to analysis_results.json), interaction_context (str: optional summary from interaction log).
    Outputs: Markdown content (str)."""

    def _run(self, analysis_path: str, interaction_context: str = "") -> str:
        try:
            # Get the base directory from analysis_path
            base_dir = os.path.dirname(analysis_path)

            # Load analysis results
            with open(analysis_path, 'r') as f:
                analysis_data = json.load(f)

            # Load simulation parameters if available
            sim_params = None
            sim_params_path = os.path.join(base_dir, "simulation_parameters.json")
            if os.path.exists(sim_params_path):
                try:
                    with open(sim_params_path, 'r') as f:
                        sim_params = json.load(f)
                except Exception as e:
                    logger.warning(f"Failed to load simulation parameters: {e}")

            # Extract polymer type from system_info
            polymer_type = "Unknown"
            if 'system_info' in analysis_data and 'polymer_type' in analysis_data['system_info']:
                polymer_type = analysis_data['system_info']['polymer_type'].capitalize()

            # Check for interaction log for additional context
            interaction_log_path = os.path.join(base_dir, "interaction_log.txt")
            workflow_context = ""
            if os.path.exists(interaction_log_path):
                try:
                    with open(interaction_log_path, 'r') as f:
                        log_content = f.read()
                        # Extract initial prompt
                        if "Initial Prompt:" in log_content:
                            for line in log_content.split('\n'):
                                if line.startswith("Initial Prompt:"):
                                    workflow_context = line.replace("Initial Prompt:", "").strip()
                                    break
                except Exception:
                    pass

            # Generate Markdown report
            report = "# Polymer Simulation Analysis Report\n\n"

            # Add workflow context if provided by agent
            if interaction_context:
                report += "## Workflow Context\n"
                report += f"{interaction_context}\n\n"
            # Otherwise try to extract from interaction log file
            elif workflow_context:
                report += "## Workflow Context\n"
                report += f"**Initial Request**: {workflow_context}\n\n"
                report += "*For detailed interaction history, see `interaction_log.txt`*\n\n"

            # Summary
            report += "## Summary\n"
            report += f"- **Polymer Type**: {polymer_type}\n"
            report += f"- **Number of Polymer Atoms**: {analysis_data.get('system_info', {}).get('n_polymer_atoms', 'N/A')}\n"
            report += f"- **Simulation Success**: {analysis_data['metadata'].get('success', False)}\n\n"

            # Simulation Parameters Section
            if sim_params:
                report += "## Simulation Parameters\n\n"

                # Polymer Configuration
                report += "### Polymer Configuration\n"
                report += f"- **Type**: {sim_params.get('polymer_type', 'N/A')}\n"

                # Add polymer-specific parameters
                poly_params = sim_params.get('polymer_params', {})
                if 'chain_length' in poly_params:
                    report += f"- **Chain Length**: {poly_params['chain_length']} beads\n"
                elif 'arm_length' in poly_params and 'num_arms' in poly_params:
                    report += f"- **Number of Arms**: {poly_params['num_arms']}\n"
                    report += f"- **Arm Length**: {poly_params['arm_length']} beads\n"
                elif 'backbone_length' in poly_params:
                    report += f"- **Backbone Length**: {poly_params['backbone_length']} beads\n"
                    report += f"- **Side Chain Length**: {poly_params.get('side_chain_length', 'N/A')} beads\n"
                    report += f"- **Grafting Density**: {poly_params.get('grafting_density', 'N/A')}\n"
                elif 'generation' in poly_params:
                    report += f"- **Generation**: {poly_params['generation']}\n"
                    report += f"- **Branching Factor**: {poly_params.get('branching_factor', 'N/A')}\n"

                report += f"- **Box Size**: {sim_params.get('box_size', 'N/A')} (LJ units)\n"
                report += f"- **Solvent Density**: {sim_params.get('solvent_density', 'N/A')}\n\n"

                # Simulation Settings
                sim_settings = sim_params.get('simulation', {})
                report += "### Simulation Settings\n"
                report += f"- **Total Steps**: {sim_settings.get('run_steps', 'N/A')}\n"
                report += f"- **Thermostat**: {sim_settings.get('thermostat', 'N/A')}\n"
                report += f"- **Temperature**: {sim_settings.get('temperature', 'N/A')} (LJ units)\n"
                report += f"- **Timestep**: {sim_settings.get('timestep', 'N/A')} (LJ units)\n\n"

                # Interaction Parameters
                interactions = sim_params.get('interactions', {})
                report += "### Interaction Parameters\n"
                report += f"- **Polymer-Polymer (ε_pp)**: {interactions.get('polymer_polymer', 'N/A')}\n"
                report += f"- **Solvent-Solvent (ε_ss)**: {interactions.get('solvent_solvent', 'N/A')}\n"
                report += f"- **Polymer-Solvent (ε_ps)**: {interactions.get('polymer_solvent', 'N/A')}\n"

                if 'description' in interactions:
                    report += f"\n*{interactions['description']}*\n"
                report += "\n"

            # Metrics
            report += "## Analysis Metrics\n"

            # Extract metrics from analysis_results section
            if 'analysis_results' in analysis_data:
                results = analysis_data['analysis_results']

                # Radius of Gyration
                if 'mean_radius_of_gyration' in results:
                    rg_mean = results['mean_radius_of_gyration']
                    report += f"- **Radius of Gyration (Rg)**: Mean = {rg_mean:.3f}\n"

                # End-to-End Distance
                if 'mean_end_to_end_distance' in results:
                    ree_mean = results['mean_end_to_end_distance']
                    report += f"- **End-to-End Distance (Ree)**: Mean = {ree_mean:.3f}\n"

                # Persistence Length
                if 'mean_persistence_length' in results:
                    lp_mean = results['mean_persistence_length']
                    report += f"- **Persistence Length (Lp)**: Mean = {lp_mean:.3f}\n"

                # Diffusion Coefficient
                if 'diffusion_coefficient' in results:
                    diff_coeff = results['diffusion_coefficient']
                    if str(diff_coeff).lower() != 'nan':
                        report += f"- **Diffusion Coefficient (D)**: {diff_coeff:.2e}\n"
                    else:
                        report += "- **Diffusion Coefficient (D)**: Not available\n"

                # Structure factor
                if 'scattering_pq' in results or 'mean_scattering_pq' in results:
                    report += "- **Structure Factor P(q)**: Computed and plotted\n"

                # Radial distribution function
                if 'radial_distribution_gr' in results or 'mean_radial_distribution_gr' in results:
                    report += "- **Radial Distribution Function g(r)**: Computed and plotted\n"

            report += "\n"

            # Plots
            report += "## Configurations and Analysis\n\n"

            # Initial configuration
            polymer = analysis_data['system_info']['polymer_type']
            initial_config_nosolvent = f"polymer_{polymer}_nosolvent.png"
            initial_config_with_solvent = f"system_{polymer}.png"

            report += "### Initial Polymer Configuration\n"
            if os.path.exists(os.path.join(base_dir, initial_config_nosolvent)):
                report += f'<img src="{initial_config_nosolvent}" width="30%" style="display: inline-block; margin: 10px;" alt="Initial Polymer Configuration">\n'
            if os.path.exists(os.path.join(base_dir, initial_config_with_solvent)):
                report += f'<img src="{initial_config_with_solvent}" width="30%" style="display: inline-block; margin: 10px;" alt="Initial System Configuration">\n'
            if not os.path.exists(os.path.join(base_dir, initial_config_nosolvent)) and not os.path.exists(os.path.join(base_dir, initial_config_with_solvent)):
                report += "Initial configuration images not found.\n"
            report += "\n"

            # Final configuration
            final_config_nosolvent = "final_state_nosolvent.png"
            final_config_with_solvent = "final_state.png"

            report += "### Final Polymer Configuration\n"
            if os.path.exists(os.path.join(base_dir, final_config_nosolvent)):
                report += f'<img src="{final_config_nosolvent}" width="30%" style="display: inline-block; margin: 10px;" alt="Final State Without Solvent">\n'
            if os.path.exists(os.path.join(base_dir, final_config_with_solvent)):
                report += f'<img src="{final_config_with_solvent}" width="30%" style="display: inline-block; margin: 10px;" alt="Final State With Solvent">\n'
            if not os.path.exists(os.path.join(base_dir, final_config_nosolvent)) and not os.path.exists(os.path.join(base_dir, final_config_with_solvent)):
                report += "Final configuration images not found.\n"
            report += "\n"

            plot_path = "conformation_analysis.png"
            # Conformation analysis plot
            if os.path.exists(os.path.join(base_dir, plot_path)):
                report += "### Conformation Analysis\n"
                report += f'<img src="{plot_path}" width="50%" style="display: block; margin: 10px 0;" alt="Conformation Analysis">\n\n'
            else:
                report += "### Conformation Analysis\n"
                report += "Conformation analysis plot not found.\n\n"

            logger.info(f"Generated report from {analysis_path}")
            return report
        except Exception as e:
            logger.error(f"Failed to generate report: {e}")
            return f"# Error Generating Report\n\nError: {str(e)}"


class HumanFeedbackTool(BaseTool):
    name: str = "HumanFeedbackTool"
    description: str = """
    Asks for human feedback. Use this to present results to the user and get their confirmation to proceed, or to get specific instructions for modifications.
    The tool will pause the execution and wait for the user to provide input from the command line.
    You must provide a clear 'message' to the user explaining what you have done and what you need from them.
    You can optionally provide 'outputs' to show the user relevant data or file paths.
    The user's response will be returned as a string.
    """

    # Class variables for automated feedback and logging
    _auto_responses = []
    _response_index = 0
    _interaction_log_path = None

    @classmethod
    def set_auto_responses(cls, responses: list):
        """Set predefined responses for automated demo mode."""
        cls._auto_responses = responses
        cls._response_index = 0

    @classmethod
    def clear_auto_responses(cls):
        """Clear automated responses and return to interactive mode."""
        cls._auto_responses = []
        cls._response_index = 0

    @classmethod
    def set_log_path(cls, log_path: str):
        """Set the path for logging interactions."""
        cls._interaction_log_path = log_path

    def _log_interaction(self, message: str, outputs: dict, response: str):
        """Log the interaction to file if log path is set."""
        if self._interaction_log_path:
            try:
                import datetime
                with open(self._interaction_log_path, 'a') as f:
                    f.write(f"\n[{datetime.datetime.now().strftime('%H:%M:%S')}] AGENT MESSAGE:\n")
                    f.write(f"{message}\n\n")
                    if outputs:
                        f.write("OUTPUTS:\n")
                        for key, value in outputs.items():
                            f.write(f"  - {key}: {value}\n")
                        f.write("\n")
                    f.write(f"USER RESPONSE: {response}\n")
                    f.write("-" * 60 + "\n")
            except Exception as e:
                logger.warning(f"Failed to log interaction: {e}")

    def _run(self, message: str, outputs: Optional[dict] = None) -> str:
        """Pauses to show a message and optional outputs, then waits for user input."""
        # Terminal mode - use Rich UI for better UX
        console = Console()

        # Create a clear visual separation
        console.print()
        console.print(Panel.fit(
            "[bold blue]🤖 ToPolyAgent - Human Feedback Required[/bold blue]",
            border_style="blue",
            title="🔄 Interactive Workflow"
        ))
        console.print()  # Extra spacing

        # Display the agent message
        console.print(Panel(
            message,
            title="[yellow]📢 Agent Message[/yellow]",
            border_style="yellow"
        ))

        # Display current outputs if provided
        if outputs:
            output_lines = []
            for key, value in outputs.items():
                output_lines.append(f"[cyan]• {key}:[/cyan] {value}")

            console.print(Panel(
                "\n".join(output_lines),
                title="[green]📋 Current Outputs[/green]",
                border_style="green"
            ))

        response = ""

        # Check if we're in auto-response mode
        if self._auto_responses and self._response_index < len(self._auto_responses):
            response = self._auto_responses[self._response_index]
            self.__class__._response_index += 1

            console.print(Panel(
                f"[italic]Simulated user input: '{response}'[/italic]",
                title="[red]🤖 Auto-Response Mode[/red]",
                border_style="red"
            ))
        else:
            # Normal interactive mode - use Rich prompt for better UX
            console.print(Panel(
                "[bold white]✓ Type 'yes' or 'proceed' to continue to next step[/bold white]\n"
                "[bold white]✓ Type 'stop' or 'quit' to end the workflow[/bold white]\n"
                "[bold white]✓ Describe changes to modify current step (e.g., 'make polymer longer')[/bold white]",
                title="[magenta]💬 Your Response Needed[/magenta]",
                border_style="magenta"
            ))

            # Use Rich Prompt for better input handling
            response = Prompt.ask("[bold cyan]Your feedback[/bold cyan]").strip()

            console.print(Panel(
                f"[italic]Received: '{response}'[/italic]",
                border_style="cyan"
            ))

        # Log the interaction
        self._log_interaction(message, outputs or {}, response)

        return response


class ReadFileTool(BaseTool):
    name: str = "ReadFile"
    description: str = """Reads the contents of a text file.
    Inputs: file_path (str).
    Outputs: File contents (str)."""

    def _run(self, file_path: str) -> str:
        """Read and return the contents of a file."""
        try:
            if not os.path.exists(file_path):
                return f"Error: File not found at {file_path}"

            with open(file_path, 'r') as f:
                content = f.read()

            logger.info(f"Successfully read file: {file_path}")
            return content
        except Exception as e:
            error_msg = f"Error reading file {file_path}: {str(e)}"
            logger.error(error_msg)
            return error_msg
