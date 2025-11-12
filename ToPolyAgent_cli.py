#!/usr/bin/env python3
"""
ToPolyAgent CLI
===============

A command-line interface for ToPolyAgent that allows you to run
polymer simulations using natural language commands.

Usage:
    python ToPolyAgent_cli.py --mode [interactive|autonomous]

Modes:
    - interactive: Agent asks for feedback at each step
    - autonomous: Runs complete workflow without feedback

After specifying the mode, you will be prompted to enter your simulation description.
"""

import os
import sys
import time
import argparse
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / "src"))

from src.agents import run_interactive, run_autonomous
from src.utils.llm_utils import setup_llm


class ToPolyAgentCLI:
    def __init__(self):
        self.llm = None
        self.api_key = None
        self.model = "openai/gpt-4o-mini"
        self.setup_complete = False

    def setup_system(self):
        """Setup the ToPolyAgent system with LLM integration."""
        print("🚀 Initializing ToPolyAgent...")
        print("=" * 50)

        # Check for API key
        self.api_key = os.getenv('OPENROUTER_API_KEY')
        if not self.api_key:
            print("❌ No OpenRouter API key found!")
            print("Please set the OPENROUTER_API_KEY environment variable:")
            print("  export OPENROUTER_API_KEY='your-api-key-here'")
            print("Get your key from: https://openrouter.ai/")
            return False

        print("✅ OpenRouter API key found")

        # Setup LLM
        try:
            self.llm = setup_llm(api_key=self.api_key, model=self.model)
            print(f"✅ LLM configured: {self.model}")
        except Exception as e:
            print(f"❌ LLM setup failed: {e}")
            return False

        # Check for required dependencies
        self.check_dependencies()

        self.setup_complete = True
        print("🎉 ToPolyAgent is ready!")
        print()
        return True

    def check_dependencies(self):
        """Check for required simulation dependencies."""
        print("🔍 Checking dependencies...")

        # Check LAMMPS (try multiple common executables and paths)
        lammps_available = self.check_lammps()
        if lammps_available:
            print("✅ LAMMPS found")
        else:
            print("⚠️  LAMMPS not found - simulations will fail")
            print("   Install with: conda install -c conda-forge lammps")
            print("   Or check if lmp/lmp_serial/lammps is in your PATH")
            print("   Or add LAMMPS bin directory to your PATH")

        # Note: Packmol is not needed - solvent packing is done algorithmically

        print()

    def check_command(self, command):
        """Check if a command is available."""
        import subprocess
        try:
            subprocess.run([command, "--version"], capture_output=True, check=True)
            return True
        except (subprocess.CalledProcessError, FileNotFoundError):
            return False

    def check_lammps(self):
        """Check for LAMMPS with multiple detection methods."""
        # First try PATH
        if self.check_command("lmp") or self.check_command("lmp_serial") or self.check_command("lammps"):
            return True

        # Check common installation paths
        common_paths = [
            "/opt/homebrew/bin/lmp_serial",
            "/opt/homebrew/bin/lmp",
            "/usr/local/bin/lmp_serial",
            "/usr/local/bin/lmp",
            "/usr/bin/lmp_serial",
            "/usr/bin/lmp",
            "~/miniconda3/bin/lmp_serial",  # conda
            "~/miniconda3/bin/lmp",
            "~/anaconda3/bin/lmp_serial",   # anaconda
            "~/anaconda3/bin/lmp",
        ]

        import os
        for path in common_paths:
            expanded_path = os.path.expanduser(path)
            if os.path.exists(expanded_path) and os.access(expanded_path, os.X_OK):
                return True

        return False

    def run_interactive_simulation(self, prompt):
        """Run an interactive simulation with human feedback loops."""
        print("\n🎯 Starting Interactive Mode")
        print("=" * 50)
        print(f"Prompt: {prompt}")
        print("\n⚠️  Interactive mode workflow:")
        print("   1. 🏗️  CONFIGURATION: Generate polymer → Agent asks for feedback")
        print("   2. 🔬 SIMULATION: Run LAMMPS → Agent asks for feedback")
        print("   3. 📊 ANALYSIS: Analyze results → Generate report")
        print("")
        print("   At each feedback point:")
        print("   • Type 'yes' or 'proceed' to continue to next step")
        print("   • Describe changes to redo the current step (e.g., 'make polymer longer')")
        print("   • The agent will ask you to confirm - stay at your terminal!")
        print("")
        print("⏳ Starting workflow... (please wait for agent feedback prompts)")
        print()

        try:
            start_time = time.time()

            # Run the interactive workflow
            report_path = run_interactive(
                prompt=prompt,
                api_key=self.api_key,
                model=self.model
            )

            end_time = time.time()
            duration = end_time - start_time

            print("✅ Interactive workflow completed!")
            print(f"⏱️  Total duration: {duration:.2f}s")
            print(f"📄 Report: {report_path}")

        except Exception as e:
            print(f"❌ Interactive workflow failed: {e}")
            print("This might be due to missing dependencies or configuration issues.")

        print()

    def run_autonomous_simulation(self, prompt):
        """Run an autonomous simulation without human feedback."""
        print("\n🤖 Starting Autonomous Mode")
        print("=" * 50)
        print(f"Prompt: {prompt}")
        print("\n💡 About Autonomous Mode:")
        print("   • The AI will parse your prompt and run the complete workflow")
        print("   • No human feedback will be requested during execution")
        print("   • The workflow: Config → Simulation → Analysis → Report")
        print("   • This may take several minutes depending on simulation size")
        print("\n⚡ Running fully automated workflow...")
        print()

        try:
            start_time = time.time()

            # Run the autonomous workflow
            report_path = run_autonomous(
                prompt=prompt,
                api_key=self.api_key,
                model=self.model
            )

            end_time = time.time()
            duration = end_time - start_time

            print("✅ Autonomous workflow completed!")
            print(f"⏱️  Total duration: {duration:.2f}s")
            print(f"📄 Report: {report_path}")

        except Exception as e:
            print(f"❌ Autonomous workflow failed: {e}")
            print("This might be due to missing dependencies or configuration issues.")
            import traceback
            traceback.print_exc()

        print()


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(description="ToPolyAgent CLI for polymer simulations")
    parser.add_argument("--mode", choices=["interactive", "autonomous"], required=True,
                        help="Mode to run: interactive (with feedback) or autonomous (no feedback)")

    args = parser.parse_args()

    cli = ToPolyAgentCLI()

    # Setup system
    if not cli.setup_system():
        print("❌ Setup failed. Please fix the issues above and try again.")
        sys.exit(1)

    # Prompt for simulation description
    if args.mode == "interactive":
        print("\n🎯 Interactive Mode Selected")
        print("You will be prompted for feedback during the simulation.")
    else:
        print("\n🤖 Autonomous Mode Selected")
        print("The simulation will run without feedback.")

    prompt = input("\nEnter your simulation description: ").strip()
    if not prompt:
        print("No prompt entered. Exiting.")
        sys.exit(1)

    # Run the simulation based on mode
    if args.mode == "interactive":
        cli.run_interactive_simulation(prompt)
    elif args.mode == "autonomous":
        cli.run_autonomous_simulation(prompt)


if __name__ == "__main__":
    main()
