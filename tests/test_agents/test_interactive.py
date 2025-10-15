#!/usr/bin/env python3
"""
Test script for the interactive mode
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'src'))

from src.agents import run_interactive

# Test the interactive mode with a simple prompt
if __name__ == "__main__":
    try:
        print("Testing interactive mode...")
        # This will run but will wait for user input
        # For testing, we'll use a simple prompt
        report_path = run_interactive(
            prompt="Simulate a linear polymer with 10 beads in good solvent",
            api_key=os.getenv('OPENROUTER_API_KEY'),
            model="openai/gpt-4o-mini"
        )
        print(f"Interactive workflow completed. Report: {report_path}")
    except Exception as e:
        print(f"Test failed: {e}")
