#!/usr/bin/env python3
"""
Test script for parameter modification in interactive mode
"""

import sys
import os

# Add src to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'src'))

from src.utils.llm_utils import setup_llm


def test_parameter_modification():
    """Test the parameter modification logic."""
    print("Testing parameter modification logic...")

    # Setup LLM
    api_key = os.getenv('OPENROUTER_API_KEY')
    if not api_key:
        print("❌ No OpenRouter API key found")
        return

    llm = setup_llm(api_key, "openai/gpt-4o-mini")
    if not llm:
        print("❌ Failed to setup LLM")
        return

    # Simulate workflow state
    workflow_state = {
        'params': {
            'polymer_type': 'linear',
            'polymer_params': {'N': 20},
            'solvent_density': 0.2,
            'thermostat': 'langevin',
            'interaction_params': {'pp': 0.5, 'ss': 0.5, 'sp': 1.5},
            'run_steps': 10000,
            'analysis_metrics': ['Rg', 'P(q)'],
            'box_size': 20.0
        }
    }

    # Test feedback
    feedback = "make the chain longer"
    step_name = "Structure Generation"

    print(f"Original params: {workflow_state['params']}")

    # Test the handle_parameter_modification function
    from src.agents.interactive_crew import handle_parameter_modification

    modifications = handle_parameter_modification(llm, step_name, feedback, workflow_state)
    print(f"Modifications returned: {modifications}")

    if modifications:
        # Apply the deep update logic
        def deep_update(base_dict, update_dict):
            for key, value in update_dict.items():
                if isinstance(value, dict) and key in base_dict and isinstance(base_dict[key], dict):
                    deep_update(base_dict[key], value)
                else:
                    base_dict[key] = value

        deep_update(workflow_state["params"], modifications)
        print(f"Updated params: {workflow_state['params']}")
    else:
        print("No modifications applied")


if __name__ == "__main__":
    test_parameter_modification()
