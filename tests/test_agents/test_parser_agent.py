# tests/test_agents/test_parser_agent.py

def test_workflow_agent_parsing():
    """Test that WorkflowAgent can be instantiated."""
    from src.agents.agents import WorkflowAgent

    # This tests that the agent can be created (LLM parsing happens during task execution)
    agent = WorkflowAgent()
    assert agent is not None
    assert hasattr(agent, 'llm')
    assert len(agent.tools) == 1  # FullWorkflowTool


def test_config_agent_parsing():
    """Test that ConfigAgent can be instantiated."""
    from src.agents.agents import ConfigAgent

    # This tests that the agent can be created (LLM parsing happens in generate_config method)
    agent = ConfigAgent()
    assert agent is not None
    assert hasattr(agent, 'llm')
    assert len(agent.tools) == 6  # 5 polymer tools + PackSolventTool
