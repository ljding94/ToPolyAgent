# src/utils/llm_utils.py
import os
from crewai import LLM


def setup_llm(api_key: str = None, model: str = None) -> LLM:
    """Setup LLM for CrewAI agents using OpenRouter with configurable API key and model"""
    # Use provided API key or fall back to environment variable
    final_api_key = api_key or os.getenv('OPENROUTER_API_KEY')
    if not final_api_key:
        raise ValueError("No API key provided. Either pass api_key parameter or set OPENROUTER_API_KEY environment variable")

    # Use provided model or fall back to default
    #final_model = model or "openai/gpt-4o"
    final_model = model or "openai/gpt-4o-mini"
    #final_model = model or "x-ai/grok-4"
    # newer model like gpt-5 seems not working, might due to crewAI compatibility

    print(f"Using LLM model: {final_model}")
    # For OpenRouter, we need to specify the provider as "openrouter" for LiteLLM
    llm = LLM(
        model=f"openrouter/{final_model}",
        api_key=final_api_key,
        base_url="https://openrouter.ai/api/v1",
        #temperature=0.1,
        #max_tokens=2000
    )
    return llm


def get_default_llm() -> LLM:
    """Get default LLM instance for the application"""
    try:
        return setup_llm()
    except ValueError:
        # If no API key, return None - agents will use their default
        return None
