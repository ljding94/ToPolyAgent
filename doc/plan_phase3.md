### Detailed Implementation Plan for Developing the Agent System in ToPolyAgent

This plan is designed for a coding agent (e.g., you or an AI-assisted developer) to follow step-by-step. It builds on the existing project structure (e.g., `src/tools/`, `src/wrappers/`, `tests/`) and assumes Phase 1 and 2 are complete. We'll add `src/agents/` for agents and related files. The plan is modular, with estimated effort (low/medium/high) and dependencies. Total estimated time: 4-7 days, as per outline.md—focus on autonomous first, then interactive.

Use Python 3.x, CrewAI (from requirements.txt), and existing wrappers. Test incrementally. Commit often with descriptive messages (e.g., "Add ParserAgent and parse_prompt utility").

#### Step 1: Setup Directory and Dependencies (Effort: Low, ~30 min)
- **Goal**: Prepare the structure for agents.
- **Actions**:
  1. Create `src/agents/` if not exists (already in structure, but confirm).
     - Add `__init__.py`.
     - Subfiles: `parser_agent.py`, `structure_agent.py`, `simulation_agent.py`, `analysis_agent.py`, `reporting_agent.py`, `main_crew.py`.
  2. In `src/utils/`, add `prompt_parser.py` (for shared parsing logic). **[UPDATED: Removed - parsing now handled by agents directly]**
  3. Update `requirements.txt` if needed (e.g., ensure `crewai` and `crewai_tools` are there; add if missing).
  4. In `src/__init__.py`, import agents for easy access.
- **Testing**: Run `python -m src.agents` to check no import errors.
- **Dependencies**: None.

#### Step 2: Implement Agent-Based Parsing (Effort: Medium, ~1-2 hours) **[UPDATED: Parsing now integrated into agents]**
- **Goal**: Agents handle natural language understanding directly using LLM capabilities.
- **Actions**:
  1. Each agent (WorkflowAgent, ConfigAgent) includes LLM-based parsing logic
  2. Parsing happens within agent methods or CrewAI task execution
  3. No separate parsing utility needed - agents understand context directly
- **Testing**:
  - Test agents can instantiate and handle natural language prompts
  - Integration tests verify end-to-end functionality
- **Dependencies**: CrewAI LLM integration
  2. Add error handling: If ambiguous, return partial dict with notes (e.g., {"error": "Missing density"}).
- **Testing**:
  - Test agents can instantiate and handle natural language prompts **[UPDATED: No separate parser tests needed]**
  - Integration tests verify end-to-end functionality
- **Dependencies**: CrewAI LLM integration

#### Step 3: Implement Custom Tools (Effort: Low, ~1 hour)
- **Goal**: Add tools not in wrappers (for reporting and HITL).
- **Actions**:
  1. In `src/wrappers/custom_tools.py` (new file):
     - `ReportTool` (CrewAI Tool subclass):
       - Name: "GenerateReport", Description: "Compiles analysis JSON and plots into Markdown report."
       - `_run(self, analysis_path: str) -> str`: Read JSON, format as Markdown (e.g., tables for metrics, embed plot paths as images).
       - Output: str (Markdown content).
     - `HumanFeedbackTool` (extend crewai_tools.HumanInputTool or custom):
       - Name: "GetHumanFeedback", Description: "Pauses to show output and get user input."
       - `_run(self, message: str, output_to_show: str = None) -> str`: Print message + output (e.g., plot path), then `input("User: ")`.
       - For plots: Print path (user opens manually); optional: Use matplotlib to show if possible.
- **Testing**:
  - Add `tests/test_wrappers/test_custom_tools.py`: Mock inputs, assert outputs.
- **Dependencies**: `json`, `input()` for human tool.

#### Step 4: Define Agents (Effort: Medium, ~2-3 hours)
- **Goal**: Create agent classes in `src/agents/`.
- **Actions** (Each agent is a CrewAI Agent subclass):
  1. `parser_agent.py`: Agent(goal="Extract params from prompt", tools=[], llm=your_choice).
     - Task: Use `parse_prompt` utility.
  2. `structure_agent.py`: Agent(goal="Generate and pack config", tools=[PolymerGeneratorTool, PackSolventTool, PlotConfigTool]).
     - Handles: Parse → Generate polymer → Pack solvent → Plot (with/without solvent).
     - Output: Dict with paths (polymer_file, system_file, plots).
  3. `simulation_agent.py`: Agent(goal="Run LAMMPS sim", tools=[RunLammpsTool]).
     - Input: system_file + params (thermostat, interaction_params, run_steps).
     - Output: Dict with dump_files, final_config, log.
  4. `analysis_agent.py`: Agent(goal="Analyze dumps", tools=[ComprehensiveAnalysisTool, PlotAnalysisTool]).
     - Input: system_file, dump_pattern, metrics (from prompt or all).
     - Dynamically: If metrics specified, subset analysis; else full.
     - Output: analysis_json_path, plot_path.
  5. `reporting_agent.py`: Agent(goal="Compile final report", tools=[ReportTool]).
     - Input: All prior outputs.
     - Output: Markdown str.
- **Common**: Each agent has verbose=True for logging; use `src/utils/logging.py`.
- **Testing**: Add `tests/test_agents/test_[agent].py`: Mock tools, assert outputs for sample inputs.

#### Step 5: Implement Autonomous Mode (Effort: Medium, ~1-2 hours)
- **Goal**: End-to-end without user input.
- **Actions**:
  1. In `main_crew.py`, add `run_autonomous(prompt: str) -> str` (as sketched earlier).
     - Create Crew with agents.
     - Sequential tasks: Parse → Structure → Sim → Analysis → Report.
     - Use memory to pass paths (e.g., context from previous task).
     - Handle errors: If parse fails, default or raise.
  2. Output: Save report to file (e.g., "report.md") in workflow dir.
- **Testing**:
  - Add `tests/test_agents/test_autonomous.py`: Run with sample prompt, assert report contains expected sections.
  - Manual: Run with example prompt, check files generated.

#### Step 6: Implement Interactive Mode (Effort: High, ~2-3 hours)
- **Goal**: Step-by-step with loops for feedback.
- **Actions**:
  1. In `main_crew.py`, add `run_interactive(initial_prompt: str = None) -> str` (as sketched earlier).
     - Use while loops for each major step (structure, solvent, sim+analysis, report).
     - In loops: Run agent/task, show output (e.g., print plot path), use HumanFeedbackTool.
     - If "adjust", reprompt for changes → Update params → Retry step.
     - For sim: Propose params first (e.g., "Using langevin, sp=0.7, 50000 steps—confirm?").
     - For report: Loop on revisions (e.g., "Add more on diffusion?").
     - If no initial_prompt, start with human_tool for full details.
  2. Integrate: Use same agents/tools, but wrap in loops instead of sequential crew.
- **Testing**:
  - Add `tests/test_agents/test_interactive.py`: Mock human_tool with predefined responses (e.g., "adjust density to 0.4" → retry).
  - Manual: Simulate session, check retries work.

#### Step 7: Create Entry Point and Polish (Effort: Low, ~1 hour)
- **Goal**: User-friendly launch.
- **Actions**:
  1. In project root, add `main.py`:
     - Use argparse: `parser.add_argument('--mode', choices=['autonomous', 'interactive'], required=True)`
     - `parser.add_argument('prompt', nargs='?')` (optional for interactive).
     - Call `run_autonomous` or `run_interactive`.
  2. Update README.md: Add usage examples for both modes, prompt formats.
  3. Logging: Ensure all agents use `logging.py` for debug/info.
  4. Error Handling: Global try/except, save partial results on failure.
- **Testing**: Run `python main.py --mode autonomous "simulate..."` and interactive variant.

#### Step 8: Full Testing and Demo Prep (Effort: Medium, ~1 day)
- **Actions**:
  1. Run end-to-end for all topologies in both modes.
  2. Check outputs: Workflow dirs, reports (e.g., tables: Params | Rg | lp | D).
  3. Edge Cases: Invalid prompt, no solvent, short runs.
  4. Update `doc/progress.md`: Mark Phase 3 complete.
- **Demo**: Prepare script for "AI-driven topology exploration" (e.g., interactive optimization of sp).

Follow this sequentially; if stuck, refer to CrewAI docs. Once done, we can iterate!