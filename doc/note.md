20250904: with JMC

polymer: star polymer, bruch polymer, linear, ring, dendrimer,
solvent: single bead solvent
thermostats: langevin, nose-hoover(NVT)

step1: take prompt to polymer strucuture:
        identify polymer type
        use corresponding tool to generate datafile
            for polymer connectivity/position
            for solvent, use packmol for solvent positions
step2: initialize/run simulation, datafile, thermostat, NPT vs NVT, (lammps script tool), can define output of the simulation, pressure/
step3: analyze dump file: Rg, P(q), Kuhn length, Diffusion constant, each will be a tool,
step4: summarize results, give a report


tech stack:
lammps for running simulation
packmol for generating solvent positions
crewai for multi agent



reference for lammps run: https://code.ornl.gov/jyw/cyber-training-summer-school/-/tree/master/2021?ref_type=heads


conference to consider:
https://callforabstracts.acs.org/acsspring2026/I&EC at the "Data Analytics and AI For Chemistry, Manufacturing, and Healthcare" for the ACS 2026 Spring meeting at Atlanta? Deadline to submit abstract is Sept. 29, 2025




idea 20250909:
to make our system non-trivial, we can consider:
agent 1: different polymer topologies, different solvent density,
agent 2: different thermostat: langivin vs nose-hoover, NVT vs NPT, and different interaction parameters (solvent-solvent solvent-polymer polymer-polymer)
agent 3: analysis tools, choose analysis based on request and need

human_tool? request further human input if needed
search agent with serper search tool?


