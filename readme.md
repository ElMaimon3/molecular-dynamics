This simulation is meant for simplified rvc files of proteins, and produces another rvc file, which can be converted to ther formats for viewing
# How to Run the Simulation
To run the simulation, use the following command:
python MDSimulation.py <input_file> <output_file> <kB> <kN> <nbCutoff> <dt> <mass>
<input_file>: The input file containing the initial coordinates, velocities, and connectivity information of the atoms in the system.
<output_file>: The output file to store the simulation results, including the energies and atom coordinates.
<kB>: Bonded force constant.
<kN>: Non-bonded force constant.
<nbCutoff>: Non-bonded cutoff distance.
<dt>: Time step size.
<mass>: Mass of the atoms.
For example, to run the simulation at default conditions, use the following command:
python MDSimulation.py 6pti-300K.rvc default.out 40000.0 400.0 0.50 0.0010 12.0
# Pseudocode
The following is the pseudocode of the main algorithm used in the simulation:

Read input_file
Create atom objects and store their coordinates, velocities, and connectivity information

For each step in the simulation:
    Update atom positions using Verlet integration (first half-step)
    Calculate forces and potential energies
    Update atom velocities using Verlet integration (second half-step)
    Write energies to output_file
    
    If step is a multiple of 10:
        Write atom coordinates to a separate .rvc file
    
Close output_file