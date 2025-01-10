# Molecular Dynamics Code

## Description

This project implements a Molecular Dynamics (MD) program to simulate the motion of atoms based on their initial positions and velocities. The program calculates interactions using the Lennard-Jones potential and employs the Verlet algorithm to update atomic positions and velocities over time.

## Directory Structure

The project follows the structure below:

- **`LICENSE`**: Contains the license under which the code is distributed.
- **`AUTHORS`**: Lists the contributors to this project.
- **`README.md`**: Provides an overview of the project and instructions.
- **`INSTALL.md`**: Includes detailed instructions for compiling and running the program.
- **`src/`**: Contains the source code files.
- **`tests/`**: Contains unit tests to validate the functionality of the code.
- **`doc/`** *(optional)*: Includes additional documentation if required.

## Key Features

- **Input Handling**: Reads atomic positions, masses, and other parameters from an input file.
- **Force Calculations**: Implements the Lennard-Jones potential to compute interatomic forces.
- **Energy Computations**: Calculates kinetic, potential, and total energies of the system.
- **Verlet Algorithm**: Updates atomic positions and velocities using the Verlet integration method.
- **Trajectory Output**: Outputs atomic trajectories in `xyz` format for visualization with tools like Jmol or Molden.

## Dependencies

- Fortran compiler.
- Math library (`libm` for C).
- An external visualization tool (e.g., Jmol or Molden) for viewing the trajectory output.

## Usage

1. Compile the program using the instructions in `INSTALL.md`.
2. Provide an input file (`inp.txt`) with atomic data and simulation parameters.
3. Run the executable to simulate molecular dynamics.
4. View the generated `xyz` trajectory file in a visualization tool.
