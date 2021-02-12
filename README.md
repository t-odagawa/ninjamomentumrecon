# Tools for momentum reconstruction evaluation in NINJA MC (v0.0.1)

This is a set of tools for evaluation of momentum reconstruction
in NINJA ECC using Monte Carlo simulated data.

## External dependencies (same as WAGASCI+BabyMIND Geant4 simulation)

Below is a copy of README.md in the MC simulation software written by Giorgio Pintaudi.

### GEANT4

[GEANT4](http://heant4.cern/ch/) is a toolkit for the simulation of
the passage of particles through matter developed by CERN.
Geant4 v10.5.0+ is recommended.

### ROOT

[ROOT](https://root.cern.ch/) is an object-oriented program and
library developed by CERN. ROOT 6.20+ is recommended.

### BOOST

[BOOST](https://www.boost.org/) is a set of libraries for the C++
programming language that provides support for tasks and structures
such as linear algebra, pseudorandom number generation,
multithreading, image processing, regular expressions, and unit
testing. Boost 1.53+ is recommended.

### Wagasci BabyMIND Monte Carlo
[Wagasci BabyMIND Monte Carlo](https://git.t2k.org/wagasci_babymind/wagasci-babymind-monte-carlo) is a developed Monte Carlo simulation software for the WAGASCI-BabyMIND experiment.
It includes necessary libraries for the analyses.
Wagasci BabyMIND Monte Carlo v0.1.12+ is recommended.

### Wagasci BabyMIND event display
[Wagasci BabyMIND event display](https://git.t2k.org/wagasci_babymind/wagasci-babymind-event-display) is an event display tool for the WAGASCI-BabyMIND experiment.
Wagasci BabyMIND event display v0.1.0+ is recommended.

### Wagasce Reconstruction

## Installation

Only the following operative systems are tested and supported:
 - CentOS 7 (KEKCC)

To install the software with all its dependencies, use a bash shell script `install.bash` in Wagasci BabyMIND Monte Carlo.
For further info about the software installation, refer to the [WAGASCI BabyMIND Monte Carlo gitlab](https://git.t2k.org/wagasci_babymind/wagasci-babymind-monte-carlo).