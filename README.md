# PFFRG.cpp

The repository contains a C++ implementation of the PFFRG method "PFFRG.cpp" for a Heisenberg model on the triangular lattice in a magnetic field. Seed fields are applied that break any continuous spin rotation symmetries and result in three symmetry-inequivalent sublattices.

In addition "FlowEquationTermGenerator.cpp" provides a C++ code that generates the terms of the flow equation for the two-particle vertex.
Similarly, "SpinCorrelationTermGenerator.cpp" generates terms of the susceptibility expressed via vertex functions.
The generated terms of "FlowEquationTermGenerator.cpp" and "SpinCorrelationTermGenerator.cpp" are inserted into the source code of "PFFRG.cpp".

"PFFRG.cpp" can be compiled by using the command "g++ -O2 -fopenmp -o PFFRG PFFRG.cpp -lgsl".
