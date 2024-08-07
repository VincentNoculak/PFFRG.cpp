# PFFRG.cpp

The repository contains a C++ implementation of the PFFRG method "PFFRG.cpp", as it is formulated in [Phys. Rev. B 109, 174414](https://doi.org/10.1103/PhysRevB.109.174414), for a Heisenberg model on the triangular lattice in a magnetic field, given by the Hamiltonian 

$\hat{\mathcal{H}} = \frac{1}{2} \sum_{ij} \sum_{\mu\nu} J^{\mu\nu}_{ij} \hat{S}^{\mu}_i \hat{S}^{\nu}_j - \sum_i \sum\_\mu h^{\mu}_i \hat{S}^{\mu}_i$. 

The code applies perturbative seed fields that break any continuous spin rotation symmetries and result in three symmetry-inequivalent sublattices. These fields regularize PFFRG flow breakdowns and allow to compute the $T=0$ magnetization curve of the XXZ model ($J^{xx}=J^{yy}_{ij} \neq J^{zz}$)

"FlowEquationTermGenerator.cpp" provides a C++ code that generates the terms of the flow equation for the two-particle vertex.
Similarly, "SpinCorrelationTermGenerator.cpp" generates terms of two-spin correlations expressed via vertex functions.
The generated terms of "FlowEquationTermGenerator.cpp" and "SpinCorrelationTermGenerator.cpp" are inserted into the source code of "PFFRG.cpp".

"PFFRG.cpp" can be compiled by using the command "g++ -O2 -fopenmp -o PFFRG PFFRG.cpp -lgsl".
