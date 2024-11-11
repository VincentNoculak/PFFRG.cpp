# PFFRG.cpp

The repository contains a C++ implementation of the PFFRG method (PFFRG.cpp), as it is formulated in [Phys. Rev. B 109, 174414](https://doi.org/10.1103/PhysRevB.109.174414), for a Heisenberg model on the triangular lattice in a magnetic field, given by the Hamiltonian 

$\hat{\mathcal{H}} = \frac{1}{2} \sum_{ij} \sum_{\mu\nu} J^{\mu\nu}_{ij} \hat{S}^{\mu}_i \hat{S}^{\nu}_j - \sum_i \sum\_\mu h^{\mu}_i \hat{S}^{\mu}_i$. 

The code is able to treat perturbative seed fields that break any continuous spin rotation symmetries and result in three symmetry-inequivalent sublattices. These fields regularize PFFRG flow breakdowns and allow the computation of the $T=0$ magnetization curve of the XXZ model ($J^{xx}=J^{yy}_{ij} \neq J^{zz}$)

"FlowEquationTermGenerator.cpp" provides a C++ code that generates the terms of the flow equation for the two-particle vertex, given by Eq. (35) in [Phys. Rev. B 109, 174414](https://doi.org/10.1103/PhysRevB.109.174414), as output.
Similarly, "SpinCorrelationTermGenerator.cpp" generates terms of two-spin correlations expressed via pseudo-fermion vertex functions, given by Eq. (43) in [Phys. Rev. B 109, 174414](https://doi.org/10.1103/PhysRevB.109.174414), as output.
The generated terms of "FlowEquationTermGenerator.cpp" and "SpinCorrelationTermGenerator.cpp" are inserted into the source code of "PFFRG.cpp".

"PFFRG.cpp" can be compiled by using the command "g++ -O2 -fopenmp -o PFFRG PFFRG.cpp -lgsl".

"jobScript.sh" contains an example of a job script that can be used to run the compiled "PFFRG.cpp" code via the slurm workload manager.
