# Closed-Form Global Optimization of Beyond Diagonal Reconfigurable Intelligent Surfaces

This code package is related to the paper:

M. Nerini, S. Shen and B. Clerckx, "[Closed-Form Global Optimization of Beyond Diagonal Reconfigurable Intelligent Surfaces](https://ieeexplore.ieee.org/document/10155675)," in IEEE Transactions on Wireless Communications, 2023.

If you use this code or any modified part of it, please cite our paper.

## Abstract

Reconfigurable intelligent surfaces (RISs) allow controlling the propagation environment in wireless networks by tuning multiple reflecting elements. RISs have been traditionally realized through single connected architectures, mathematically characterized by a diagonal scattering matrix. Recently, beyond diagonal RISs (BD-RISs) have been proposed as a novel branch of RISs whose scattering matrix is not limited to be diagonal, which creates new benefits and opportunities for RISs. Efficient BD-RIS architectures have been realized based on group and fully connected reconfigurable impedance networks. However, a closed-form solution for the global optimal scattering matrix of these architectures is not yet available. In this paper, we provide such a closed-form solution proving that the theoretical performance upper bounds can be exactly achieved for any channel realization. We first consider the received signal power maximization in single-user single-input single-output (SISO) systems aided by a BD-RIS working in reflective or transmissive mode. Then, we extend our solution to single-user multiple-input multiple-output (MIMO) and multi-user multiple-input single-output (MISO) systems. We show that our algorithm is less complex than the iterative optimization algorithms employed in the previous literature. The complexity of our algorithm grows linearly (resp. cubically) with the number of RIS elements in the case of group (resp. fully) connected architectures.

## Content of Code Package

The file `main_Fig*i*.m` returns the results needed to reproduce Fig. *i* in the paper, for *i* = 2, ..., 7.
