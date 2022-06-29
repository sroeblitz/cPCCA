README

The robust Perron Cluster Analysis (PCCA+) has become a popular algorithm
for coarse-graining transition matrices of nearly decomposable Markov
chains with transition states. Such matrices arise, for example, from molecular
dynamics simulations or gene-regulatory networks, and they have been the
target of Markov state model building. Though originally developed for reversible
Markov chains, it has been shown that PCCA+ can also be applied to cluster
non-reversible Markov chains. However, the algorithm was implemented by
assuming the dominant (target) eigenvalues are real numbers.
To overcome this limitation, cPCCA+ separates the real and imaginary parts
of complex eigenvectors before handing them over to the PCCA+ algorithm. 
This makes PCCA+ applicable to transition matrices that have complex
eigenvalues, including matrices with a circular transition patterns.

# General information

The code implements the cPCCA+ method and applies it to a number of example
matrices that are presented in Examples 1 and 2 in the paper 

"Robust Perron Cluster analysis for coarse-graining of
non-reversible stochastic matrices with complex eigenvalues" 
by Anna-Simone Frank, Alexander Sikorski and Susanna Röblitz

Contributer: Anna-Simone Frank, Susanna Röblitz

Maintainer: Susanna Röblitz

# Included code files  and their description


|File names |			Description |
|---------------|-----------------------------------------------------------|
|`main_cpcca_Exp?.m` |  Main file for reproducing the results presented in Example? in the manuscript|
|`Exp1.mat`    |        Input matrix for `main_cpcca_Exp1.m`|
|`Exp2_x01.mat`   |     Input matrix for `main_cpcca_Exp2iii.m` and`main_cpcca_Exp2iv.m` |
|`Exp2_x09.mat`  |      Input matrix for `main_cpcca_Exp2i.m` and `main_cpcca_Exp2ii.m` |
|`compute_subspace.m` |  Computes the eigenvectors of the input matrix for a given target eigenvalue |
| `preprocessEVS.m`  |  Separates real and imaginary parts of complex eigenvectors |
| `pcca.m`		|	 	Main file for the pcca+ analysis. It depends on  the following files: </br> `fillA.m` </br> `indexsearch.m` </br> `objective.m` </br> `orthogon.m` </br> `opt_soft.m` |
| `orthogon.m`     |     Orthonormalization of eigenvectors w.r.t. a given density|
| `opt_soft.m`     |     Call to optimization routine for finding the optimal transformation matrix A |
|`indexsearch.m`    |   Calculation of initial guess by inner simplex algorithm |
|`objective.m` |        Objective function for the optimization |
|`fillA.m` |            Function for reconstruction of A during optimization |
| `main_nlscon.m`  </br> `problem_pcca-nlscon.m` </br> `nlscon.m` | Subroutines for Gauss-Newton method |

# Code outputs

1. List of states and the cluster they belong to with largest membership
2. Coarse-grained transition matrix	
3. Plot of membership vectors (Figure 1)


# How to run the code

Run one of the files main_cpcca_Exp?.m

