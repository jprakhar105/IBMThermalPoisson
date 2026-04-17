// =================================================
// File: matrixAssembly.h
// assemble FEM matrices
// =================================================

#ifndef MATRIXASSEMBLY_H
#define MATRIXASSEMBLY_H
#include <vector>
#include <functional>

#include "mesh.h"
#include "flowParameters.h"

// Assemble K_tau, M_u and f_g
void assemble_matrices(const Mesh &mesh,
	       		     Flow &flow);

// Apply Dirichlet boundary conditions
void apply_dirichlet( const Mesh &msh,
			    SparseRowMap &Amap,
		      	    std::vector<double> &b,
			    std::vector<std::vector<int>> &DirichletNodes,
        		    std::vector<std::vector<double>> &DirichletValues);

#endif
