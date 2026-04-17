// ==============================================
// File: flowFunctions.cpp
// 1. Initialize Flow Field
// ==============================================

#include <iostream>
#include <vector>
#include <algorithm>
#include <numbers>
#include <cmath>

#include "flowParameters.h"
#include "matrixAssembly.h"
#include "linSolve.h"
#include "writeOutput.h"

double Flow::g(double x,double y, int dir){
                double gval = (2*pi*pi) * sin(pi*x) * sin(pi*y);

                return gval;
}


void Flow::setParameters() {
	// Global constants
	pi = 3.141;

	// Fluid properties
        Re = 1;

        // Time stepping details
        dt = 0.01;
        nSteps = 1000;
        theta = 0.5;
}


void Flow::calculateDirichletNodes(const Mesh &msh) {

	DirichletNodes.resize(2);
	DirichletValues.resize(2);
	nDirichletNodes.resize(2);

	for (int dir=0; dir<2; ++dir) {
		int nEdge = DirichletEdges[dir].size();
		int nElems = 0;
		int k=0;	
		for (int iEdge=0; iEdge<nEdge; ++iEdge) {
			int EdgeID = DirichletEdges[dir][iEdge];
			nElems += msh.EdgeNodeConnectivity[EdgeID].size();
			DirichletNodes[dir].resize(nElems*msh.Ordr+1, -1);
			for (const auto& EdgeCon : msh.EdgeNodeConnectivity[iEdge]) {
				for (int i=0; i<msh.Ordr+1; ++i) {
					auto it = std::find(DirichletNodes[dir].begin(), DirichletNodes[dir].end(), EdgeCon[i]);
					if (it != DirichletNodes[dir].end()) {
					}
					else {
						DirichletNodes[dir][k] = EdgeCon[i];
						++k;
					}
				}
			}
		}
		nDirichletNodes[dir] = k;
		DirichletNodes[dir].resize(nDirichletNodes[dir]);
		DirichletValues[dir].resize(nDirichletNodes[dir]);

		/*	
		std::cout << "\n" << "Dirichlet node list" << "\n";
		for (int i=0; i<k; ++i) {
			std::cout << DirichletNodes[dir][i] << "\t";
		}
		std::cout << "\n";
		*/
		
	}
	
	
}


void Flow::calculateAnalyticalSolution(const Mesh &msh) {
	int N = msh.n_nodes;

	for (int iNode=0; iNode<N; ++iNode) {
		for (int dir=0; dir<2; ++dir) {
			UA[iNode][dir] = sin(pi * msh.Coords[iNode].first) * sin(pi * msh.Coords[iNode].second);
		}
	}
}


void Flow::calculateDirichletValues(const Mesh &msh) {

	for (int dir=0; dir<2; ++dir) {
		for (int i=0; i<nDirichletNodes[dir]; ++i) {
			DirichletValues[dir][i] = 0;	// Homogeneous boundary
		}
	}
}


void Flow::initializeFlow(const Mesh &msh,	
			  	int init) {

	int N = msh.n_nodes;
	int NEdges = msh.n_edges;

	if (init==0) {

		// uniform zero velocity in whole domain
		for (int i=1; i<N; ++i) {
			for (int j=0; j<2; ++j) {
				Ui[i][j] = 0;
			}
		}
	}
	else {
		// Solve for steady state
		// Ktau*U_i = fg_i
		// Ax = b system

			
		// A: Ktau
		SparseRowMap Amap;
		Amap = SparseRowMap(2*N);
		for (int i=0; i<N; ++i) {
			// Copy K_tau
			for (int k=K_tau.row_ptr[i]; k<K_tau.row_ptr[i+1]; ++k) {
				Amap.add(i  , K_tau.col_idx[k]  , K_tau.values[k]);
				Amap.add(i+N, K_tau.col_idx[k]+N, K_tau.values[k]);
			}
		}
		
		
		// calculating RHS
		// b: fg
		std::vector<double> b;
        	b.resize(2*N);
                for (int i=0; i<N; ++i) {
                       	// Copy fg
                        b[i  ] += fg[i][0];
			b[i+N] += fg[i][1];
                }


				
		// Apply Dirichlet BCs
		apply_dirichlet(msh,Amap,b,DirichletNodes,DirichletValues);


		/*
		// Display Global Matrix
		for (int i=0; i<2*N; ++i) {
			for (int j=0; j<2*N; ++j) {
				std::cout << Amap.rows[i][j] << "\t";
			}
			std::cout << "\n";
		}
		std::cout << "\n";
		*/
		
		// initialize unknown vector
		std::vector<double> x;
                x.resize(2*N);

		//Solve Ax=b system

		// Gauss Elimination
        	//GaussElimination(Amap, b, x);	

		
		// Conjugate Gradient Methods
		double tol = 1e-12;
		int maxIter = 10000;
		CSR A = convert_to_csr(Amap);
		// Deallocate SparseRowMap
		Amap.rows.clear();
		Amap.rows.shrink_to_fit();
		//CG
		//int iter = ConjugateGradient(A,b,x,tol,maxIter);
		//BiCGStab
		int iter = BiCGSTAB(A,b,x,tol,maxIter);



		std::cout << "Solution obtained in " << iter << " Iterations" << "\n";
		

		// Extract Output Vectors
		for (int iNode=1; iNode<N; ++iNode) {
                        for (int dir=0; dir<2; ++dir) {
                                Ui[iNode][dir] = x[iNode + dir*N];
                        }
                }

	}

}
