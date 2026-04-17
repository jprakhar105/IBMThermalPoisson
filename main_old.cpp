// ========================================
// File: main.cpp
// ========================================
#include <iostream>

#include "mesh.h"
#include "flowParameters.h"
#include "matrixAssembly.h"
#include "writeOutput.h"

int main(){
	std::cout << "Starting SATEEK Solver" << "\n";

	// ----------------------------------------
	// Simulation settings
	// ----------------------------------------
	Flow flow;
	
	// Setting Flow Parameters
	flow.setParameters();


	// ----------------------------------------
	// Mesh Settings
	// ----------------------------------------
	Mesh msh, msh1;
	
	
	// define domain extents
	msh.X[0] = -1;	// X domain start
	msh.X[1] =  1;	// X domain end
	msh.Y[0] = -1;	// Y domain start
	msh.Y[1] =  1;	// Y domain end
	// Mesh Resolution (Number of Elements)
	msh.nx = 10;	// X-direction
	msh.ny = 10;	// Y-direction
	// Number of Edges
	msh.n_edges = 4;
	// Element Order (1-Linear, 2-Quadratic)
	msh.Ordr = 1;	// Order of FEM elements
	


	// ----------------------------------------
	// Generate Mesh
	// ----------------------------------------
	msh.generate();

	msh.checkNormals();

	int N = msh.n_nodes;

	// ----------------------------------------
	// Defining Dirichlet BCs
	// ----------------------------------------
	flow.DirichletEdges.resize(2);

	flow.DirichletEdges[0] = {0,1,2,3};	// in x-direction
	flow.DirichletEdges[1] = {0,1,2,3};	// in y-direction

	flow.calculateDirichletNodes(msh);
	flow.calculateDirichletValues(msh);



	// ----------------------------------------
	// Assemble Static FEM Matrices
	// ----------------------------------------

	assemble_matrices(msh, flow);


	/*	
	// print matrix
	for (int i=0; i<msh.n_nodes; ++i) {
		for (int j=0; j<msh.n_nodes; ++j) {
			std::cout << flow.K_tau.rows[i][j] << "\t";
		}
		std::cout << "\n";
	}
	std::cout << "\n";
	*/

	/*
	// print matrix
        for (int i=0; i<msh.n_nodes; ++i) {
                for (int j=0; j<2; ++j) {
                        std::cout << flow.fg[i][j] << "\t";
                }
                std::cout << "\n";
	}
	*/	


	// ----------------------------------------
        // Initialize Flow Field
        // ----------------------------------------
	flow.Ui.resize(N,std::vector<double>(2));

        int init = 1;           // 0 - given value, 1 - calculate steady state field
        flow.initializeFlow(msh, init);


	// Obtain Analytical Solution
	flow.UA.resize(N,std::vector<double>(2));
	flow.calculateAnalyticalSolution(msh);


	std::vector<std::vector<double>> UErr;
	UErr.resize(N,std::vector<double>(2));
	for (int i=0; i<msh.n_nodes; ++i) {
		for (int dir=0; dir<2; ++dir) {
			UErr[i][dir] = flow.Ui[i][dir] - flow.UA[i][dir];
		}
	}


	char fname[256];
        std::sprintf(fname, "solution_%04d.dat", 0);
        write_output(fname, msh, flow.Ui, flow.UA, UErr);

	double max = 0;
	for (int i=0; i<msh.n_nodes; ++i) {
		for (int dir=0; dir<2; ++dir) {
			if (UErr[i][dir]>max) {
				max = UErr[i][dir];
			}
		}
	}

	std::cout << "L-inf Norm = " << max << "\n";

	return 0;
}
