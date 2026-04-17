// ==============================================
// File: flowFunctions.cpp 
// ==============================================

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>

#include "flowParameters.h"
#include "matrixAssembly.h"
#include "linSolve.h"
#include "writeOutput.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Source Term g = 1.0
double Flow::g(double x, double y, int dir) {
    return 1.0; 
}

void Flow::setParameters() {
    pi = M_PI;
    Re = 1.0;
    dt = 0.01;
    nSteps = 1000;
    theta = 0.5;
}

// -------------------------------------------------------------
// 1. ROBUST DIRICHLET DETECTION (DYNAMIC EDGE BASED)
// -------------------------------------------------------------
void Flow::calculateDirichletNodes(const Mesh &msh) {
    DirichletNodes.resize(2);
    DirichletValues.resize(2);
    nDirichletNodes.resize(2);

    // Clear previous lists
    DirichletNodes[0].clear();
    DirichletNodes[1].clear();

    double tol = 1e-6; // Tolerance for floating point comparison

    // Loop over the edges requested in main.cpp for Variable 0 (Temperature)
    for (int edge : DirichletEdges[0]) {
        for (int i = 0; i < msh.n_nodes; ++i) {
            double x = msh.Coords[i].first;
            double y = msh.Coords[i].second;

            bool is_on_edge = false;
            // Map Edge IDs to Physical Coordinates
            if (edge == 0 && std::abs(y - (-1.0)) < tol) is_on_edge = true;      // South / Bottom
            else if (edge == 1 && std::abs(x - 1.0) < tol) is_on_edge = true;    // East / Right
            else if (edge == 2 && std::abs(y - 1.0) < tol) is_on_edge = true;    // North / Top
            else if (edge == 3 && std::abs(x - (-1.0)) < tol) is_on_edge = true; // West / Left

            if (is_on_edge) {
                if (std::find(DirichletNodes[0].begin(), DirichletNodes[0].end(), i) == DirichletNodes[0].end()) {
                    DirichletNodes[0].push_back(i);
                    DirichletNodes[1].push_back(i);
                }
            }
        }
    }

    // Finalize sizes and set values to 0.0
    nDirichletNodes[0] = DirichletNodes[0].size();
    DirichletValues[0].assign(nDirichletNodes[0], 0.0); // T = 0
    
    nDirichletNodes[1] = DirichletNodes[1].size();
    DirichletValues[1].assign(nDirichletNodes[1], 0.0);
    
    if (msh.nx == 5) { 
         std::cout << "Dirichlet Nodes Detected: " << nDirichletNodes[0] << "\n";
    }
}

// -------------------------------------------------------------
// 2. ANALYTICAL SOLUTION (Fallback)
// -------------------------------------------------------------
void Flow::calculateAnalyticalSolution(const Mesh &msh) {
    int N = msh.n_nodes;
    UA.resize(N, std::vector<double>(2));

    for (int iNode=0; iNode<N; ++iNode) {
        double x = msh.Coords[iNode].first;
        double val = 0.5 * (1.0 - x*x);
        UA[iNode][0] = val;
        UA[iNode][1] = 0.0; // Unused variable
    }
}

void Flow::calculateDirichletValues(const Mesh &msh) {
    for (int dir=0; dir<2; ++dir) {
        for (int i=0; i<nDirichletNodes[dir]; ++i) DirichletValues[dir][i] = 0;
    }
}

void Flow::initializeFlow(const Mesh &msh, int init) {
    int N = msh.n_nodes;
    if (init==0) {
        for (int i=1; i<N; ++i) 
            for (int j=0; j<2; ++j) Ui[i][j] = 0;
    }
    else {
        SparseRowMap Amap(2*N);
        std::vector<double> b(2*N, 0.0);

        for (int i=0; i<N; ++i) {
            for (int k=K_tau.row_ptr[i]; k<K_tau.row_ptr[i+1]; ++k) {
                Amap.add(i, K_tau.col_idx[k], K_tau.values[k]);
                Amap.add(i+N, K_tau.col_idx[k]+N, K_tau.values[k]);
            }
            b[i]   += fg[i][0];
            b[i+N] += fg[i][1];
        }

        apply_dirichlet(msh, Amap, b, DirichletNodes, DirichletValues);

        std::vector<double> x(2*N);
        CSR A = convert_to_csr(Amap);
        Amap.rows.clear(); 
        
        BiCGSTAB(A, b, x, 1e-12, 10000);
        
        for (int iNode=0; iNode<N; ++iNode) {
            for (int dir=0; dir<2; ++dir) Ui[iNode][dir] = x[iNode + dir*N];
        }
    }
}