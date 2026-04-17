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

// Define PI if not available
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// ------------------------------------------------
// 1. Source Term (qs = 1)
// Equation: Grad^2 T + qs = 0  =>  -Grad^2 T = qs
// Weak form R.H.S integrates this function.
// ------------------------------------------------
double Flow::g(double x, double y, int dir) {
    return 1.0; 
}

void Flow::setParameters() {
    pi = M_PI;
    Re = 1.0; // Diffusion coefficient 1.0
    
    // Iterative solver settings
    dt = 0.01;
    nSteps = 1000;
    theta = 0.5;
}

// ------------------------------------------------
// 2. Dirichlet Nodes (Unchanged logic)
// ------------------------------------------------
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
    }
}

// ------------------------------------------------
// 3. Analytical Solution for Square Domain
// Solves: Grad^2 T = -1, T=0 on boundaries [-1,1]x[-1,1]
// Method: Fourier Series Expansion
// ------------------------------------------------
void Flow::calculateAnalyticalSolution(const Mesh &msh) {
    int N = msh.n_nodes;
    double a = 1.0; // Half-width of the square (Domain is -1 to 1)

    for (int iNode=0; iNode<N; ++iNode) {
        double x = msh.Coords[iNode].first;
        double y = msh.Coords[iNode].second;
        double sum = 0.0;
        
        // Truncate series at 50 terms for high precision
        for (int n = 0; n < 50; ++n) {
            double term_k = 2.0 * n + 1.0;
            double constant = std::pow(-1.0, n) / std::pow(term_k, 3.0);
            
            double cos_term = std::cos((term_k * M_PI * x) / (2.0 * a));
            
            double cosh_num = std::cosh((term_k * M_PI * y) / (2.0 * a));
            double cosh_den = std::cosh((term_k * M_PI * a) / (2.0 * a)); // cosh((2n+1)pi/2)
            
            sum += constant * (1.0 - (cosh_num / cosh_den)) * cos_term;
        }

        double val = (16.0 * a * a / std::pow(M_PI, 3.0)) * sum;

        // Apply to both DOF components (though usually scalar T is 1 DOF, code uses 2)
        for (int dir=0; dir<2; ++dir) {
            UA[iNode][dir] = val;
        }
    }
}

void Flow::calculateDirichletValues(const Mesh &msh) {
    for (int dir=0; dir<2; ++dir) {
        for (int i=0; i<nDirichletNodes[dir]; ++i) {
            DirichletValues[dir][i] = 0; // T = 0 at boundaries
        }
    }
}

void Flow::initializeFlow(const Mesh &msh, int init) {
    int N = msh.n_nodes;

    if (init==0) {
        for (int i=1; i<N; ++i) {
            for (int j=0; j<2; ++j) Ui[i][j] = 0;
        }
    }
    else {
        // Assemble Global Matrix
        SparseRowMap Amap(2*N);
        for (int i=0; i<N; ++i) {
            for (int k=K_tau.row_ptr[i]; k<K_tau.row_ptr[i+1]; ++k) {
                Amap.add(i, K_tau.col_idx[k], K_tau.values[k]);
                Amap.add(i+N, K_tau.col_idx[k]+N, K_tau.values[k]);
            }
        }
        
        // Assemble RHS
        std::vector<double> b(2*N, 0.0);
        for (int i=0; i<N; ++i) {
            b[i]   += fg[i][0];
            b[i+N] += fg[i][1];
        }

        // Apply BCs
        apply_dirichlet(msh, Amap, b, DirichletNodes, DirichletValues);

        // Solve
        std::vector<double> x(2*N);
        double tol = 1e-12;
        int maxIter = 10000;
        CSR A = convert_to_csr(Amap);
        Amap.rows.clear(); 
        
        // Use BiCGSTAB
        int iter = BiCGSTAB(A, b, x, tol, maxIter);
        
        // Store Result
        for (int iNode=0; iNode<N; ++iNode) {
            for (int dir=0; dir<2; ++dir) {
                Ui[iNode][dir] = x[iNode + dir*N];
            }
        }
    }
}