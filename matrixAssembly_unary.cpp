// ====================================================
// File: matrixAssembly.cpp (Q4 Version - Corrected)
// ====================================================
#include <cmath>
#include <array>
#include <vector>
#include <iostream>
#include "matrixAssembly_unary.h"

// GLOBAL VARIABLES FOR QUADRATURE
std::vector<std::vector<double>> gaussPts;
std::vector<double> gaussWts;

// 1. Bilinear Shape Functions 
// N = 0.25 * (1 +/- xi) * (1 +/- eta)
static inline void shapeFunction(double xi, double eta, int Ordr,
             std::vector<double>& N,
             std::vector<std::array<double,2>>& dN)
{
    if (N.size() != 4) N.resize(4);
    if (dN.size() != 4) dN.resize(4);

    // Node 1: Bottom-Left (-1, -1)
    N[0] = 0.25 * (1 - xi) * (1 - eta);
    // Node 2: Bottom-Right (1, -1)
    N[1] = 0.25 * (1 + xi) * (1 - eta);
    // Node 3: Top-Right (1, 1)
    N[2] = 0.25 * (1 + xi) * (1 + eta);
    // Node 4: Top-Left (-1, 1)
    N[3] = 0.25 * (1 - xi) * (1 + eta);

    // Derivatives dN/dxi (index 0), dN/deta (index 1)
    // Node 1
    dN[0][0] = -0.25*(1-eta); 
    dN[0][1] = -0.25*(1-xi);

    // Node 2
    dN[1][0] =  0.25*(1-eta); 
    dN[1][1] = -0.25*(1+xi);

    // Node 3
    dN[2][0] =  0.25*(1+eta); 
    dN[2][1] =  0.25*(1+xi);

    // Node 4
    dN[3][0] = -0.25*(1+eta); 
    dN[3][1] =  0.25*(1-xi);
}

// 2. Gauss Quadrature (2x2 Rule)
static inline void gaussQuad(int n) {
    gaussPts.resize(4);
    gaussWts.resize(4);
    
    // 1/sqrt(3) is the sampling point for 2-point Gaussian rule
    double pt = 0.57735026919; 

    gaussPts[0] = {-pt, -pt}; 
    gaussPts[1] = { pt, -pt}; 
    gaussPts[2] = { pt,  pt}; 
    gaussPts[3] = {-pt,  pt}; 
    
    // Weights are 1.0 for all points in 2x2 rule
    gaussWts = {1.0, 1.0, 1.0, 1.0};
}

void assemble_matrices(const Mesh &msh, Flow &flow)
{
    double ReInv = 1.0/flow.Re;
    int n = msh.n_nodes;
    
    // Initialize Global Arrays
    flow.Ktau.resize(n, std::vector<double>(n, 0.0));
    flow.fg.resize(n, std::vector<double>(2, 0.0));
    
    SparseRowMap Kmap(n);
    flow.K_tau = CSR(n);

    // Initialize Gauss Points
    gaussQuad(2); 

    // Loop over Elements
    for (int e=0; e<msh.n_elements; ++e) {
        
        auto conn = msh.Elements[e];
        int nElemNodes = 4; // Q4 Fixed

        // Get Element Coordinates
        std::vector<std::pair<double,double>> xy(nElemNodes);
        for (int a=0; a<nElemNodes; ++a) {
            xy[a] = msh.Coords[conn[a]];
        }

        // Element Matrices
        std::vector<std::vector<double>> Ke_tau(nElemNodes, std::vector<double>(nElemNodes, 0.0));
        std::vector<std::vector<double>> fe_g(nElemNodes, std::vector<double>(2, 0.0));

        // Loop over Gauss Points
        for (size_t gp=0; gp<gaussPts.size(); ++gp) {
            double xi = gaussPts[gp][0];
            double eta = gaussPts[gp][1];

            std::vector<double> N(nElemNodes);
            std::vector<std::array<double,2>> dN(nElemNodes);
            shapeFunction(xi, eta, 1, N, dN);

            // Calculate Jacobian
            double J[2][2] = {{0}, {0}};
            double x_phys = 0, y_phys = 0;
            
            for (int a=0; a<nElemNodes; ++a) {
                J[0][0] += dN[a][0] * xy[a].first;  // dx/dxi
                J[0][1] += dN[a][0] * xy[a].second; // dy/dxi
                J[1][0] += dN[a][1] * xy[a].first;  // dx/deta
                J[1][1] += dN[a][1] * xy[a].second; // dy/deta
                
                x_phys += N[a] * xy[a].first;
                y_phys += N[a] * xy[a].second;
            }

            double detJ = J[0][0]*J[1][1] - J[0][1]*J[1][0];
            
            // Inverse Jacobian
            double invJ[2][2] = {{ J[1][1]/detJ, -J[0][1]/detJ}, 
                                 {-J[1][0]/detJ,  J[0][0]/detJ}};

            // Convert local derivatives to physical derivatives (dN/dx, dN/dy)
            std::vector<std::array<double, 2>> dNdx(nElemNodes);
            for (int a=0; a<nElemNodes; ++a) {
                // FIXED: Properly aligning chain rule for physical derivatives
                dNdx[a][0] = invJ[0][0]*dN[a][0] + invJ[0][1]*dN[a][1]; // dN/dx
                dNdx[a][1] = invJ[1][0]*dN[a][0] + invJ[1][1]*dN[a][1]; // dN/dy
            }

            double wDetJ = gaussWts[gp] * detJ;
            double src = flow.g(x_phys, y_phys, 0); 

            // Assemble Element Matrix
            for (int i=0; i<nElemNodes; ++i) {
                for (int j=0; j<nElemNodes; ++j) {
                    double val = ReInv * (dNdx[i][0]*dNdx[j][0] + dNdx[i][1]*dNdx[j][1]);
                    Ke_tau[i][j] += val * wDetJ;
                }
                fe_g[i][0] += N[i] * src * wDetJ;
                fe_g[i][1] += N[i] * src * wDetJ;
            }
        }

        // Scatter to Global Matrix
        for (int i=0; i<nElemNodes; ++i) {
            int gi = conn[i];
            for (int j=0; j<nElemNodes; ++j) {
                int gj = conn[j];
                double val = Ke_tau[i][j];
                flow.Ktau[gi][gj] += val;
                if (std::abs(val) > 1e-16) {
                    Kmap.add(gi, gj, val);
                }
            }
            flow.fg[gi][0] += fe_g[i][0];
            flow.fg[gi][1] += fe_g[i][1];
        }
    }
    
    flow.K_tau = convert_to_csr(Kmap);
}

// Apply Dirichlet BCs (Unchanged)
// =========================================================
// MATLAB-style Symmetric Dirichlet Boundary Application
// =========================================================
void apply_dirichlet(const Mesh &msh, SparseRowMap &K, std::vector<double> &F,
                     std::vector<std::vector<int>> &DirichletNodes,
                     std::vector<std::vector<double>> &DirichletValues)
{
    // First, flatten the DOF structures into 1D arrays just like MATLAB
    std::vector<int> bc_dofs;
    std::vector<double> bc_vals;

    int N = msh.n_nodes; // Total number of spatial nodes

    // Assuming 2 DOFs per node (T is dir=0, vector is dir=1)
    for (int dir = 0; dir < 2; ++dir) {
        for (size_t j = 0; j < DirichletNodes[dir].size(); ++j) {
            int global_dof = DirichletNodes[dir][j] + N * dir;
            bc_dofs.push_back(global_dof);
            bc_vals.push_back(DirichletValues[dir][j]);
        }
    }

    int n_bc = bc_dofs.size();

    // 1. Modify global force vector: F = F - K(:, bc_dofs) * bc_vals
    // 2. Zero out the columns: K(:, bc_dofs) = 0
    for (int k = 0; k < n_bc; ++k) {
        int dof = bc_dofs[k];
        double val = bc_vals[k];

        for (int i = 0; i < K.n; ++i) {
            auto it = K.rows[i].find(dof); // Check if column 'dof' exists in row 'i'
            if (it != K.rows[i].end()) {
                // Subtract known value from force vector
                F[i] -= it->second * val; 
                
                // Erase the column entry to maintain matrix symmetry
                K.rows[i].erase(it); 
            }
        }
    }

    // 3. Overwrite Force vector: F(bc_dofs) = bc_vals
    // 4. Zero out the rows: K(bc_dofs, :) = 0
    // 5. Set diagonal to 1: K = K + K_bc
    for (int k = 0; k < n_bc; ++k) {
        int dof = bc_dofs[k];
        double val = bc_vals[k];

        F[dof] = val;             // Step 3
        K.rows[dof].clear();      // Step 4 (Clears all column entries in this row)
        K.rows[dof][dof] = 1.0;   // Step 5
    }
}