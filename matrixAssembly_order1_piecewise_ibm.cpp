// ====================================================
// File: matrixAssembly.cpp 
// 2D Gauss Quadrature & Piecewise IBM Allocation
// ====================================================
#include <cmath>
#include <array>
#include <vector>
#include <iostream>
#include <algorithm>
#include "matrixAssembly.h"

// GLOBAL VARIABLES FOR QUADRATURE
std::vector<std::vector<double>> gaussPts;
std::vector<double> gaussWts;

// 1. Bilinear Shape Functions 
static inline void shapeFunction(double xi, double eta, int Ordr,
             std::vector<double>& N,
             std::vector<std::array<double,2>>& dN)
{
    if (N.size() != 4) N.resize(4);
    if (dN.size() != 4) dN.resize(4);

    N[0] = 0.25 * (1 - xi) * (1 - eta);
    N[1] = 0.25 * (1 + xi) * (1 - eta);
    N[2] = 0.25 * (1 + xi) * (1 + eta);
    N[3] = 0.25 * (1 - xi) * (1 + eta);

    dN[0][0] = -0.25*(1-eta); dN[0][1] = -0.25*(1-xi);
    dN[1][0] =  0.25*(1-eta); dN[1][1] = -0.25*(1+xi);
    dN[2][0] =  0.25*(1+eta); dN[2][1] =  0.25*(1+xi);
    dN[3][0] = -0.25*(1+eta); dN[3][1] =  0.25*(1-xi);
}

static inline void gaussQuad(int n) {
    gaussPts.resize(4);
    gaussWts.resize(4);
    double pt = 0.57735026919; 
    gaussPts[0] = {-pt, -pt}; 
    gaussPts[1] = { pt, -pt}; 
    gaussPts[2] = { pt,  pt}; 
    gaussPts[3] = {-pt,  pt}; 
    gaussWts = {1.0, 1.0, 1.0, 1.0};
}

void solve_CG(const SparseRowMap &A, const std::vector<double> &b, std::vector<double> &x) {
    int n = b.size();
    x.assign(n, 0.0);
    
    // Set boundary nodes to their exact values to speed up convergence
    for (int i = 0; i < n; ++i) {
        if (A.rows[i].size() == 1 && A.rows[i].count(i)) {
            x[i] = b[i];
        }
    }

    std::vector<double> r(n, 0.0);
    double rsold = 0.0;
    for (int i = 0; i < n; ++i) {
        double Ax = 0.0;
        for (auto const& [col, val] : A.rows[i]) Ax += val * x[col];
        r[i] = b[i] - Ax;
        rsold += r[i] * r[i];
    }

    if (std::sqrt(rsold) < 1e-12) return;

    std::vector<double> p = r;

    for (int iter = 0; iter < n * 10; ++iter) { 
        std::vector<double> Ap(n, 0.0);
        double pAp = 0.0;
        for (int i = 0; i < n; ++i) {
            for (auto const& [col, val] : A.rows[i]) {
                Ap[i] += val * p[col];
            }
            pAp += p[i] * Ap[i];
        }

        double alpha = rsold / std::max(pAp, 1e-16);

        std::vector<double> r_new(n, 0.0);
        double rsnew = 0.0;
        for (int i = 0; i < n; ++i) {
            x[i] += alpha * p[i];
            r_new[i] = r[i] - alpha * Ap[i];
            rsnew += r_new[i] * r_new[i];
        }

        if (std::sqrt(rsnew) < 1e-10) {
            std::cout << "    -> Solver Converged in " << iter << " iterations.\n";
            break;
        }

        for (int i = 0; i < n; ++i) {
            p[i] = r_new[i] + (rsnew / rsold) * p[i];
        }
        r = r_new;
        rsold = rsnew;
    }
}

// 3. Piecewise IBM Assembly Function
void assemble_matrices_ibm(const Mesh &msh, Flow &flow, 
                           const std::vector<double>& element_k, 
                           std::function<double(double, double)> dynamic_source)
{
    int n = msh.n_nodes;
    flow.Ktau.assign(n, std::vector<double>(n, 0.0));
    flow.fg.assign(n, std::vector<double>(2, 0.0));
    SparseRowMap K_total(n);

    gaussQuad(2); 

    for (int e=0; e<msh.n_elements; ++e) {
        double k_elem = element_k[e]; 
        auto conn = msh.Elements[e];
        int nElemNodes = 4;

        std::vector<std::pair<double,double>> xy(nElemNodes);
        for (int a=0; a<nElemNodes; ++a) xy[a] = msh.Coords[conn[a]];

        std::vector<std::vector<double>> Ke_tau(nElemNodes, std::vector<double>(nElemNodes, 0.0));
        std::vector<std::vector<double>> fe_g(nElemNodes, std::vector<double>(2, 0.0));

        for (size_t gp=0; gp<gaussPts.size(); ++gp) {
            double xi = gaussPts[gp][0];
            double eta = gaussPts[gp][1];

            std::vector<double> N(nElemNodes);
            std::vector<std::array<double,2>> dN(nElemNodes);
            shapeFunction(xi, eta, 1, N, dN);

            double J[2][2] = {{0}, {0}};
            double x_phys = 0, y_phys = 0;
            
            for (int a=0; a<nElemNodes; ++a) {
                J[0][0] += dN[a][0] * xy[a].first;  J[0][1] += dN[a][0] * xy[a].second; 
                J[1][0] += dN[a][1] * xy[a].first;  J[1][1] += dN[a][1] * xy[a].second; 
                x_phys += N[a] * xy[a].first;       y_phys += N[a] * xy[a].second;
            }

            double detJ = J[0][0]*J[1][1] - J[0][1]*J[1][0];
            double invJ[2][2] = {{ J[1][1]/detJ, -J[0][1]/detJ}, {-J[1][0]/detJ,  J[0][0]/detJ}};

            std::vector<std::array<double, 2>> dNdx(nElemNodes);
            for (int a=0; a<nElemNodes; ++a) {
                dNdx[a][0] = invJ[0][0]*dN[a][0] + invJ[0][1]*dN[a][1]; 
                dNdx[a][1] = invJ[1][0]*dN[a][0] + invJ[1][1]*dN[a][1]; 
            }

            double wDetJ = gaussWts[gp] * detJ;
            double src = dynamic_source(x_phys, y_phys); 

            for (int i=0; i<nElemNodes; ++i) {
                for (int j=0; j<nElemNodes; ++j) {
                    double val = k_elem * (dNdx[i][0]*dNdx[j][0] + dNdx[i][1]*dNdx[j][1]);
                    Ke_tau[i][j] += val * wDetJ;
                }
                fe_g[i][0] += N[i] * src * wDetJ;
            }
        }

        for (int i=0; i<nElemNodes; ++i) {
            int gi = conn[i];
            for (int j=0; j<nElemNodes; ++j) {
                int gj = conn[j];
                double val = Ke_tau[i][j];
                if (std::abs(val) > 1e-16) K_total.add(gi, gj, val);
            }
            flow.fg[gi][0] += fe_g[i][0];
        }
    } 
    
    // --- APPLY DIRICHLET BCs ---
    std::vector<double> F_temp(n, 0.0);
    for(int i=0; i<n; ++i) F_temp[i] = flow.fg[i][0]; 
    apply_dirichlet(msh, K_total, F_temp, flow.DirichletNodes, flow.DirichletValues);

    // --- EXECUTE INTERNAL LINEAR SOLVER ---
    std::vector<double> T_solution;
    solve_CG(K_total, F_temp, T_solution);

    // --- STORE SECURELY IN FLOW STRUCT ---
    flow.Ui.assign(n, std::vector<double>(2, 0.0));
    for(int i=0; i<n; ++i) {
        flow.Ui[i][0] = T_solution[i]; 
    }
}

// 4. Apply Dirichlet BCs
void apply_dirichlet(const Mesh &msh, SparseRowMap &K, std::vector<double> &F,
                     std::vector<std::vector<int>> &DirichletNodes,
                     std::vector<std::vector<double>> &DirichletValues)
{
    std::vector<int> bc_dofs;
    std::vector<double> bc_vals;
    
    int dir = 0;
    if (dir < DirichletNodes.size() && dir < DirichletValues.size()) {
        for (size_t j = 0; j < DirichletNodes[dir].size(); ++j) {
            bc_dofs.push_back(DirichletNodes[dir][j]);
            bc_vals.push_back(DirichletValues[dir][j]);
        }
    }

    int n_bc = bc_dofs.size();
    for (int k = 0; k < n_bc; ++k) {
        int dof = bc_dofs[k];
        double val = bc_vals[k];
        for (int i = 0; i < K.n; ++i) {
            auto it = K.rows[i].find(dof); 
            if (it != K.rows[i].end()) {
                F[i] -= it->second * val; 
                K.rows[i].erase(it); 
            }
        }
    }
    for (int k = 0; k < n_bc; ++k) {
        int dof = bc_dofs[k];
        F[dof] = bc_vals[k];             
        K.rows[dof].clear();      
        K.rows[dof][dof] = 1.0;   
    }
}