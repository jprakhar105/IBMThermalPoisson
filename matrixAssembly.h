// =======================================
// File: matrixAssembly.h
// =======================================
#ifndef MATRIXASSEMBLY_H
#define MATRIXASSEMBLY_H

#include "mesh.h"
#include "flowParameters.h"
#include "sparse.h"
#include <vector>
#include <functional>

// Updated signature for piecewise constant IBM
void assemble_matrices_ibm(const Mesh &msh, Flow &flow, 
                           const std::vector<double>& node_phi, 
                           double k1, double k2,
                           std::function<double(double, double)> dynamic_source);
                           
// MATLAB-style Symmetric Dirichlet Application
void apply_dirichlet(const Mesh &msh, SparseRowMap &K, std::vector<double> &F,
                     std::vector<std::vector<int>> &DirichletNodes,
                     std::vector<std::vector<double>> &DirichletValues);

#endif