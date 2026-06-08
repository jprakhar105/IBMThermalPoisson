// =======================================
// File: flowParameters.h
// =======================================

#ifndef FLOWPARAMETERS_H
#define FLOWPARAMETERS_H
#include <vector>
#include <iostream>
#include <functional>

#include "mesh.h"
#include "sparse.h"

struct Flow {

    // Global constants
    double pi;

    // Fluid properties
    double Re;

    double g(double x, double y, int dir);  // gravity

    // Time settings
    double dt;
    double nSteps;
    double theta;

    void setParameters();

    // Global matrices
    CSR K_tau;

    // Global Vectors
    std::vector<std::vector<double>> fg;

    // Initial Velocity/Temperature Vector
    std::vector<std::vector<double>> Ui;

    void initializeFlow(const Mesh &msh, int init);

    // Dirichlet BC vectors
    std::vector<std::vector<int>>    DirichletEdges;
    std::vector<std::vector<int>>    DirichletNodes;
    std::vector<std::vector<double>> DirichletValues;
    std::vector<int> nDirichletNodes;

    void calculateDirichletNodes(const Mesh &msh);
    void calculateDirichletValues(const Mesh &msh);

    // Analytical Velocity/Temperature Vector
    std::vector<std::vector<double>> UA;
    void calculateAnalyticalSolution(const Mesh &msh);
};

#endif