// =======================================
// File: mesh.cpp (Q4 Version)
// =======================================
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include "mesh.h"

void Mesh::generate() {
    double Lx = X[1] - X[0];
    double Ly = Y[1] - Y[0];

    // 1. Generate Nodes
    n_nodes = (nx + 1) * (ny + 1);
    Coords.resize(n_nodes);

    double xk, yk;
    for(int j=0; j<=ny; ++j){
        for(int i=0; i<=nx; ++i){
            int k = j*(nx + 1) + i;
            xk = X[0] + Lx * double(i) / nx;
            yk = Y[0] + Ly * double(j) / ny;
            Coords[k] = {xk, yk};
        }
    }

    // 2. Generate Q4 Elements (1 per cell)
    n_elements = nx * ny;
    Elements.resize(n_elements);

    int e = 0;
    int row_len = nx + 1;

    for(int j=0; j<ny; ++j){
        for(int i=0; i<nx; ++i){
            // Nodes: Counter-Clockwise from Bottom-Left
            int n1 = j * row_len + i;           
            int n2 = n1 + 1;                    
            int n3 = n1 + row_len + 1;          
            int n4 = n1 + row_len;              

            Elements[e++] = {n1, n2, n3, n4};
        }
    }

    // 3. Define Edge Connectivity
    int n_edges = 4;
    EdgeNodeConnectivity.resize(n_edges);

    // South (j=0)
    EdgeNodeConnectivity[0].resize(nx);
    for (int i=0; i<nx; ++i) EdgeNodeConnectivity[0][i] = {i, i+1};

    // East (i=nx)
    EdgeNodeConnectivity[1].resize(ny);
    for (int j=0; j<ny; ++j) {
        int n1 = j*row_len + nx;
        int n2 = (j+1)*row_len + nx;
        EdgeNodeConnectivity[1][j] = {n1, n2};
    }

    // North (j=ny)
    EdgeNodeConnectivity[2].resize(nx);
    for (int i=0; i<nx; ++i) {
        int n1 = ny*row_len + (nx-i);
        int n2 = ny*row_len + (nx-i-1);
        EdgeNodeConnectivity[2][i] = {n1, n2};
    }

    // West (i=0)
    EdgeNodeConnectivity[3].resize(ny);
    for (int j=0; j<ny; ++j) {
        int n1 = (ny-j)*row_len;
        int n2 = (ny-j-1)*row_len;
        EdgeNodeConnectivity[3][j] = {n1, n2};
    }
}

void Mesh::checkNormals() {
    for (int iElem=0; iElem<n_elements; ++iElem) {
        const auto& el = Elements[iElem];
        double x1 = Coords[el[0]].first, y1 = Coords[el[0]].second;
        double x2 = Coords[el[1]].first, y2 = Coords[el[1]].second;
        double x3 = Coords[el[2]].first, y3 = Coords[el[2]].second;
        if (((x2 - x1)*(y3 - y1) - (x3 - x1)*(y2 - y1)) <= 0) {
            std::cerr << "Warning: Element " << iElem << " might be inverted.\n";
        }
    }
    std::cout << "Mesh Generation OK (Q4).\n";
}