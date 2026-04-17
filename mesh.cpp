// =======================================
// File: mesh.cpp (Q2 / 9-Node Version)
// =======================================
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include "mesh.h"

void Mesh::generate() {
    double Lx = X[1] - X[0];
    double Ly = Y[1] - Y[0];

    // 1. Generate Nodes (2*nx + 1 and 2*ny + 1 for Q2 elements)
    int n_nodes_x = 2 * nx + 1;
    int n_nodes_y = 2 * ny + 1;
    n_nodes = n_nodes_x * n_nodes_y;
    Coords.resize(n_nodes);

    for(int j=0; j < n_nodes_y; ++j){
        for(int i=0; i < n_nodes_x; ++i){
            int k = j * n_nodes_x + i;
            double xk = X[0] + Lx * double(i) / (2.0 * nx);
            double yk = Y[0] + Ly * double(j) / (2.0 * ny);
            Coords[k] = {xk, yk};
        }
    }

    // 2. Generate Q2 Elements (9 nodes per element)
    n_elements = nx * ny;
    Elements.resize(n_elements);

    for(int j=0; j < ny; ++j){
        for(int i=0; i < nx; ++i){
            // Standard VTK_BIQUADRATIC_QUAD (Type 28) Node Ordering:
            // Corners (0-3), Mid-sides (4-7), Center (8)
            int n0 = (2 * j)     * n_nodes_x + (2 * i);         // Bottom-Left
            int n1 = (2 * j)     * n_nodes_x + (2 * i + 2);     // Bottom-Right
            int n2 = (2 * j + 2) * n_nodes_x + (2 * i + 2);     // Top-Right
            int n3 = (2 * j + 2) * n_nodes_x + (2 * i);         // Top-Left
            
            int n4 = (2 * j)     * n_nodes_x + (2 * i + 1);     // Bottom-Mid
            int n5 = (2 * j + 1) * n_nodes_x + (2 * i + 2);     // Right-Mid
            int n6 = (2 * j + 2) * n_nodes_x + (2 * i + 1);     // Top-Mid
            int n7 = (2 * j + 1) * n_nodes_x + (2 * i);         // Left-Mid
            
            int n8 = (2 * j + 1) * n_nodes_x + (2 * i + 1);     // Center

            Elements[j * nx + i] = {n0, n1, n2, n3, n4, n5, n6, n7, n8};
        }
    }

    // 3. Define Edge Connectivity (Now 3 nodes per edge segment)
    int n_edges = 4;
    EdgeNodeConnectivity.resize(n_edges);

    // South (j=0)
    EdgeNodeConnectivity[0].resize(nx);
    for (int i=0; i<nx; ++i) EdgeNodeConnectivity[0][i] = {2*i, 2*i+1, 2*i+2};

    // East (i=2*nx)
    EdgeNodeConnectivity[1].resize(ny);
    for (int j=0; j<ny; ++j) {
        int base = 2*nx;
        EdgeNodeConnectivity[1][j] = {(2*j)*n_nodes_x + base, (2*j+1)*n_nodes_x + base, (2*j+2)*n_nodes_x + base};
    }

    // North (j=2*ny)
    EdgeNodeConnectivity[2].resize(nx);
    for (int i=0; i<nx; ++i) {
        int base = (2*ny) * n_nodes_x;
        EdgeNodeConnectivity[2][i] = {base + 2*nx - 2*i, base + 2*nx - (2*i+1), base + 2*nx - (2*i+2)};
    }

    // West (i=0)
    EdgeNodeConnectivity[3].resize(ny);
    for (int j=0; j<ny; ++j) {
        EdgeNodeConnectivity[3][j] = {(2*ny - 2*j)*n_nodes_x, (2*ny - (2*j+1))*n_nodes_x, (2*ny - (2*j+2))*n_nodes_x};
    }
}

void Mesh::checkNormals() {
    for (int iElem=0; iElem<n_elements; ++iElem) {
        const auto& el = Elements[iElem];
        // Checking corners to ensure counter-clockwise orientation
        double x1 = Coords[el[0]].first, y1 = Coords[el[0]].second;
        double x2 = Coords[el[1]].first, y2 = Coords[el[1]].second;
        double x3 = Coords[el[2]].first, y3 = Coords[el[2]].second;
        if (((x2 - x1)*(y3 - y1) - (x3 - x1)*(y2 - y1)) <= 0) {
            std::cerr << "Warning: Element " << iElem << " might be inverted.\n";
        }
    }
    std::cout << "Mesh Generation OK (Q2 - 9 Node).\n";
}