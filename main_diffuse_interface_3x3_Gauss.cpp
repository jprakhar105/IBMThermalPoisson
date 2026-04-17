// ========================================
// File: main.cpp 
// Immersed Boundary Method Benchmark (Q2 / Diffuse Interface)
// ========================================

#include <iostream>
#include <fstream>
#include <cmath>
#include <chrono>
#include <vector>
#include <iomanip>
#include <functional>
#include <cstdlib>
#include <string>

#include "mesh.h"
#include "flowParameters.h"
#include "matrixAssembly.h"
#include "writeOutput.h"

double exact_T_IBM(double x, double y) {
    double R = 0.5;
    double k1 = 1.0;
    double k2 = 10.0;
    
    double r = std::sqrt(x*x + y*y);
    double a = (k1 + k2) / (2.0 * k2);           
    double b = (R*R * (k2 - k1)) / (2.0 * k2);   

    if (r < R) return x; 
    else       return a * x + (b * x) / (r * r); 
}

void plot_ibm_results() {
    std::cout << "Generating 2D contour comparison...\n";
    std::ofstream pyFile("plot_ibm.py");
    
    pyFile << "import csv\nimport numpy as np\nimport matplotlib.pyplot as plt\nimport matplotlib.tri as tri\n\n";
    pyFile << "def load_data(filename):\n    x, y, t, exact, err = [], [], [], [], []\n    with open(filename, 'r') as f:\n        reader = csv.reader(f)\n        next(reader)\n        for row in reader:\n            x.append(float(row[0]))\n            y.append(float(row[1]))\n            t.append(float(row[2]))\n            exact.append(float(row[3]))\n            err.append(float(row[4]))\n    return x, y, t, exact, err\n\n";
    
    pyFile << "x1, y1, t1, ex1, err1 = load_data('results_ibm.csv')\n";
    pyFile << "triang = tri.Triangulation(x1, y1)\n";
    
    pyFile << "fig, axs = plt.subplots(1, 3, figsize=(18, 5))\n\n";
    
    pyFile << "tpc1 = axs[0].tricontourf(triang, t1, levels=50, cmap='jet')\nfig.colorbar(tpc1, ax=axs[0])\naxs[0].set_title('FEM Temperature')\naxs[0].set_aspect('equal')\n";
    pyFile << "tpc2 = axs[1].tricontourf(triang, ex1, levels=50, cmap='jet')\nfig.colorbar(tpc2, ax=axs[1])\naxs[1].set_title('Exact Analytical Solution')\naxs[1].set_aspect('equal')\n";
    
    pyFile << "ea1 = np.array(err1)\nif np.ptp(ea1) < 1e-12: ea1[0] += 1e-12\n";
    pyFile << "epc1 = axs[2].tricontourf(triang, ea1, levels=50, cmap='inferno')\nfig.colorbar(epc1, ax=axs[2])\naxs[2].set_title('Absolute Error |T_fem - T_exact|')\naxs[2].set_aspect('equal')\n\n";
    
    pyFile << "plt.tight_layout()\nplt.savefig('ibm_comparison.png', dpi=300)\nplt.show()\n";
    pyFile.close();

    std::system("python plot_ibm.py");
}

int main(){
    std::cout << "Starting IBM Benchmark Solver (Order 2, Nodal Phase)...\n";

    for (int n_elem = 20; n_elem <= 100; n_elem += 80) { 
        Flow flow;
        flow.setParameters();
        Mesh msh;
        msh.X[0] = -1; msh.X[1] = 1; msh.Y[0] = -1; msh.Y[1] = 1;
        msh.nx = n_elem; msh.ny = n_elem; msh.n_edges = 4; 
        msh.Ordr = 2; 
        msh.generate(); 
        msh.checkNormals();

        // ---------------------------------------------------------
        // IBM Step 1: Nodal Phase Fraction Assignment (Diffuse Interface)
        // ---------------------------------------------------------
        std::vector<double> node_phi(msh.n_nodes, 0.0);
        double R = 0.5; 
        double k1 = 1.0;
        double k2 = 10.0;

        for (int i = 0; i < msh.n_nodes; ++i) {
            double x = msh.Coords[i].first;
            double y = msh.Coords[i].second;
            double r = std::sqrt(x*x + y*y);
            
            if (r < R) node_phi[i] = 1.0;  // Inside phase
            else       node_phi[i] = 0.0;  // Outside phase
        }

        auto dynamic_source = [](double x, double y) { return 0.0; };

        flow.DirichletEdges.resize(2); 
        flow.DirichletEdges[0] = {0, 1, 2, 3}; 
        flow.DirichletEdges[1] = {};           
        flow.calculateDirichletNodes(msh);
        
        flow.DirichletValues[0].resize(flow.DirichletNodes[0].size());
        for(size_t i = 0; i < flow.DirichletNodes[0].size(); ++i) {
            int node_id = flow.DirichletNodes[0][i];
            double x = msh.Coords[node_id].first;
            double y = msh.Coords[node_id].second;
            flow.DirichletValues[0][i] = exact_T_IBM(x, y); 
        }

        assemble_matrices_ibm(msh, flow, node_phi, k1, k2, dynamic_source);

        int N = msh.n_nodes;
        flow.UA.resize(N, std::vector<double>(2));
        double numerator_sum = 0.0;
        double denominator_sum = 0.0;
        
        for (int i=0; i<N; ++i) {
            flow.UA[i][0] = exact_T_IBM(msh.Coords[i].first, msh.Coords[i].second);
            flow.UA[i][1] = 0.0; 
            
            double diff = flow.Ui[i][0] - flow.UA[i][0];
            numerator_sum += (diff * diff);
            denominator_sum += (flow.UA[i][0] * flow.UA[i][0]);
        }
        
        double L2_error = std::sqrt(numerator_sum / denominator_sum);

        if (n_elem == 100) {
            std::cout << "Mesh N=" << n_elem << " | Relative L2 Error = " << L2_error << "\n";
            
            std::ofstream plotFile("results_ibm.csv");
            plotFile << "X,Y,T_FEM,T_Exact,Error\n";
            for (int i = 0; i < N; ++i) {
                plotFile << msh.Coords[i].first << "," << msh.Coords[i].second << "," 
                         << flow.Ui[i][0] << "," << flow.UA[i][0] << "," 
                         << std::abs(flow.Ui[i][0] - flow.UA[i][0]) << "\n";
            }
            plotFile.close();
        }
    }
    
    plot_ibm_results();
    
    std::cout << "IBM Benchmark Testing Complete.\n";
    return 0;
}