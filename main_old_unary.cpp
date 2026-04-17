// ========================================
// File: main.cpp (Mixed BC + 2D Plotting)
// The "Opposite Constant" Flux Test
// Case 1: 4 Dirichlet
// Case 2: 3 Dirichlet + 1 Constant Neumann (q = -1.0)
// ========================================

#include <iostream>
#include <fstream>
#include <cmath>
#include <chrono>
#include <vector>
#include <iomanip>
#include <functional>
#include <cstdlib>

#ifdef _WIN32
#include <windows.h>
#include <psapi.h>
#else
#include <sys/resource.h>
#endif

#include "mesh.h"
#include "flowParameters.h"
#include "matrixAssembly_unary.h"
#include "writeOutput.h"

// ---------------------------------------------------------
// 1. Exact Analytical Solution (4 Dirichlet Walls)
// ---------------------------------------------------------
const double PI = 3.14159265358979323846;

double exact_T_4Dirichlet(double x, double y) {
    double sum = 0.0;
    for (int k = 0; k < 40; ++k) {
        int n = 2 * k + 1;
        double sign = (k % 2 == 0) ? 1.0 : -1.0;
        
        double term1 = sign / (n * n * n);
        double cosh_ratio = std::cosh(n * PI * y / 2.0) / std::cosh(n * PI / 2.0);
        double term2 = 1.0 - cosh_ratio;
        double term3 = std::cos(n * PI * x / 2.0);
        
        sum += term1 * term2 * term3;
    }
    return (16.0 / (PI * PI * PI)) * sum;
}

// ---------------------------------------------------------
// 2. Utility Functions
// ---------------------------------------------------------
void apply_neumann_flux(Flow &flow, const Mesh &msh, int edgeID, std::function<double(double, double)> q_flux) {
    if (edgeID >= msh.EdgeNodeConnectivity.size()) return;
    const auto& segments = msh.EdgeNodeConnectivity[edgeID];
    const double inv_sqrt3 = 1.0 / std::sqrt(3.0);
    std::vector<double> xi = {-inv_sqrt3, inv_sqrt3};
    std::vector<double> weight = {1.0, 1.0};

    for (const auto& seg : segments) {
        int n1 = seg[0]; int n2 = seg[1];
        double x1 = msh.Coords[n1].first;  double y1 = msh.Coords[n1].second;
        double x2 = msh.Coords[n2].first;  double y2 = msh.Coords[n2].second;
        double Fe[2] = {0.0, 0.0};

        for (size_t q = 0; q < weight.size(); ++q) {
            double N1 = 0.5 * (1.0 - xi[q]); double N2 = 0.5 * (1.0 + xi[q]);
            double dN1_dxi = -0.5;           double dN2_dxi = 0.5;
            double dx_dxi = dN1_dxi * x1 + dN2_dxi * x2;
            double dy_dxi = dN1_dxi * y1 + dN2_dxi * y2;
            double detJ = std::sqrt(dx_dxi * dx_dxi + dy_dxi * dy_dxi);

            double x_gp = N1 * x1 + N2 * x2;
            double y_gp = N1 * y1 + N2 * y2;

            double q_val = q_flux(x_gp, y_gp);
            Fe[0] += N1 * q_val * detJ * weight[q];
            Fe[1] += N2 * q_val * detJ * weight[q];
        }
        flow.fg[n1][0] += Fe[0];
        flow.fg[n2][0] += Fe[1];
    }
}

void plot_results() {
    std::cout << "Generating side-by-side 2D contour comparisons...\n";
    std::ofstream pyFile("plot_comparison.py");
    
    pyFile << "import csv\nimport numpy as np\nimport matplotlib.pyplot as plt\nimport matplotlib.tri as tri\n\n";
    pyFile << "def load_data(filename):\n";
    pyFile << "    x, y, t, exact, err = [], [], [], [], []\n";
    pyFile << "    with open(filename, 'r') as f:\n";
    pyFile << "        reader = csv.reader(f)\n        next(reader)\n";
    pyFile << "        for row in reader:\n";
    pyFile << "            x.append(float(row[0]))\n            y.append(float(row[1]))\n";
    pyFile << "            t.append(float(row[2]))\n            exact.append(float(row[3]))\n            err.append(float(row[4]))\n";
    pyFile << "    return x, y, t, exact, err\n\n";
    
    pyFile << "x1, y1, t1, ex1, err1 = load_data('results_case1.csv')\n";
    pyFile << "x2, y2, t2, ex2, err2 = load_data('results_case2.csv')\n";
    pyFile << "triang = tri.Triangulation(x1, y1)\n";
    
    pyFile << "fig, axs = plt.subplots(2, 2, figsize=(14, 12))\n\n";
    
    // Case 1 Plots
    pyFile << "tpc1 = axs[0,0].tricontourf(triang, t1, levels=50, cmap='jet')\nfig.colorbar(tpc1, ax=axs[0,0])\n";
    pyFile << "axs[0,0].set_title('Case 1: 4 Dirichlet (FEM Temp)')\naxs[0,0].set_aspect('equal')\n";
    pyFile << "ea1 = np.array(err1)\nif np.ptp(ea1) < 1e-12: ea1[0] += 1e-12\n";
    pyFile << "epc1 = axs[0,1].tricontourf(triang, ea1, levels=50, cmap='inferno')\nfig.colorbar(epc1, ax=axs[0,1])\n";
    pyFile << "axs[0,1].set_title('Case 1: Error (FEM vs 4-Dirichlet Exact)')\naxs[0,1].set_aspect('equal')\n\n";

    // Case 2 Plots
    pyFile << "tpc2 = axs[1,0].tricontourf(triang, t2, levels=50, cmap='jet')\nfig.colorbar(tpc2, ax=axs[1,0])\n";
    pyFile << "axs[1,0].set_title('Case 2: 3 Dir + 1 Constant Neumann q=-1.0')\naxs[1,0].set_aspect('equal')\n";
    pyFile << "ea2 = np.array(err2)\nif np.ptp(ea2) < 1e-12: ea2[0] += 1e-12\n";
    pyFile << "epc2 = axs[1,1].tricontourf(triang, ea2, levels=50, cmap='inferno')\nfig.colorbar(epc2, ax=axs[1,1])\n";
    pyFile << "axs[1,1].set_title('Case 2: Error (FEM vs 4-Dirichlet Exact)')\naxs[1,1].set_aspect('equal')\n\n";

    pyFile << "plt.tight_layout()\nplt.savefig('solution_comparison.png', dpi=300)\nplt.show()\n";
    pyFile.close();

    std::system("python plot_comparison.py");
}

// ---------------------------------------------------------
// 3. Main Loop
// ---------------------------------------------------------
int main(){
    std::cout << "Starting The Opposite Constant Validation Test...\n";

    for (int case_id = 1; case_id <= 2; ++case_id) {
        std::cout << "\n--- Running Case " << case_id << " ---\n";
        
        for (int n_elem = 5; n_elem <= 100; n_elem += 95) { // Skip straight to N=100
            Flow flow;
            flow.setParameters();
            Mesh msh;
            msh.X[0] = -1; msh.X[1] = 1; msh.Y[0] = -1; msh.Y[1] = 1;
            msh.nx = n_elem; msh.ny = n_elem; msh.n_edges = 4; msh.Ordr = 1;
            msh.generate();
            msh.checkNormals();

            flow.DirichletEdges.resize(2);
            
            if (case_id == 1) {
                // CASE 1: True 4-Dirichlet
                flow.DirichletEdges[0] = {0, 1, 2, 3}; 
                flow.DirichletEdges[1] = {0, 1, 2, 3}; 
            } else {
                // CASE 2: 3 Dirichlet, Top is Neumann
                flow.DirichletEdges[0] = {0, 1, 3}; 
                flow.DirichletEdges[1] = {0, 1, 3}; 
            }

            flow.calculateDirichletNodes(msh);
            flow.calculateDirichletValues(msh);
            assemble_matrices(msh, flow);

            if (case_id == 2) {
                // THE OPPOSITE CONSTANT: Extracting heat at exactly q = -1.0
                apply_neumann_flux(flow, msh, 2, [](double x, double y) { return 1.0; });
            }

            int N = msh.n_nodes;
            flow.Ui.resize(N, std::vector<double>(2));
            flow.initializeFlow(msh, 1);

            // Calculate Exact Error (always using the 4-Dirichlet analytical solution)
            flow.UA.resize(N, std::vector<double>(2));
            double Linf = 0.0;
            for (int i=0; i<N; ++i) {
                flow.UA[i][0] = exact_T_4Dirichlet(msh.Coords[i].first, msh.Coords[i].second);
                flow.UA[i][1] = 0.0; 
                double err = std::abs(flow.Ui[i][0] - flow.UA[i][0]);
                if (err > Linf) Linf = err;
            }

            if (n_elem == 100) {
                std::cout << "Case " << case_id << " | N=100 | Max Error = " << Linf << "\n";
                
                std::string fname = (case_id == 1) ? "results_case1.csv" : "results_case2.csv";
                std::ofstream plotFile(fname);
                plotFile << "X,Y,T_FEM,T_Exact,Error\n";
                for (int i = 0; i < N; ++i) {
                    plotFile << msh.Coords[i].first << "," << msh.Coords[i].second << "," 
                             << flow.Ui[i][0] << "," << flow.UA[i][0] << "," 
                             << flow.Ui[i][0] - flow.UA[i][0] << "\n";
                }
                plotFile.close();
            }
        }
    }
    
    plot_results();
    std::cout << "Testing Complete.\n";
    return 0;
}