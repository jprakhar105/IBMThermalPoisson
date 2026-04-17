// ========================================
// File: main.cpp (Windows Compatible)
// ========================================

#include <iostream>
#include <fstream>
#include <cmath>
#include <chrono>
#include <vector>
#include <iomanip>

// Windows-specific headers for memory tracking
#include <windows.h>
#include <psapi.h>

#include "mesh.h"
#include "flowParameters.h"
#include "matrixAssembly.h"
#include "writeOutput.h"

// Helper to get peak memory usage in KB (Windows version)
long getPeakMemory() {
    PROCESS_MEMORY_COUNTERS pmc;
    if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc))) {
        // PeakWorkingSetSize returns bytes, convert to KB
        return (long)(pmc.PeakWorkingSetSize / 1024);
    }
    return 0;
}

int main(){
    std::cout << "Starting Simulation Loop (N=5 to 100)..." << "\n";

    // Output file for stats
    std::ofstream statsFile("simulation_stats.csv");
    statsFile << "N,Nodes,L1_Error,L2_Error,Linf_Error,Time_ms,Memory_KB\n";

    // Loop through Mesh resolutions
    for (int n_elem = 5; n_elem <= 100; n_elem++) {
        
        // Start Timer
        auto start = std::chrono::high_resolution_clock::now();

        // ----------------------------------------
        // 1. Setup Simulation
        // ----------------------------------------
        Flow flow;
        flow.setParameters();

        Mesh msh;
        msh.X[0] = -1; msh.X[1] = 1;
        msh.Y[0] = -1; msh.Y[1] = 1;
        msh.nx = n_elem; 
        msh.ny = n_elem;
        msh.n_edges = 4;
        msh.Ordr = 1; // Linear elements

        // ----------------------------------------
        // 2. Generate Mesh & Matrices
        // ----------------------------------------
        msh.generate();
        msh.checkNormals();

        // Dirichlet BCs (All 4 edges)
        flow.DirichletEdges.resize(2);
        flow.DirichletEdges[0] = {0,1,2,3}; 
        flow.DirichletEdges[1] = {0,1,2,3}; 
        flow.calculateDirichletNodes(msh);
        flow.calculateDirichletValues(msh);

        assemble_matrices(msh, flow);

        // ----------------------------------------
        // 3. Solve System
        // ----------------------------------------
        int N = msh.n_nodes;
        flow.Ui.resize(N, std::vector<double>(2));
        flow.initializeFlow(msh, 1); // 1 = Calculate Steady State

        // ----------------------------------------
        // 4. Analytical Solution & Error
        // ----------------------------------------
        flow.UA.resize(N, std::vector<double>(2));
        flow.calculateAnalyticalSolution(msh);

        double L1 = 0.0, L2 = 0.0, Linf = 0.0;
        
        for (int i=0; i<N; ++i) {
            // Using direction 0 (scalar field T)
            double err = std::abs(flow.Ui[i][0] - flow.UA[i][0]);
            
            L1 += err;
            L2 += err * err;
            if (err > Linf) Linf = err;
        }
        L1 /= N;
        L2 = std::sqrt(L2 / N);

        // Stop Timer
        auto end = std::chrono::high_resolution_clock::now();
        double elapsed = std::chrono::duration<double, std::milli>(end - start).count();
        long memory = getPeakMemory();

        // ----------------------------------------
        // 5. Report & Save
        // ----------------------------------------
        std::cout << "N=" << n_elem 
                  << " | Nodes=" << N 
                  << " | L_inf=" << Linf 
                  << " | Time=" << elapsed << "ms\n";

        statsFile << n_elem << "," << N << "," 
                  << L1 << "," << L2 << "," << Linf << "," 
                  << elapsed << "," << memory << "\n";

        // Save the detailed field only for N=100 (optional)
        if (n_elem == 100) {
            std::vector<std::vector<double>> UErr(N, std::vector<double>(2));
            for(int i=0; i<N; ++i) UErr[i][0] = flow.Ui[i][0] - flow.UA[i][0];
            
            char fname[256];
            std::sprintf(fname, "solution_N%03d.dat", n_elem);
            write_output(fname, msh, flow.Ui, flow.UA, UErr);
        }
    }

    statsFile.close();
    std::cout << "Simulation Complete. Stats saved to simulation_stats.csv\n";
    return 0;
}