// ========================================
// File: main.cpp (Q4 Loop Version)
// ========================================

#include <iostream>
#include <fstream>
#include <cmath>
#include <chrono>
#include <vector>
#include <iomanip>

// Use Windows header if on Windows, otherwise sys/resource.h
#ifdef _WIN32
#include <windows.h>
#include <psapi.h>
#else
#include <sys/resource.h>
#endif

#include "mesh.h"
#include "flowParameters.h"
#include "matrixAssembly.h"
#include "writeOutput.h"

// Helper for memory usage
long getPeakMemory() {
#ifdef _WIN32
    PROCESS_MEMORY_COUNTERS pmc;
    if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc))) {
        return (long)(pmc.PeakWorkingSetSize / 1024);
    }
    return 0;
#else
    struct rusage r_usage;
    getrusage(RUSAGE_SELF, &r_usage);
    return r_usage.ru_maxrss; 
#endif
}

int main(){
    std::cout << "Starting Q4 Simulation Loop (N=5 to 100)..." << "\n";

    std::ofstream statsFile("simulation_stats_Q4.csv");
    statsFile << "N,Nodes,L1_Error,L2_Error,Linf_Error,Time_ms,Memory_KB\n";

    for (int n_elem = 5; n_elem <= 100; n_elem++) {
        
        auto start = std::chrono::high_resolution_clock::now();

        // 1. Setup Simulation
        Flow flow;
        flow.setParameters();

        Mesh msh;
        msh.X[0] = -1; msh.X[1] = 1;
        msh.Y[0] = -1; msh.Y[1] = 1;
        msh.nx = n_elem; 
        msh.ny = n_elem;
        msh.n_edges = 4;
        msh.Ordr = 1; // Linear Quads

        // 2. Generate Mesh (Q4) & Matrices
        msh.generate();
        msh.checkNormals();

        // 3. Define BCs
        flow.DirichletEdges.resize(2);
        flow.DirichletEdges[0] = {0,1,2,3}; 
        flow.DirichletEdges[1] = {0,1,2,3}; 
        flow.calculateDirichletNodes(msh);
        flow.calculateDirichletValues(msh);

        assemble_matrices(msh, flow);

        // 4. Solve System
        int N = msh.n_nodes;
        flow.Ui.resize(N, std::vector<double>(2));
        flow.initializeFlow(msh, 1);

        // 5. Analytical Solution & Error
        flow.UA.resize(N, std::vector<double>(2));
        flow.calculateAnalyticalSolution(msh);

        double L1 = 0.0, L2 = 0.0, Linf = 0.0;
        for (int i=0; i<N; ++i) {
            double err = std::abs(flow.Ui[i][0] - flow.UA[i][0]);
            L1 += err;
            L2 += err * err;
            if (err > Linf) Linf = err;
        }
        L1 /= N;
        L2 = std::sqrt(L2 / N);

        auto end = std::chrono::high_resolution_clock::now();
        double elapsed = std::chrono::duration<double, std::milli>(end - start).count();
        long memory = getPeakMemory();

        std::cout << "N=" << n_elem << " | Nodes=" << N << " | L_inf=" << Linf << "\n";
        statsFile << n_elem << "," << N << "," << L1 << "," << L2 << "," << Linf << "," << elapsed << "," << memory << "\n";

        // Output N=100 specific file
        if (n_elem == 100) {
            std::vector<std::vector<double>> UErr(N, std::vector<double>(2));
            for(int i=0; i<N; ++i) UErr[i][0] = flow.Ui[i][0] - flow.UA[i][0];
            
            char fname[256];
            std::sprintf(fname, "solution_Q4_N%03d.dat", n_elem);
            write_output(fname, msh, flow.Ui, flow.UA, UErr);
        }
    }

    statsFile.close();
    std::cout << "Simulation Complete.\n";
    return 0;
}