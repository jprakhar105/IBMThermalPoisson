// ======================================================
// File: writeOutput.cpp (Q4 Update)
// ======================================================
#include "writeOutput.h"
#include <fstream>
#include <iomanip>
#include <iostream>

void write_output(const std::string &fname, const Mesh &msh,
                  const std::vector<std::vector<double>> &U,
                  const std::vector<std::vector<double>> &UA,
                  const std::vector<std::vector<double>> &UErr)
{
    std::ofstream fout(fname);
    fout << "TITLE = Tecplot Data\n";
    fout << "VARIABLES = \"X\", \"Y\", \"Z\", \"Ux\", \"Uy\", \"UAx\", \"UAy\", \"UErrx\", \"UErry\"\n";
    fout << "ZONE T=\"Zone 1\", N=" << msh.n_nodes << ", E=" << msh.n_elements << ", F=FEPOINT  ET=QUADRILATERAL\n";

    for (int i=0; i<msh.n_nodes; ++i) {
        fout << msh.Coords[i].first << "\t" << msh.Coords[i].second << "\t0.0\t"
             << U[i][0] << "\t" << U[i][1] << "\t"
             << UA[i][0] << "\t" << UA[i][1] << "\t"
             << UErr[i][0] << "\t" << UErr[i][1] << "\n";
    }

    for (int i=0; i<msh.n_elements; ++i) {
        fout << msh.Elements[i][0]+1 << "\t" << msh.Elements[i][1]+1 << "\t" 
             << msh.Elements[i][2]+1 << "\t" << msh.Elements[i][3]+1 << "\n";
    }
    fout.close();
    std::cout << "Tecplot file exported: " << fname << std::endl;
}