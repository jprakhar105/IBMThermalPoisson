// ======================================================
// File: vtkio.cpp (Q4 Update)
// ======================================================
#include "vtkio.h"
#include <fstream>
#include <iomanip>

void write_vtk(const std::string &fname,
               const Mesh &msh,
               const std::vector<std::vector<double>> &U)
{
    std::ofstream ofs(fname);

    ofs << "# vtk DataFile Version 2.0\n";
    ofs << "Vector field (ux, uy)\nASCII\n";
    ofs << "DATASET UNSTRUCTURED_GRID\n";

    ofs << "POINTS " << msh.n_nodes << " double\n";
    ofs << std::scientific << std::setprecision(6);
    for (int i = 0; i < msh.n_nodes; ++i)
        ofs << msh.Coords[i].first << " " << msh.Coords[i].second << " 0.0\n";

    int nc = msh.n_elements;
    // Format: [count=4] [n1] [n2] [n3] [n4] -> total size 5
    ofs << "\nCELLS " << nc << " " << nc * 5 << "\n";
    for (int e = 0; e < nc; ++e) {
        auto c = msh.Elements[e];
        ofs << 4 << " " << c[0] << " " << c[1] << " " << c[2] << " " << c[3] << "\n";
    }

    // VTK_QUAD = 9
    ofs << "\nCELL_TYPES " << nc << "\n";
    for (int e = 0; e < nc; ++e) ofs << 9 << "\n";

    ofs << "\nPOINT_DATA " << msh.n_nodes << "\n";
    ofs << "VECTORS u double\n";
    for (int i = 0; i < msh.n_nodes; ++i)
        ofs << U[i][0] << " " << U[i][1] << " 0.0\n";
}