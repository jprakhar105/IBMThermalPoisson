// ======================================================
// File: vtkio.cpp (Q2 / 9-Node Update)
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
    // Format: [count=9] [n1] [n2] ... [n9] -> total size 10 per cell
    ofs << "\nCELLS " << nc << " " << nc * 10 << "\n";
    for (int e = 0; e < nc; ++e) {
        auto c = msh.Elements[e];
        ofs << 9 << " " << c[0] << " " << c[1] << " " << c[2] << " " 
                 << c[3] << " " << c[4] << " " << c[5] << " " 
                 << c[6] << " " << c[7] << " " << c[8] << "\n";
    }

    // VTK_BIQUADRATIC_QUAD = 28
    ofs << "\nCELL_TYPES " << nc << "\n";
    for (int e = 0; e < nc; ++e) ofs << 28 << "\n";

    ofs << "\nPOINT_DATA " << msh.n_nodes << "\n";
    ofs << "VECTORS u double\n";
    for (int i = 0; i < msh.n_nodes; ++i)
        ofs << U[i][0] << " " << U[i][1] << " 0.0\n";
}