// ======================================================
// File: vtkio.cpp
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

    // Write node coordinates
    ofs << "POINTS " << msh.n_nodes << " double\n";
    ofs << std::scientific << std::setprecision(6);
    for (int i = 0; i < msh.n_nodes; ++i)
        ofs << msh.Coords[i].first << " "
            << msh.Coords[i].second << " 0.0\n";

    // Triangular cells: 3 nodes per cell
    int nc = msh.n_elements;
    ofs << "\nCELLS " << nc << " " << nc * 4 << "\n";
    for (int e = 0; e < nc; ++e) {
        auto c = msh.Elements[e];
        ofs << 3 << " " << c[0] << " " << c[1]
                     << " " << c[2] << "\n";
    }

    // VTK cell type for triangle = 5
    ofs << "\nCELL_TYPES " << nc << "\n";
    for (int e = 0; e < nc; ++e) ofs << 5 << "\n";

    // Write vector field at nodes
    ofs << "\nPOINT_DATA " << msh.n_nodes << "\n";
    ofs << "VECTORS u double\n";

    // VTK expects 3 components -> (ux, uy, 0.0)
    for (int i = 0; i < msh.n_nodes; ++i)
        ofs << U[i][0] << " " << U[i][1] << " 0.0\n";
}

