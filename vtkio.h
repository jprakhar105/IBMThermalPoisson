// ======================================================
// File: vtkio.h — Write VTK output
// ======================================================
#ifndef VTKIO_H
#define VTKIO_H
#include <string>
#include <vector>
#include "mesh.h"

void write_vtk(const std::string &fname,
               const Mesh &mesh,
               const std::vector<std::vector<double>> &U);

#endif
