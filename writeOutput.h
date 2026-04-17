// ======================================================
// File: vtkio.h — Write VTK output
// ======================================================
#ifndef WRITE_OUTPUT_H
#define WRITE_OUTPUT_H
#include <string>
#include <vector>
#include "mesh.h"

void write_output(const std::string &fname,
                  const Mesh &msh,
                  const std::vector<std::vector<double>> &U,
                  const std::vector<std::vector<double>> &UA,
                  const std::vector<std::vector<double>> &UErr );

#endif
