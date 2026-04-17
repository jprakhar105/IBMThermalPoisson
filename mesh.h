// ========================================
// File: mesh.h
// ========================================

#ifndef MESH_H
#define MESH_H
#include <vector>
struct Mesh {
	int nx, ny, Ordr;
	double X[2], Y[2];

	int n_nodes, n_elements, n_edges;
	std::vector<std::pair<double,double>> Coords;
	std::vector<std::vector<int>> Elements;
	std::vector<std::vector<std::vector<int>>> EdgeNodeConnectivity;
	void generate();
	void checkNormals();
};
#endif
