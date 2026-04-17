// =======================================
// File: mesh.cpp
// =======================================

#include <iostream>
#include <vector>
#include <algorithm> // For std::copy
#include <iterator>  // For std::begin and std::end
#include <fstream>
#include <unordered_map>

#include "mesh.h"
void Mesh::generate() {

	// Calculating size of rectangular domain
	double Lx = X[1] - X[0];
	double Ly = X[1] - X[0];

	// Initializing List of Coordinates
	double xk, yk;

	n_nodes = (nx*Ordr+1)*(ny*Ordr+1);
	Coords.resize(n_nodes);
	for(int j=0; j<=ny*Ordr; ++j){
		for(int i=0; i<=nx*Ordr; ++i){
			int k = j*(nx*Ordr+1)+i;
			xk = X[0] + Lx*double(i)/(nx*Ordr);
			yk = Y[0] + Ly*double(j)/(ny*Ordr);
			Coords[k] = {xk, yk};
		}
	}


	// Initializing Element Connectivity Table
	n_elements = nx*ny*2;
	Elements.resize(n_elements, std::vector<int>(3*Ordr));
	int e = 0;
	for(int j=0; j<ny; j++){
		for(int i=0; i<nx; i++){
			if(Ordr>=1){
				//Elements[e].resize(3*Ordr);
				// lower left triangle
				int l0 = (j*Ordr)*(nx*Ordr+1) + i;
				int l1 = l0 + 1*Ordr;
				int l2 = l0 + Ordr*(nx*Ordr+1);
				if(Ordr==1) Elements[e++] = {l0, l1, l2};
				if(Ordr==2){
					int l3 = l0 + 1;
					int l4 = l3 + (nx*Ordr+1);
					int l5 = l4 - 1;
					Elements[e++] = {l0, l1, l2, l3, l4, l5};
				}

				// upper right triangle
                                int u0 = l2 + 1*Ordr;
                                int u1 = l2;
                                int u2 = l1;
				if(Ordr==1) Elements[e++] = {u0, u1, u2};
                                if(Ordr==2){
                                        int u3 = u0 - 1;
                                        int u4 = u3 - (nx*Ordr+1);
                                        int u5 = u4 + 1;
					Elements[e++] = {u0, u1, u2, u3, u4, u5};
                                }
                        }
			
		}
	}


	// Defining Nodal Connectivity for Edges
	int n_edges = 4;	// number of edges in the rectangular domain
	EdgeNodeConnectivity.resize(n_edges);
	for (int i=0; i<n_edges; ++i) {
		if (i==0) { // South Edge
			EdgeNodeConnectivity[i].resize(nx);
			for (int j=0; j<nx; ++j) {
				EdgeNodeConnectivity[i][j].resize(Ordr+1);
				for (int k=0; k<=Ordr; ++k) {
					EdgeNodeConnectivity[i][j][k] = j*Ordr+k;
				}
			}
		}
		if (i==1) { // East Edge
                        EdgeNodeConnectivity[i].resize(ny);
                        for (int j=0; j<ny; ++j) {
                                EdgeNodeConnectivity[i][j].resize(Ordr+1);
                                for (int k=0; k<=Ordr; ++k) {
                                        EdgeNodeConnectivity[i][j][k] = (j*Ordr+1+k)*(nx*Ordr+1)-1;
                                }
                        }
                }
		if (i==2) { // North Edge
                        EdgeNodeConnectivity[i].resize(nx);
                        for (int j=0; j<nx; ++j) {
                                EdgeNodeConnectivity[i][j].resize(Ordr+1);
                                for (int k=0; k<=Ordr; ++k) {
                                        EdgeNodeConnectivity[i][j][k] = (ny*Ordr+1)*(nx*Ordr+1)-(j*Ordr+1)-k;
                                }
                        }
                }
		if (i==3) { // West Edge
                        EdgeNodeConnectivity[i].resize(ny);
                        for (int j=0; j<ny; ++j) {
                                EdgeNodeConnectivity[i][j].resize(Ordr+1);
                                for (int k=0; k<=Ordr; ++k) {
                                        EdgeNodeConnectivity[i][j][k] = ((ny-j)*Ordr-k)*(nx*Ordr+1);
                                }
                        }
                }
	}


	/*		
	// display nodal connectivity of edges
	 for (int i=0; i<n_edges; ++i) {
		std::cout << i << "\n";
        	for (int j=0; j<nx; ++j) {
                	for (int k=0; k<=Ordr; ++k) {
                                std::cout << EdgeNodeConnectivity[i][j][k] << "\t";
                        }
                        std::cout << "\n";
                }
      	}
	*/
	
	
	
}



void Mesh::checkNormals() {

	for (int iElem=0; iElem<n_elements; ++iElem) {
		double area = 0.0;
		for (int ind=0; ind<3*Ordr; ++ind) {
			int i1 = ind;
			int i2 = (ind+1)%(3*Ordr);
			
			area += 0.5*( Coords[Elements[iElem][i1]].first*Coords[Elements[iElem][i2]].second
				    - Coords[Elements[iElem][i2]].first*Coords[Elements[iElem][i1]].second );
		}
		if (area<0) {
			int tp = Elements[iElem][2];
			Elements[iElem][2] = Elements[iElem][1];
			Elements[iElem][1] = tp;
		}
		else if (area == 0) {
			throw std::runtime_error("Mesh Error: Zero area element found");
		}
	}
	std::cout << "Mesh OK \n";
}
