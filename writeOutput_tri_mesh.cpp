// ======================================================
// File: vtkio.cpp
// ======================================================
#include "vtkio.h"
#include <fstream>
#include <iomanip>
#include <iostream>

// Write Output in TecPlot format
void write_output(const std::string &fname,
                  const Mesh &msh,
                  const std::vector<std::vector<double>> &U,
		  const std::vector<std::vector<double>> &UA,
		  const std::vector<std::vector<double>> &UErr)
{
    	// Count mesh node and element
    	size_t nNodes = msh.Coords.size();
	size_t nElems = msh.Elements.size();

	if (msh.n_nodes!=nNodes || msh.n_elements!=nElems) {
		throw std::runtime_error("Inconsistent Mesh");
	}

    	std::ofstream fout(fname);
    	if (!fout) {
        	throw std::runtime_error("Cannot open output file.");
    	}

    	// Write header
    	fout << "TITLE = Tecplot Data" << "\n";
    	fout << "VARIABLES = \"X\", \"Y\", \"Z\", \"Ux\", \"Uy\", \"UAx\", \"UAy\", \"UErrx\", \"UErry\"" << "\n";

    	// Write zone
    	fout << "ZONE T=\"Zone 1\", N=" << nNodes << ", E=" << nElems << ", F=FEPOINT  ET=TRIANGLE\n";

    	// Write data
    	for (int i=0; i<nNodes; ++i) {
        	fout << msh.Coords[i].first << "\t" << msh.Coords[i].second << "\t" << 0.0 << "\t"
             		<< U[i][0] << "\t" << U[i][1] << "\t"
		       	<< UA[i][0] << "\t" << UA[i][1] << "\t"
			<< UErr[i][0] << "\t" << UErr[i][1] << "\t"
			<< "\n";
    	}

    	// Write mesh connectivity
    	for (int i=0; i<nElems; ++i) {
		if (msh.Ordr == 1) {
			fout << msh.Elements[i][0]+1 << "\t" << msh.Elements[i][1]+1 << "\t" << msh.Elements[i][2]+1 << "\n";
		}
		else if (msh.Ordr == 2) {
			fout << msh.Elements[i][0]+1 << "\t" << msh.Elements[i][3]+1 << "\t" << msh.Elements[i][5]+1 << "\n";
			fout << msh.Elements[i][3]+1 << "\t" << msh.Elements[i][1]+1 << "\t" << msh.Elements[i][4]+1 << "\n";
			fout << msh.Elements[i][5]+1 << "\t" << msh.Elements[i][4]+1 << "\t" << msh.Elements[i][2]+1 << "\n";
			fout << msh.Elements[i][4]+1 << "\t" << msh.Elements[i][5]+1 << "\t" << msh.Elements[i][3]+1 << "\n";
		}
    	}

    	fout.close();
    	std::cout << "Tecplot file exported: " << fname << std::endl;
}
