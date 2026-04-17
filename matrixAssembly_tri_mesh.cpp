// ====================================================
// File: matrixAssembly
// 1. Shape Function and derivatives
// 2. Gauss Quadrature Points and Weights
// 3. Assemble Global Matrices
// 4. Apply Dirichlet BCs
// ====================================================

#include <cmath>
#include <array>
#include <vector>
#include <iostream>

#include "matrixAssembly.h"
#include "shapeFunction.h"

//1. Shape Functions and derivatives
static inline void shapeFunction(double xi, double eta, int Ordr,
			 std::vector<double>& N,
			 std::vector<std::array<double,2>>& dN)
{
	double tau = 1-xi-eta;

	// Shape Functions
	if (Ordr==1) {
		N[0] = tau;
		N[1] = xi;
		N[2] = eta;
	}
	else {
		N[0] = tau*(2*tau-1);
        	N[1] = xi*(2*xi-1);
       		N[2] = eta*(2*eta-1);
        	N[3] = 4.0*tau*xi;
        	N[4] = 4.0*xi*eta;
        	N[5] = 4.0*eta*tau;
	}


	// Derivatives wrt xi and eta
	if (Ordr==1) {
		dN[0][0] = -1.0;	dN[0][1] = -1.0;
		dN[1][0] =  1.0;	dN[1][1] =  0.0;
		dN[2][0] =  0.0;	dN[2][1] =  1.0;
	}
	else {
		dN[0][0] =  1.0-4.0*tau;   dN[0][1] =  1.0-4.0*tau;
        	dN[1][0] =  4.0*xi-1.0;    dN[1][1] =  0.0;
        	dN[2][0] =  0.0;           dN[2][1] =  4.0*eta-1.0;
        	dN[3][0] =  4.0*(tau-xi);  dN[3][1] = -4.0*xi;
        	dN[4][0] =  4.0*eta;       dN[4][1] =  4.0*xi;
        	dN[5][0] = -4.0*eta;       dN[5][1] =  4.0*(tau-eta);
	}

}



//2. Gauss Quadrature Points and Weights (for n-order polynomial integration)
static inline void gaussQuad(int n) {

	// Adjusting size of Gauss Quadrature  Arrays
	gaussPts.resize(n);
        gaussWts.resize(n);

	// Linear polynomial (n=1)
	if (n<=1) {
		gaussPts = {{1.0/3.0, 1.0/3.0}};
		gaussWts = {1.0/2.0};
	}

	// Quadratic polynomial (n=2)
        if (n==2) {
                gaussPts = {{1.0/6.0, 1.0/6.0}, 
			    {1.0/6.0, 2.0/3.0}, 
			    {2.0/3.0, 1.0/6.0}};
                gaussWts = {1.0/6.0, 
			    1.0/6.0, 
			    1.0/6.0};
        }

	// Cubic polynomial (n=3)
        if (n==3) {
                gaussPts = {{1.0/3.0, 1.0/3.0}, 
			    {  0.6, 0.2  }, 
			    {  0.2, 0.6  }, 
			    {  0.2, 0.2  }};
                gaussWts = {-27.0/96.0, 
			     25.0/96.0, 
			     25.0/96.0, 
			     25.0/96.0};
	}

	// BiQuadratic polynomial (n=4)
        if (n==4) {
                gaussPts = {{0.108103018, 0.445948491},
   			    {0.445948491, 0.108103018},
   			    {0.445948491, 0.445948491},
   			    {0.816847573, 0.091576214},
   			    {0.091576214, 0.816847573},
   			    {0.091576214, 0.091576214}};
                gaussWts = {0.111690794839005,
                	    0.111690794839005,
                	    0.111690794839005,
                	    0.054975871827661,
                	    0.054975871827661,
                	    0.054975871827661};
	}

}



//3. Assemble Global Matrices
void assemble_matrices(const Mesh &msh,
			     Flow &flow)
{

	double ReInv = 1.0/flow.Re;

	int n = msh.n_nodes;
	SparseRowMap Kmap;

	// initialize Matrices
	flow.Ktau.resize(n,std::vector<double>(n));
        flow.fg.resize(n,std::vector<double>(2));

	// initialize sparse matrix
	Kmap = SparseRowMap(n);
	flow.K_tau = CSR(n);

	// Selecting higher order polynomial for Gauss Quadrature (M_u)
	int gaussOrdr = msh.Ordr;
	gaussQuad(gaussOrdr);

	// Looping over Element Matrices
	for (int e=0; e<msh.n_elements; ++e) {
		auto conn = msh.Elements[e];	// Element nodal connectivity

		// Physical Coordinates of element nodes
		std::vector<std::pair<double,double>> xy;
		xy.resize(3*msh.Ordr);
		for (int a=0; a<3*msh.Ordr; ++a) {
			xy[a] = msh.Coords[conn[a]];
		}


		// Local matrices (defined for higher order)
		double Ke_tau[6][6] 	= {{0}};
		double fe_g[6][2]	= {{0}};

		// Quadrature Loop
		int ng = gaussPts.size();
		for (int gp=0; gp<ng; ++gp) {
			double xi	= gaussPts[gp][0];
			double eta	= gaussPts[gp][1];

			// Evaluate Shape Functions
			std::vector<double> N;
			std::vector<std::array<double,2>> dN;
			N.resize(3*msh.Ordr);
			dN.resize(3*msh.Ordr);
			shapeFunction(xi, eta, msh.Ordr, N, dN);

			// Gauss Coordinates
                        double x=0, y=0;
                        for (int a=0; a<3*msh.Ordr; ++a) {
                                x += N[a]*xy[a].first;
				y += N[a]*xy[a].second;
                        }

			// Compute Jacobian Matrix
			double J[2][2] = {{0.0,0.0},{0.0,0.0}};
			for (int a=0; a<3*msh.Ordr; ++a) {
				J[0][0] += dN[a][0] * xy[a].first;
				J[0][1] += dN[a][0] * xy[a].second;
				J[1][0] += dN[a][1] * xy[a].first;
				J[1][1] += dN[a][1] * xy[a].second;
			}
			

			// Determinant of Jacobian
			double detJ = J[0][0]*J[1][1] - J[0][1]*J[1][0];

			// Inverse of Jacobian
			double invJ[2][2];
			invJ[0][0] =  J[1][1]/detJ;
			invJ[1][0] = -J[1][0]/detJ;
			invJ[0][1] = -J[0][1]/detJ;
			invJ[1][1] =  J[0][0]/detJ;

			// Natural Gradients of Shape Functions dN/dx_i
			std::vector<std::array<double, 2>> dNdx;
			dNdx.resize(3*msh.Ordr);
			for (int a=0; a<3*msh.Ordr; ++a) {
				dNdx[a][0] = invJ[0][0]*dN[a][0] + invJ[0][1]*dN[a][1];
				dNdx[a][1] = invJ[1][0]*dN[a][0] + invJ[1][1]*dN[a][1];
			}


			double wDetJ = gaussWts[gp]*detJ;

			// Assemble Element Matrices and Vectors
			for (int i=0; i<3*msh.Ordr; ++i) {
				for (int j=0; j<3*msh.Ordr; ++j) {
				
					// Ke_tau Matrix
					// Int(dNdx'*dNdx)dV
					Ke_tau[i][j] 	+= ReInv * (dNdx[i][0]*dNdx[j][0] + dNdx[i][1]*dNdx[j][1]) * wDetJ;
				}

				// fe_g Vector
				// Int(N'*g)dV
				fe_g[i][0] += (N[i]*flow.g(x,y,0)) * wDetJ;
				fe_g[i][1] += (N[i]*flow.g(x,y,1)) * wDetJ;
			}

		} // Quadrature Loop end
		
		// Scatter into global maps
                for (int i=0; i<3*msh.Ordr; ++i) {
                    	int gi = conn[i];
                      	for (int j=0; j<3*msh.Ordr; j++) {
                             	int gj = conn[j];
                            	flow.Ktau[gi][gj] +=  Ke_tau[i][j];
				if (Ke_tau[i][j] != 0) {
					Kmap.add(gi,gj,Ke_tau[i][j]);
				}
                      	}
                    	flow.fg[gi][0] += fe_g[i][0];
                    	flow.fg[gi][1] += fe_g[i][1];
		}

	} // Element matrix-vector loop end
	
	flow.K_tau = convert_to_csr(Kmap);

}



// Apply Dirichlet boundary conditions
void apply_dirichlet( const Mesh &msh,
                            SparseRowMap &Amap,
                            std::vector<double> &b,
                            std::vector<std::vector<int>> &DirichletNodes,
                            std::vector<std::vector<double>> &DirichletValues)
{
	int n = Amap.n;

	// Loop over individual Dirichlet Nodes
	int index;
        for (int dir=0; dir<2; ++dir) {
              	for (int j=0; j<DirichletNodes[dir].size(); ++j) {
                     	int iNode = DirichletNodes[dir][j];
                    	index = iNode + (n/2)*dir;
                     	for (int i=0; i<n; ++i) {
				if (i==index) continue;
				auto it = Amap.rows[i].find(index);
				if ( it!=Amap.rows[i].end() ) {
					double val = it->second;
					b[i] -= val*DirichletValues[dir][j];
					Amap.rows[i].erase(it);
				}
                     	}
			Amap.rows[index].clear();
                     	Amap.rows[index][index] = 1.0;
                      	b[index] = DirichletValues[dir][j];
             	}
       	}


	/*
	// Loop over individual Dirichlet Nodes
    	for (int k=0; k<nDirichlet; ++k) {
        	int index = DirichletIndex[k];
		// Zero out column k in all other rows
		for (int j=0; j<nDirichlet; ++j) {
			if (j==index) continue;
			auto it = Amap.rows[j].find(index);
			if (it!=Amap.rows[j].end() ) {
				double val = it->second;
				b[j] -= val*DirichletValue[k];
				Amap.rows[j].erase(it);
			}
		}

		// replace row k by identity
		Amap.rows[index].clear();
		Amap.rows[index][index] = 1.0;
		b[index] = DirichletValue[k];

        }
	*/
}
