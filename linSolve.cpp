// ===========================================
// File: linSolve.cpp
// ===========================================

#include <iostream>
#include <numeric>
#include <cmath>
#include <vector>
#include <random>

#include "linSolve.h"


//===========================================================
//  DOT PRODUCT
//===========================================================
double dot(const std::vector<double>& a,
           const std::vector<double>& b)
{
    double sum = 0.0;
    for (size_t i = 0; i < a.size(); ++i)
        sum += a[i] * b[i];
    return sum;
}

//===========================================================
//  Matrix Operation y = y + alpha*x
//===========================================================
void OperateYPlusAlphaX(      std::vector<double>& y,
          	   	const std::vector<double>& x, 
		   	      double a)
{
    for (size_t i = 0; i < y.size(); ++i)
        y[i] += a * x[i];
}


// ============================
// APPLY JACOBI PRECONDITIONER
// M^{-1} r = z
// ============================
void applyJacobi(const std::vector<double>& Minv,
                 const std::vector<double>& r,
                       std::vector<double>& z)
{
    for (size_t i=0; i<r.size(); i++)
        z[i] = Minv[i] * r[i];
}


// =====================================================
// Check if matrix is symmetric
// =====================================================
bool isSymmetric(const CSR& A, 
		       double tol = 1e-12)
{
    int n = A.n;

    // For each row i
    for (int i = 0; i < n; i++) {

        // For each entry A(i, j)
        for (int k = A.row_ptr[i]; k < A.row_ptr[i+1]; k++) {

            int j = A.col_idx[k];
            double Aij = A.values[k];

            // Search for A(j, i)
            bool found = false;
            for (int kk = A.row_ptr[j]; kk < A.row_ptr[j+1]; kk++) {
                if (A.col_idx[kk] == i) {
                    found = true;
                    double Aji = A.values[kk];
                    if (std::fabs(Aij - Aji) > tol)
                        return false;  // values differ
                    break;
                }
            }

            if (!found) return false; // missing symmetric entry
        }
    }

    return true;
}


// ---------------------------------------
// Gauss Elimination Solver
// ---------------------------------------
bool GaussElimination(const SparseRowMap &Amap,
                      const std::vector<double> &b_,
                            std::vector<double> &x)
{

	SparseRowMap A = Amap;
	std::vector<double> b = b_;

	int n = A.n;

	// Forward Elimination
	for (int i=0; i<n; ++i) {
		for (int j=0; j<i; ++j) {
			if (A.rows[i][j] != 0) {
				double factor = A.rows[i][j]/A.rows[j][j];
				for (int k=0; k<n; ++k) {
					A.rows[i][k] -= A.rows[j][k]*factor;
				}
				b[i] -= b[j]*factor;
			}
		}
	}


	// Backward Substitution
	for (int i=n-1; i>=0; --i) {
		x[i] = b[i];
		for (int j=i+1; j<n; ++j) {
			x[i] -= A.rows[i][j]*x[j];
		}
		x[i] = x[i]/A.rows[i][i];
	}

	return false;
}



// ==============================================
// Conjugate Gradient Solver (CG)
// A must be symmetric & SPD
// ==============================================
int ConjugateGradient(const CSR& A,
                      const std::vector<double>& b,
                       	    std::vector<double>& x,
                       	    double tol = 1e-12,
                       	    int maxIter = 10000)
{

	bool symCheck = isSymmetric(A, 1e-12);
       	if (!symCheck)	{
		std::cerr << "Matrix not Symmetric, use BiCGStab instead\n";
		return 0;
	}
    
	int n = A.n;

    	std::vector<double> r(n), p(n), Ap(n);

    	// r = b - A*x
    	r = A.mul(x);
    	for (int i = 0; i < n; i++)
        	r[i] = b[i] - r[i];

    	p = r;    // initial search direction

    	double rs_old = dot(r, r);

    	for (int iter = 1; iter <= maxIter; iter++) {

        	// Ap = A*p
        	Ap = A.mul(p);

        	double alpha = rs_old / dot(p, Ap);

        	// Update x
        	for (int i = 0; i < n; i++)
            		x[i] += alpha * p[i];

        	// Update r
        	for (int i = 0; i < n; i++)
            		r[i] -= alpha * Ap[i];

        	double rs_new = dot(r, r);

        	// Convergence check
        	if (std::sqrt(rs_new) < tol) {
            		return iter;
        	}

        	double beta = rs_new / rs_old;

        	// Update search direction
        	for (int i = 0; i < n; i++)
            		p[i] = r[i] + beta * p[i];

        	rs_old = rs_new;
    	}

    	return maxIter;
}



//===========================================================
//  BiCGSTAB SOLVER FOR CSR MATRIX
//===========================================================
//
// Solves A*x = b
// A : CSR matrix
// b : RHS vector
// x : initial guess (input) + solution (output)
// tol : tolerance for convergence
// maxIter : max number of iterations
//
// Returns: number of iterations used
//
//===========================================================
int BiCGSTAB(const CSR& A,
             const std::vector<double>& b,
             std::vector<double>& x,
             double tol = 1e-12,
             int maxIter = 10000)
{
    	int n = A.n;

    	std::vector<double> r(n), r0(n), p(n), v(n), s(n), t(n);
	std::vector<double> z(n), z_t(n);

    	// r0 = b - A*x
    	r = A.mul(x);
    	for (int i = 0; i < n; ++i)
        	r[i] = b[i] - r[i];


	// initial shadow residual
	r0 = r;

	// ---------------------------------------
    	// BUILD JACOBI PRECONDITIONER
    	// Minv[i] = 1 / A(i,i)
    	// ---------------------------------------
    	std::vector<double> Minv(n, 0.0);

    	for (int i=0; i<n; i++) {
        	for (int k = A.row_ptr[i]; k < A.row_ptr[i+1]; k++) {
            		if (A.col_idx[k] == i) {
                		Minv[i] = 1.0 / A.values[k];
                		break;
            		}
        	}
    	}

    	double rho_old = 1, alpha = 1, omega = 1;

    	p.assign(n, 0.0);
    	v.assign(n, 0.0);

    	double normb = std::max(1e-16, std::sqrt(dot(b,b)));
    	double resid = std::sqrt(dot(r,r)) / normb;
    	if (resid < tol) return 0;

    	//-------------------------------------------------------
    	// MAIN ITERATION LOOP
    	//-------------------------------------------------------
    	for (int iter = 1; iter <= maxIter; ++iter)
    	{
        	double rho = dot(r0, r);
		int update = 0;
        	while (std::fabs(rho) < 1e-30) {
            		std::cout << update++ << " - BiCGSTAB: breakdown (rho=0) - Updating Shadow residual\n";
			// update shadow residual
                	// Create a random number generator
                	std::random_device rd;
                	std::mt19937 gen(rd());
                	std::uniform_real_distribution<> dis(-1.0, 1.0); // Uniform distribution between -1 and 1
                	// randomize shadow residual
                	for (int i = 0; i < n; ++i) {
                        	r0[i] += dis(gen);
                	}
			if (update>5) {
            			std::cerr << "BiCGSTAB: breakdown (rho=0)\n";
            			return iter;
			}
        	}

        	double beta = (rho / rho_old) * (alpha / omega);

        	// p = r + beta*(p - omega*v)
        	for (int i = 0; i < n; ++i)
            		p[i] = r[i] + beta * (p[i] - omega * v[i]);

		// Apply Jacobi: p_hat = M^{-1} p
        	applyJacobi(Minv, p, z);

        	// v = A*p
        	v = A.mul(z);

        	alpha = rho / dot(r0, v);

        	// s = r - alpha*v
        	for (int i = 0; i < n; ++i)
            		s[i] = r[i] - alpha * v[i];

        	// Check early convergence
        	if (std::sqrt(dot(s, s)) / normb < tol) {
            		OperateYPlusAlphaX(x, z, alpha);
            		return iter;
        	}

		// Apply Jacobi: s_hat = M^{-1} s
        	applyJacobi(Minv, s, z_t);

        	// t = A*s
        	t = A.mul(z_t);

        	omega = dot(t, s) / dot(t, t);
        	if (std::fabs(omega) < 1e-16) {
            		std::cerr << "BiCGSTAB: breakdown (omega=0)\n";
            		return iter;
        	}

        	// x = x + alpha*p + omega*s
        	for (int i = 0; i < n; ++i)
            		x[i] += alpha * z[i] + omega * z_t[i];

        	// r = s - omega*t
        	for (int i = 0; i < n; ++i)
            		r[i] = s[i] - omega * t[i];

        	// compute residual
        	resid = std::sqrt(dot(r, r)) / normb;
        	if (resid < tol)
            		return iter;

        	rho_old = rho;
    	}

    	std::cerr << "BiCGSTAB: max iterations reached\n";
    	return maxIter;
}
