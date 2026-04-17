// ============================================
// File: linSolve.h
// ============================================

#ifndef LINSOLVE_H
#define LINSOLVE_H
#include <vector>
#include <unordered_map>

#include "sparse.h"

// Dot product calculator
double dot(const std::vector<double>& a,
           const std::vector<double>& b);

// Matrix Operator
void OperateYPlusAlphaX(      std::vector<double> &y,
          	   	const std::vector<double> &x, 
		   	      double a);

// Jacobi Preconditioner
void applyJacobi(const std::vector<double>& Minv,
                 const std::vector<double>& r,
                       std::vector<double>& z);

// Check matrix symmetry
bool isSymmetric(const CSR& A, 
		       double tol);

// Gauss Elimination solver
bool GaussElimination(const SparseRowMap &A,
                      const std::vector<double> &b,
                            std::vector<double> &x);

// Conjugate Gradient
int ConjugateGradient(const CSR& A,
                      const std::vector<double>& b,
                            std::vector<double>& x,
                            double tol,
                            int maxIter);

// BiCGStab Solver
int BiCGSTAB(const CSR& A,
             const std::vector<double>& b,
             	   std::vector<double>& x,
             	   double tol,
             	   int maxIter);

#endif
