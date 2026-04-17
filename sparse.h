// ======================================================
// File: sparse.h — Sparse storage + CSR converter
// ======================================================
#ifndef SPARSE_H
#define SPARSE_H
#include <vector>
#include <unordered_map>

// A preliminary sparse storage using row-wise hash maps.
// This is good for FEM assembly where entries appear in arbitrary order.
struct SparseRowMap {
    int n;
    std::vector<std::unordered_map<int,double>> rows;

    SparseRowMap(int n_=0): n(n_), rows(n_) {}

    // Accumulate A(i,j) += v
    void add(int i,int j,double v) {
        rows[i][j] += v;
    }
};


// CSR = Compressed Sparse Row format
// (row_ptr, col_idx, values)
struct CSR {
    int n;                         // matrix dimension
    std::vector<int> row_ptr;     // row_ptr[i] = index in col_idx where row i begins
    std::vector<int> col_idx;     // column index for each nonzero
    std::vector<double> values;   // values for each nonzero

    CSR() : n(0) {}
    CSR(int n_): n(n_), row_ptr(n_+1,0) {}

    
    // Compute y = A*x
    std::vector<double> mul(const std::vector<double>& x) const {
        std::vector<double> y(n,0.0);
        for (int i=0;i<n;++i) {
            for (int k=row_ptr[i]; k<row_ptr[i+1]; ++k)
                y[i] += values[k] * x[col_idx[k]];
        }
        return y;
    }
    
};


// Convert SparseRowMap (easy to assemble) → CSR (efficient for solve)
CSR convert_to_csr(const SparseRowMap &Amap);

#endif
