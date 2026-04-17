// ======================================================
// File: sparse.cpp
// ======================================================
#include "sparse.h"
#include <algorithm>
#include <iostream>

CSR convert_to_csr(const SparseRowMap &Amap) {
    int n = Amap.n;
    CSR A(n);
    A.row_ptr.assign(n+1, 0);

    // Count total nonzeros
    int nnz = 0;
    for (int i=0;i<n;++i) { 
	    nnz += (int)Amap.rows[i].size();
    }
    //std::cout << nnz << "\n";


    A.col_idx.reserve(nnz);
    A.values.reserve(nnz);

    int ptr = 0;
    A.row_ptr[0] = 0;

    // Convert each row
    for (int i=0;i<n;++i) {
        // Copy unordered_map entries into vector for sorting
        std::vector<std::pair<int,double>> entries;
        entries.reserve(Amap.rows[i].size());
        for (auto &kv : Amap.rows[i]) {
	     	entries.push_back(kv);
	}

	// sort entries
        std::sort(entries.begin(), entries.end(),
                 [](auto &a, auto &b){ return a.first < b.first; });

        // Push into CSR storage
        for (auto &e : entries) {
            A.col_idx.push_back(e.first);
            A.values.push_back(e.second);
            ++ptr;
        }

        A.row_ptr[i+1] = ptr;
    }

    return A;
}
