
#include "sim_networkop.h"

// returns the first k eigenvectors via kempe's method
  // algorithm DecentralizedOI(k):
  // 1. choose a random k-dimensional vector Qi
  // 2. loop
  // 3.   set Vi = sum_{j\in nbrs} aij Qj
  // 4.   compute Ki = Vi' * Vi
  // 5.   K = pushsum(B, Ki)
  // 6.   Cholesky factorization K = R' * R
  // 7.   Set Qi = Vi*R^{-1}
  // 8. end loop
  // return Qi as the ith component of each eigenvector

sim_matrix sim_networkop::eigenvectors_kempe(int k, int maxloops) {
  
  // algorithm DecentralizedOI(k):
  // 1. choose a random k-dimensional vector Qi
  sim_matrix Q(_top->size(), k);
  Q.init_random();
  
  // 2. loop
  for (int i=0; i<maxloops; ++i) {
    
    // 3.   set Vi = sum_{j\in nbrs} aij Qj
    sim_matrix V(_top->size(), k);
    V = nbr_sum(Q);
    
    // 4.   compute Ki = Vi' * Vi
    
    
    // 5.   K = pushsum(B, Ki)
    // 6.   Cholesky factorization K = R' * R
    // 7.   Set Qi = Vi*R^{-1}
  
  // 8. end loop
  }
  
  // return Qi as the ith component of each eigenvector
  
  
  
  
  sim_matrix res(k);
  return res;
}



sim_matrix sim_networkop::nbr_sum(sim_matrix m) {
  
  // assertion - network size and matrix size need to be combined somehow, now they are duplicated
  if (m.nrows()!=_top->size()) {
    cout << "sim_networkop::nbr_sum: mismatched network size and matrix size" << endl;
    exit(0);
  }
  
  sim_matrix res(_top->size(),1);
  res.init_const(0);
  
  for (int i=0; i<_top->size(); ++i) {
    list<int> nbrs = _top->getnbrsnode(i);
    for (list<int>::iterator it=nbrs.begin(); it!=nbrs.end(); ++it) {
      for (int j=0; j<m.ncols(); ++j) {
        res[i][j] += m[*it][j];
      }
    }
  }
  
  return res;
}
