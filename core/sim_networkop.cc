
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
  
/*  // algorithm DecentralizedOI(k):
  // 1. choose a random k-dimensional vector Qi
  MatrixXf Q;
  Q = MatrixXf::Random(_top->size(), k);
  
  // 2. loop
  for (int i=0; i<maxloops; ++i) {
    
    // 3.   set Vi = sum_{j\in nbrs} aij Qj
    MatrixXf V = mt_sum(Q);
    
    // 4.   compute Ki = Vi' * Vi
//    vector<simmatrix> Ki;
//    Ki = map(product_rows, V, top);
    vector<MatrixXf> Ki = map_rowprod(V); 
    
    // 5.   K = pushsum(B, Ki)
    MatrixXf B = MatrixXf::Zero(k,k);  // replace with something meaningful
    MatrixXk K = Ki.sum();             // replace with something meaningful
    
    // 6.   Cholesky factorization K = R' * R
    
    
    // 7.   Set Qi = Vi*R^{-1}
    
  // 8. end loop
  }
  
  // return Qi as the ith component of each eigenvector
*/  
  
  sim_matrix res(k);
  return res;
}



// matrix and vector based operations
//template <typename InType, typename OutType>
//vector<OutType> map(OutType (*fcnptr)(InType, InType), vector<InType> inv) {
//  vector<OutType> res;
//  for (int i=0; i<inv.size(); ++i) {
//    res.push_back((*fcnptr)(inv[i], inv[0]));
//  }  
//  return res;
//}



vector<MatrixXf> sim_networkop::map_rowprod(MatrixXf m) {
  vector<MatrixXf> res;
  for (int i=0; i<m.rows(); ++i) {
    MatrixXf temp = m.row(i).transpose() * m.row(i);
    res.push_back(temp);
  }
  return res;
}



MatrixXf sim_networkop::mt_sum(MatrixXf m) {
  
  // assertion - network size and matrix size need to be combined somehow, now they are duplicated
  if (m.rows()!=_top->size()) {
    cout << "sim_networkop::nbr_sum: mismatched network size and matrix size" << endl;
    exit(0);
  }
  
  MatrixXf res = MatrixXf::Zero(m.rows(), m.cols());
  for (int i=0; i<m.rows(); ++i) {
    list<int> nbrs = _top->getnbrsnode(i);
    for (list<int>::iterator it=nbrs.begin(); it!=nbrs.end(); ++it) {
      for (int j=0; j<m.cols(); ++j) {
        res(i,j) += m(*it,j);
      }
    }
  }
  
  return res;
}



