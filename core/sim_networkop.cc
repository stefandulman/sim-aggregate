
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

MatrixXf sim_networkop::eigenvectors_kempe(int k, int maxloops) {
  
  // algorithm DecentralizedOI(k):
  // 1. choose a random k-dimensional vector Qi
  MatrixXf Q;
  Q = MatrixXf::Random(_top->size(), k);
  //Q = MatrixXf::Constant(_top->size(), k, 1);
    
  // 2. loop
  for (int i=0; i<maxloops; ++i) {
    
    // 3.   set Vi = sum_{j\in nbrs} aij Qj
    MatrixXf V = mt_sum(Q);
    
    // 4.   compute Ki = Vi' * Vi
    vector<MatrixXf> Ki = map_rowprod(V); 
    
    // 5.   K = pushsum(B, Ki)
    MatrixXf temp = _top->getnnbrs();  
    //MatrixXf B(_top->size(),1);
    //for (int j=0; j<_top->size(); ++j) {
    //  B(j,0) = 1 / (temp(j,0)+1);
    //}     
    //MatrixXf K = mt_sum(Ki, B);
    MatrixXf K = mtv_sum(Ki);
    
    //cout << "matrix K" << endl;
    //cout << K << endl;
    
    // 6.   Cholesky factorization K = R' * R
    MatrixXf R = K.llt().matrixU();

    //cout << "matrix R" << endl;
    //cout << R << endl;
    
    //exit(0);

    
    // 7.   Set Qi = Vi*R^{-1}
    Q = V * R.inverse();
    
  // 8. end loop
  }
  
  // return Qi as the ith component of each eigenvector
  return Q;
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


MatrixXf sim_networkop::mtv_sum(vector<MatrixXf> k) {
  MatrixXf res = MatrixXf::Zero(k[0].rows(), k[0].cols());
  for (int i=0; i<k.size(); ++i) {
    res += k[i];
  }
  return res;
}


// todo - the nbrs returns a list including the node, while the counting of neighbors returns a list 
// excluding the node in question - make it uniform!!!
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



