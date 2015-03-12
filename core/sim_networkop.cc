
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
    //MatrixXf temp = _top->getnnbrs();  
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



// returns the fiedler vector via bertrand's method
MatrixXf sim_networkop::fiedlervector_bertrand(int maxloops) {

  // 1. choose a random k-dimensional vector x
  MatrixXf x = MatrixXf::Random(_top->size(), 1);
  
  // 2. compute alfa (as a factor of the maximum connectivity
  double alfa = max(_top->getnnbrs());
  
  // loop
  for (int i=0; i<maxloops; ++i) {
  
    // 3. compute the multiplication with the transformed laplacean
    MatrixXf v = x - lap_mult(x)/alfa;
  
    // 4. (in parallel with 4) compute the norm of the previous result
    MatrixXf v1 = norm(v);
  
    // 5. compute the average of v
    MatrixXf v2 = avg(v);
  
    // 6. update x
    x = v.cWiseProduct(v1) - v2/v1*MatrixXf::Constant(_top->size(),1,1);
  }
  
  return x;
}



MatrixXf sim_networkop::lap_mult(MatrixXf x) {
  MatrixXf res = MatrixXf::Zero(x.rows(), x.cols());
  for (int i=0; i<x.rows(); ++i) {
    list<int> nbrs = _top->getnbrsnode(i);
    for (list<int>::iterator it=nbrs.begin(); it!=nbrs.end(); ++it) {
      for (int j=0; j<x.cols(); ++j) {
        if (i==*it) {
          res(i,j) += nbrs.size() * x(*it,j);
        } else {
          res(i,j) += -x(*it,j);
        }
      }
    }
  }
  return res;
}

double sim_networkop::norm(MatrixXf x) {
  return x.norm();
}



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


MatrixXf sim_networkop::min(MatrixXf x) {
  MatrixXf res = MatrixXf::Constant(x.rows(), x.cols(),1);
  for (int i=0; i<x.cols(); ++i) {
    res.col(i) *= x.col(i).minCoeff();
  }
  return res;
}

MatrixXf sim_networkop::max(MatrixXf x) {
  MatrixXf res = MatrixXf::Constant(x.rows(), x.cols(),1);
  for (int i=0; i<x.cols(); ++i) {
    res.col(i) *= x.col(i).maxCoeff();
  }
  return res;
}

MatrixXf sim_networkop::sum(MatrixXf x) {
  MatrixXf res = MatrixXf::Constant(x.rows(), x.cols(),1);
  for (int i=0; i<x.cols(); ++i) {
    res.col(i) *= x.col(i).sum();
  }
  return res;
}

MatrixXf sim_networkop::avg(MatrixXf x) {
  MatrixXf res = MatrixXf::Constant(x.rows(), x.cols(),1);
  for (int i=0; i<x.cols(); ++i) {
    res.col(i) *= x.col(i).mean();
  }
  return res;
}


