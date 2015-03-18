

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
    MatrixXf V = mult_cm(Q);
    
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
    
  // 2. compute alfa (as a factor of the maximum connectivity)
  MatrixXf alfa = max(_top->getnnbrs());
  
  // loop
  for (int i=0; i<maxloops; ++i) {
  
    // 3. compute the multiplication with the transformed laplacean
    MatrixXf v = x - elwise_div(mult_lap(x), alfa);
  
    // 4. (in parallel with 5) compute the norm of the previous result
    MatrixXf v1 = norm(v);
  
    // 5. compute the average of v
    MatrixXf v2 = avg(v);
  
    // 6. update x
    x = elwise_div(v - v2, v1);
  }
  
  return x;
}


MatrixXf sim_networkop::mult_cm(MatrixXf x) {
  MatrixXf res = MatrixXf::Zero(x.rows(), x.cols());
  
  for (int i=0; i<x.rows(); ++i) {
    list<int> nbrs = _top->getnbrsnode(i);
    
    for (list<int>::iterator it=nbrs.begin(); it!=nbrs.end(); ++it) {
      for (int j=0; j<x.cols(); ++j) {
          res(i,j) += x(*it,j);
      }
    }    
  }
  
  return res;
}


MatrixXf sim_networkop::mult_lap(MatrixXf x) {
  MatrixXf res = MatrixXf::Zero(x.rows(), x.cols());
  
  for (int i=0; i<x.rows(); ++i) {
    list<int> nbrs = _top->getnbrsnode(i);
    
    for (list<int>::iterator it=nbrs.begin(); it!=nbrs.end(); ++it) {
      for (int j=0; j<x.cols(); ++j) {
        if (i==*it) {
          res(i,j) += (nbrs.size()-1) * x(*it,j);
        } else {
          res(i,j) += -x(*it,j);
        }
      }
    }    
  }
  
  return res;
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

MatrixXf sim_networkop::power(MatrixXf x, double p) {
  MatrixXf res(x.rows(), x.cols());
  for (int i=0; i<x.rows(); ++i) {
    for (int j=0; j<x.cols(); ++j) {
      res(i,j) = pow(x(i,j),p);
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

MatrixXf sim_networkop::elwise_prod(MatrixXf a, MatrixXf b) {
  MatrixXf res(a.rows(), a.cols());
  for (int i=0; i<a.rows(); ++i) {
    for (int j=0; j<a.cols(); ++j) {
      res(i,j) = a(i,j) * b(i,j);
    }
  }
  return res;
}

MatrixXf sim_networkop::elwise_div(MatrixXf a, MatrixXf b) {
  MatrixXf res(a.rows(), a.cols());
  for (int i=0; i<a.rows(); ++i) {
    for (int j=0; j<a.cols(); ++j) {
      res(i,j) = a(i,j) / b(i,j);
    }
  }
  return res;
}


MatrixXf sim_networkop::norm(MatrixXf x) {
  return power(sum(power(x, 2)), 0.5);
}

void sim_networkop::print(string s, MatrixXf x) {
  cout << s << endl << x << endl;
}

