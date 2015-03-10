
#include "../../core/sim_header.h"
#include "../../Eigen/Sparse"

int main(void) {

  MatrixXf a(3,3);
  a << 2, -1, 0, -1, 2, -1, 0, -1, 2;
  MatrixXf l = a.llt().matrixL();
  cout << l << endl;

/*  
  vector<double> inv(10,1);
  vector<double> outv = map(myfcn, inv);
  for (int i=0; i<inv.size(); ++i) {
    cout << outv[i] << " ";
  }
  cout << endl;
*/

/*
  int nnodes = 10;
  sim_networkdata x(10), temp(10);
  
  temp = x.sum() / x.count();
  double avg = temp.value();
  
  cout << "average: " << avg << endl;
  
  x = x.sum(2) / x.sum();
  cout << "min val: " << x.min().value() << endl;
  
  cout << "done!" << endl;
*/

/*  
  MatrixXf m(3,3), q, v;

  m(0,0) = 0; m(0,1) = 1; m(0,2) = 2; 
  m(1,0) = 1; m(1,1) = 1; m(1,2) = 1;
  m(2,0) = 2; m(2,1) = 1; m(2,2) = 2;

  q = MatrixXf::Random(3,3);
  
  for (int i=0; i<100; ++i) {
    v = m*q;
    Eigen::HouseholderQR<MatrixXf> temp(v);
    q = temp.householderQ();
  }
  cout << "Iterated Q matrix is:\n" << q << endl;

  Eigen::SelfAdjointEigenSolver<MatrixXf> mf(m);
  cout << "The ideal eigenvector matrix is:\n" << mf.eigenvectors() << endl;
  cout << "The eigenvalues are:\n" << mf.eigenvalues() << endl;
*/

  return 0;  
}

