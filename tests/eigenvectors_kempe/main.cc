
#include "../../core/sim_header.h"



int main(void) {

  cout << "test: eigenvectors_kempe" <<endl;

  // generate known topology
  int n = 7;
  MatrixXf pos(7,2);
  
  pos(0,0) = 0.1; pos(0,1) = 0.1;
  pos(1,0) = 0.5; pos(1,1) = 0.1;
  pos(2,0) = 0.9; pos(2,1) = 0.1;
  pos(3,0) = 0.1; pos(3,1) = 0.5;
  pos(4,0) = 0.1; pos(4,1) = 0.9;
  pos(5,0) = 0.1; pos(5,1) = 0.9;
  pos(6,0) = 0.1; pos(6,1) = 0.9;
  
  sim_topology top(7,3,pos);
  
  // assert the connectivity matrix
  cout << "  connectivity matrix... ";
  MatrixXf temp = top.getnnbrs();
  MatrixXf exp_temp(7,1);
  exp_temp << 2, 3, 1, 5, 3, 3, 3;
  for (int i=0; i<7; ++i) { 
    if (temp(i,0)!=exp_temp(i,0)) {
      cout << "eigenvectors_kempe error: wrong connectivity matrix" << endl;
      cout << "expected: " << endl << exp_temp << endl;
      cout << "got:      " << endl << temp << endl;
      return(0);
    }
  }
  cout << "passed" << endl;
 
  // run the eigen decomposition according to kempe
  cout << "  decomposition... ";
  sim_networkop sim(&top);
  MatrixXf res = sim.eigenvectors_kempe(3,100);  
  MatrixXf exp_res(7,3);
  exp_res << 0.2499,  0.4601, -0.5333, \
             0.2691,  0.6085,  0.1317, \
             0.0820,  0.3837,  0.7777, \
             0.5512,  0.1211, -0.2221, \
             0.4300, -0.2922,  0.1213, \
             0.4300, -0.2922,  0.1213, \
             0.4300, -0.2922,  0.1213;
  // compute maximum difference
  for (int i=0; i<7; ++i) {
    for (int j=0; j<3; ++j) {
      double d = fabs(res(i,j) - exp_res(i,j));
      if (d>0.1) {
        cout << "eigenvectors_kempe error: wrong eigenvalues matrix" << endl;
        cout << "expected:" << endl << exp_res << endl;
        cout << "received:" << endl << res << endl;
      }
    }
  }
  cout << "passed" << endl;
  
  return 0;  
}

