
#include "../../core/sim_header.h"



int main(void) {

  cout << "test: fiedlervector_bertrand" <<endl;

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
      cout << "fiedlervector_kempe error: wrong connectivity matrix" << endl;
      cout << "expected: " << endl << exp_temp << endl;
      cout << "got:      " << endl << temp << endl;
      return 0;
    }
  }
  cout << "passed" << endl;
 
  // run the eigen decomposition according to kempe
  cout << "  decomposition... ";
  sim_networkop sim(&top);
  MatrixXf res = sim.fiedlervector_bertrand(100);  
  MatrixXf exp_res(7,1);
  exp_res << 0.1145, 0.2968, 0.7349, -0.1360, -0.3367, -0.3367, -0.3367;
  
  
  // compute maximum difference
  for (int i=0; i<7; ++i) {
    double d = fabs(res(i,0) - exp_res(i,0));
    if (d>0.1) {
      cout << "fiedlervector_bertrand error: wrong fiedler vector" << endl;
      cout << "expected:" << endl << exp_res << endl;
      cout << "received:" << endl << res << endl;
      return 0;
    }
  }
  cout << "passed" << endl;
  
  return 0;  
}

