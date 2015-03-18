/** 
  
  \file
  
  \todo add tuples - chenge the mtv_sum and map_rowprod functions with template functions 
  \todo check if eigen provides alternatives for element-wise division on matrixes

 */
 
#include "sim_header.h"

#ifndef __SIM_NETWORKOP_H
#define __SIM_NETWORKOP_H


/** network operation class. 
 * this should be the base class for any network optimization operations.
 */
 
class sim_networkop {

    sim_topology *_top;

  public:

    // basic functionality
    sim_networkop(sim_topology *top = NULL) : _top(top) {};
    ~sim_networkop() {};
    
    // state-based network operations
    MatrixXf eigenvectors_kempe(int k=1, int maxloops=10);   // returns the first k eigenvectors via kempe's method
    MatrixXf fiedlervector_bertrand(int maxloops=30);        // returns the fiedler's vector via bertrand method
    
    MatrixXf eigenvectors_wave_sahai(int k=1);             // returns the first k eigenvectors via sahai's method
    MatrixXf eigenvectors_wave_franchescelli(int k=1);     // returns the first k eigenvectors via franchescelli's method
    MatrixXf eigenvectors_pia_korada(int k=3);             // returns the first k eigenvectors via korada's method (mdsmap)
    
    MatrixXf mult_lap(MatrixXf x);  // multiplication with laplacean
    MatrixXf mult_cm(MatrixXf x);   // multiplication with connectivity matrix
    
    // todo: change these to template functions to support somehow any type of local data structures
    MatrixXf mtv_sum(vector<MatrixXf> k);        // summing algorithm - input are sets of matrices and weights
    vector<MatrixXf> map_rowprod(MatrixXf m);    // returns a vector of matrixes constructed from products of rows 

    // stateless network-based operations (always return matrices)
    MatrixXf min(MatrixXf x);
    MatrixXf max(MatrixXf x);
    MatrixXf avg(MatrixXf x);
    MatrixXf sum(MatrixXf x);
    
    // local arithmetic operations - provided for convenience
    MatrixXf power(MatrixXf x, double p);
    
    // element-wise helper functions (check if Eigen provides alternatives)
    MatrixXf elwise_prod(MatrixXf a, MatrixXf b);
    MatrixXf elwise_div(MatrixXf a, MatrixXf b);
    
    // combined functions
    MatrixXf norm(MatrixXf x);
    
    // other helpers
    void print(string s, MatrixXf x);
};

#endif

