/** 
  
  \file
  
  \todo add touples
  
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
    
    // matrix and vector based operations
//    template <typename InType, typename OutType>
//    vector<OutType> map(OutType (*fcnptr)(InType, InType), vector<InType> inv);
  
    vector<MatrixXf> map_rowprod(MatrixXf m);                // returns a vector of matrixes constructed from products of rows 
        
    MatrixXf mt_sum(MatrixXf m);                             // matrix topology sum - sum of the neighboring values
    MatrixXf mtv_sum(vector<MatrixXf> k);         // summing algorithm - input are sets of matrices and weights
    
    // state-based network operations
    MatrixXf eigenvectors_kempe(int k=1, int maxloops=10);   // returns the first k eigenvectors via kempe's method
    sim_matrix eigenvectors_wave_sahai(int k=1);             // returns the first k eigenvectors via sahai's method
    sim_matrix eigenvectors_wave_franchescelli(int k=1);     // returns the first k eigenvectors via franchescelli's method
    
    sim_matrix eigenvectors_pia_korada(int k=3);             // returns the first k eigenvectors via korada's method (mdsmap)
    sim_matrix fiedlervector(void);                          // returns the fiedler's vector via bertrand method
        
};

#endif

