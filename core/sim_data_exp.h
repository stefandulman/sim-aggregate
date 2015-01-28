
/** 
  
  \file
  
 */

#include "sim_header.h"

#ifndef __sim_data_exp_h
#define __sim_data_exp_h

/** class holds "new" gossiping algorithm. data type is a vector of exponentially distributed random variables. 
  mixing computes the position-wise minimum between the vectors. 
 */
class sim_data_exp : public sim_data {

    vector<double> vals;
    double origlambda;
    int nvals;

  public:
  
    sim_data_exp() {};
    void init(double _val, int _nvals);
    
    #ifdef __debug_expfunctionality
    void mix(sim_data *d2, sim_data *res1);
    #endif
        
    bool mix(sim_data *d2, sim_data *res1, sim_data *res2);

    /** empty implementation as this is class is not using time steps */
    void step(void) {}; 

    double getvalue(void);                      // returns a computed expected value
    double getlambda(void);                     // returns the orig lambda used to generate the exponential data
    int getnv(void);                            // returns the number of values
    vector<double> getvector(void);        // returns the stored vector
    double getvector(int pos);                  // returns the value of a certain position in the vector
    void setvalue(double _val);                 // changes the original value and triggers recomputation of exp variables
    void setvector(vector<double> v);     // forces the stored vector to change completely - the origval is not touched
    void setvector(int pos, double val);        // changes a position in the vector
};

#endif
