
/** 
  
  \file
  
 */

#include "sim_header.h"

#ifndef __sim_data_exp_fast_h
#define __sim_data_exp_fast_h

/** class adds self-stabilization properties to the "new" gossiping algorithm. 
  data type is a vector of exponentially distributed random variables seconded by a time-to-live vector. 
  mixing computes the position-wise minimum between the vectors. upin changing the global minimum value, it is marked as negative and a value removal mechanism is put into place.  
 */
class sim_data_exp_fast : public sim_data_exp {

    bool checkexpired(int index);
    
    vector<double> origvals;
    vector<int> ttl;
    
    int maxttl;
    
    void defaultinit(double _val, int _nvals, int _maxttl);

  public:
  
    sim_data_exp_fast() {};                                     
    void init(double _val, int _nvals, int _maxttl);            

    bool mix(sim_data *d2, sim_data *res1, sim_data *res2 = NULL);
    
    
    /** implementation of the time stepping mechanism */
    void step(void);                                          

    double getvalue(void);
    
    int getttl(int pos);
    vector<int> getttl(void);
    void setttl(vector<int> newttl);
    void setttl(int pos, int val);
    
    int getmaxttl(void) { return maxttl; };
    void updatevalttl(int pos, double val, int ttl);
    vector<double> getorigvector(void);
    double getorigvector(int pos);
};

#endif
