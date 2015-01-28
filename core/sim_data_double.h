
/** 
  
  \file
  
 */

#include "sim_header.h"

#ifndef __sim_data_double_h
#define __sim_data_double_h

/** class mapping on basic "old" gossip algorithms. basic value is a double value. mixing is based on simple averaging. 
 */
class sim_data_double : public sim_data {

    double val;

  public:
  
    sim_data_double();
    void init(double _val);
    bool mix(sim_data *d2, sim_data *res1, sim_data *res2);
    
    /** empty implementation as this is class is not using time steps */
    void step(void) {}; 

    double getvalue(void);
    void setvalue(double _val);
};

#endif

