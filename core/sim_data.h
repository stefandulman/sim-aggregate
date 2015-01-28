
/** 
  
  \file
  
 */
 
#include "sim_header.h"

#ifndef __sim_data_h
#define __sim_data_h

/** virtual base class for various data types
 */
class sim_data {

    int initialized;


  public:
    sim_data() { 
      initialized = 0; 
    };
    
    virtual bool mix(sim_data *d2, sim_data *res1, sim_data *res2) = 0;  // mix two values 

    virtual double getvalue(void) = 0;
    virtual void setvalue(double _val) = 0;
    
    /** explicitely advance the time with one step. makes sense for algorithms which expire 
      data structures with time. leave an empty implementation otherwise.
     */
    virtual void step(void) = 0;
    
    int isinitialized(void) { return initialized; };
    void initialize(void) { initialized = 1; };
};

#endif

