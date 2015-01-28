
/** 
  
  \file
  
 */

#include "sim_header.h"

#ifndef __sim_core_h
#define __sim_core_h

/** template class for overall simulation. start your simulation by deriving a class from this one.
 */ 
template <class DataType, class AlgType>
class sim_core {

  protected:

  public:
    sim_core() {};
    
    virtual vector<double> run(int nrounds) = 0;     // function returns a double value on the nodes
};

#endif
