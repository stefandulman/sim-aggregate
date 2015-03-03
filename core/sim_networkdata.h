/** 
  
  \file
  
  \todo add touples
  
 */
 
#include "sim_header.h"

#ifndef __SIM_NETWORKDATA_H
#define __SIM_NETWORKDATA_H

/** network data storage class. 
 * think of it as a generalized matlab vector, on which operations can be defined.
 * this should be the base class for any network optimization operations.
 */
class sim_networkdata {

  public:
    
    
    
    // basic functionality
    
    sim_networkdata();
    ~sim_networkdata();
    
    void setup();
    
    std::vector<double> get_data();           // getting data out
    double              get_data(int _pos);
    
    void set_data(std::vector<double> _data); // setting data
    void set_data(int _pos, double _data);
    
    void reinit(void);                        // force reinitialization of the data structures
    
    void advance_time(int _nsteps);           // advance time with so many steps
    
    double value(int _pos = -1);      // returns a double value from a certain position in the vector
    std::vector<double> vvalue(void); // returns a vector of values 
    
    
    
    // basic math operations
    
    sim_networkdata     sum(int _power = 1, int _nsteps = -1);
    sim_networkdata     min(int _nsteps = -1);
    sim_networkdata     count(int _nsteps = -1);
    


    // overloaded operators
    sim_networkdata& operator/=(const sim_networkdata& rhs) {
        return *this;
    }
    sim_networkdata operator/(const sim_networkdata& rhs) {
        return *this /= rhs;
    }
    
};

#endif

