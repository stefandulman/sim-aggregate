
/** 
  
  \file
  
 */

#include "sim_header.h"

#ifndef __sim_commalg_h
#define __sim_commalg_h

/** template class for communication algorithms
 */
template <class DataType>
class sim_commalg {
  
  public:
  
    sim_topology *t;
    vector<DataType> v;
  
//  public:
  
    /** basic class constructor */
    sim_commalg(sim_topology *_t) : t(_t) {};             
    
    /** basic class desctructor */
    virtual ~sim_commalg() {};
    
    /** function performs one computation step */
    virtual bool runround(void) = 0;
        
    /** convenience function to retrieve a value on a node 
      \param nodeid node identifier
      \return a double value from obtained from the data stored on the node 
     */
    virtual double getvalue(int nodeid) = 0;
    
    /** convenience function to retrieve all values in the network as a vector
      \return vector of double values
     */
    virtual vector<double> getvalue(void) = 0;

};

#endif

