
/** 
  
  \file
  
 */

#include "sim_header.h"

#ifndef __sim_gossip_pushpull_h
#define __sim_gossip_pushpull_h

/** push-pull gossip communication algorithm
 */
template <class DataType> 
class sim_gossip_pushpull : public sim_commalg<DataType> {

  public:
    
    /** basic constructor 
      \param _t pointer to a topology object */
    sim_gossip_pushpull(sim_topology *_t) : sim_commalg<DataType>(_t) {};
    
    /** basic destructor */
    ~sim_gossip_pushpull() {};

    /** updates the topology information with a different topology object 
      \param _t pointer to a topology object */
    void changetopology(sim_topology *_t) {
      // change topology
      this->t = _t;
      // change the vector size
      this->v.resize(_t->size());
    }
    
    
    /** runs one round of the algorithm.
      \return true if any value on any node changed
     */
    bool runround(void) {
      
      bool flag = false;
      
      // for each node
      for (int i=0; i<sim_commalg<DataType>::t->size(); ++i) {
        // get random neighbor
        int nbr = sim_commalg<DataType>::t->randnbr(i);
        // perform data exchange
        if ((this->v[i]).mix(&this->v[nbr], &this->v[i], &this->v[nbr])) {
          flag = true;
        }
      }
      
      // advance time
      for (int i=0; i<sim_commalg<DataType>::t->size(); ++i) {
          (this->v[i]).step();
      }
      
      return flag;
    }


    double getvalue(int nodeid) { return sim_commalg<DataType>::v[nodeid].getvalue(); };

    vector<double> getvalue(void) {
      vector<double> res(sim_commalg<DataType>::v.size(),0);
      for (int i=0; i<sim_commalg<DataType>::v.size(); ++i) {
        res[i] = sim_commalg<DataType>::v[i].getvalue();
      }
      return res;
    }
};

#endif

