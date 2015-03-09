
/** 
  
  \file
  
 */
 
#include "sim_header.h"

#ifndef __SIM_TOPOLOGY_H
#define __SIM_TOPOLOGY_H

/** topology storage class
 */
class sim_topology {

    int nsq;                               // number of squares per line
    vector<int> nbrs;                      // array of number of neighbors 
    vector<int> ninzone;                   // number of nodes in a zone

    vector<double> px;                     // position x
    vector<double> py;                     // position y

    unordered_map<int, vector<int>> htn;   // nodes storage per zone - hash table nodes

    int nnodes;                            // number of nodes 
    int nzones;                            // number of zones

    int getzone(int nodeid);               // returns the zone for a node
    int getzone(int x, int y);             // returns an index from an x, y coordinate
    
    int getx(int zone) { return zone/(nsq+2) - 1; };
    int gety(int zone) { return zone%(nsq+2) - 1; };
    
    list<int> getnbrszone(int zoneid);

    int getnnbr(int nodeid, int posid);

  public:

    // returns a list with the neighbors of this node
    list<int> getnbrsnode(int nodeid);

    /** basic constructor
      \param _nnodes number of nodes
      \param _nsq number of grid steps on one direction
     */
    sim_topology(int _nnodes, int _nsq);

    /** function returns a random neighbor or the node itself if it is alone.
      \param nodeid node identifier 
      \return the id of a random neighbor or the node itself if alone */
    int randnbr(int nodeid);  

    /** getter for the number of nodes in the topology 
      \return the number of nodes in the topology */
    int size(void) { return nnodes; };
    
    /** getter for the number of grid steps in one direction 
      \return the number of grid steps in one direction
     */
    int getnzones(void) { return nzones; }; 
    
    /** getter for the x position of a node
      \param nodeid node identifier
      \return the x position of the specified node
     */
    double getpx(int nodeid) { return px[nodeid]; };
    
    /** getter for the y position of a node
      \param nodeid node identifier 
      \return the y position of the specified node
     */
    double getpy(int nodeid) { return py[nodeid]; };
    
    /** \todo - move all data saving functions somewhere else or redefine them to be consistent */
    
    /** functions saves the positions of all nodes to a file
      \param name file path */
    void savepos(string name);
    
    /** functions saves all the links from the topology to a file
      \param name file path */
    void savelinks(string name);
    
    // debug functions
    #ifdef __debug_topology_printall
    void debug(void);
    #endif
};

#endif

