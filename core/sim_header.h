
// main header file


/** 
  
  \file

  \todo doxygen - add text description to all base members
  
  \todo doxygen - add text desciption to all classes
  
  \todo unify notation   get_value, getValue
  
  \todo split data types and mixing into two different classes under the data types
  
  \todo figure out how to store the vectors of data types in connection to the above point
  
  \todo figure out how to add algorithms on top of the basic data types

  \todo transform topology class into a base class and allow derivations of various topologies

  \todo in function sim_topology::randnbr() - define a common function (template) for printing error text and exiting
  
  \todo modify all the error messages to look nice
  
  \todo add ifdef guards for debug
  
  \todo modify topology to not use a cm matrix at all
  
  \todo modify topology to do computation of cm faster than o(n^2)
  
  \todo optimize cm away in topology
  
  \todo move data collection in sim_core
  
  \todo fix the public declaration in sim_alg
  
  \todo fix initialization and access to v in algorithm
  
  \todo replace all the double comparisons with a delta comparison
  
  \todo gossipexpttl example outputs a damaged pos.eps file
  
  \todo change the position generatos in topology to use the proper random functions
  
  \todo add paramater protection to all functions in topology_fast
  
  \todo optimization idea - detect zones with no changes and ignore those nodes for the alg.step() function
  
  \todo possible bug in the topology recomputation - on a network with diamter 14, values converge in 15 steps... secondly, after killing half of the network, the global sum does not jump to the proper value

 */


#ifndef __sim_header_h
#define __sim_header_h

  // debug messages
  
  // standard checking of ranges in all function calls
//  #define __debug_rangecheck
  
  //#define __debug_matrix_printall
  //#define __debug_topology_printall
  //#define __debug_topology_printdist
  //#define __debug_simdata_randomsequence
  //#define __debug_pushpull_trackminimum
  //#define __debug_dataexp_checkblending
  //#define __debug_expfunctionality
  //#define __debug_nbrcomm

  // includes - standard
  #include <stdlib.h>
  #include <iostream>
  #include <sstream>
  #include <fstream>
  #include <cstdlib>
  #include <algorithm>
  #include <vector>
  #include <unordered_map>
  #include <list>
  #include <random>
  #include <chrono>
  #include <ctime>
  #include <float.h>

  using namespace std;
  
  // includes Eigen library
  #include "../Eigen/Dense"
  using Eigen::MatrixXf;
  
  // includes - local
  #include "sim_topology.h"
  #include "sim_data.h"
  #include "sim_commalg.h"
  #include "sim_core.h"
  #include "sim_data_double.h"
  #include "sim_gossip_pushpull.h"
  #include "sim_data_exp.h"
  #include "sim_data_exp_fast.h"
  #include "sim_matrix.h"
  #include "sim_networkop.h"

  // algorithms customization
  #define __custom_dataexpttl_smooth

  // random number generator - visible from any file, does not repeat
  namespace {
      default_random_engine generator;
  }
  
  // print a vector to screen
  inline void printv(const vector<double> &v) {
    for (unsigned int i=0; i<v.size(); ++i) {
      cout << v[i] << " ";
    }
    cout << endl;
  }
  
  // execute a command in the shell and return result
  inline std::string exec(string cmd) {

      FILE* pipe = popen(cmd.c_str(), "r");
      if (!pipe) return "ERROR";
      char buffer[128];
      std::string result = "";
      while(!feof(pipe)) {
      	if(fgets(buffer, 128, pipe) != NULL)
      		result += buffer;
      }
      pclose(pipe);
      return result;
  }


    /** convenience function computing the sum of values in the network
      \param v vector holding the values in the network
      \return sum of the values in the network
     */
    inline double csum(const vector<double> &v) {
      double m = 0;
      for (int i=0; i<v.size(); ++i) {
        m += v[i];
      }
      return m;
    }
    
    inline double csum(const list<double> &v) {
      double m = 0;
      for (list<double>::const_iterator it=v.begin(); it!=v.end(); ++it) {
        m += *it;
      }
      return m;
    }
    
    
    /** convenience function computing the mean value of values in the network 
      \param v vector holding the values in the network
      \return mean value of all values in the network */
    inline double cmean(const vector<double> &v) {
      return csum(v) / v.size();
    }
    
    inline double cmean(const list<double> &v) {
      return csum(v) / v.size();
    }
    
    /** convenience function computing the variance of the values in the network
      \param v vector holding the values in the network
      \return variance for all values in the network */
    inline double cvariance(const vector<double> &v) {
      double m = cmean(v);
      double m2 = 0;
      for (int i=0; i<v.size(); ++i) {
        m2 += v[i]*v[i];
      }
      m2 = m2 / v.size();
      return (m2 - m*m);
    }
    
    inline double cvariance(const list<double> &v) {
      double m = cmean(v);
      double m2 = 0;
      for (list<double>::const_iterator it=v.begin(); it!=v.end(); ++it) {
        m2 += (*it) * (*it);
      }
      m2 = m2 / v.size();
      return (m2 - m*m);
    }
    
    
    /** convenience function computing the standard deviation of the values in the network
      \param v vector holding the values in the network
      \return standard deviation for all values in the network */
    inline double cstd(const vector<double> &v) {
      return sqrt(cvariance(v));
    }
    
    inline double cstd(const list<double> &v) {
      return sqrt(cvariance(v));
    }

    /** convenience function for saving a list to a file */
    inline void csave(const std::list<double> &l, const std::string &name) {
      ofstream f;
      f.open(name);
      int i=0;
      for (std::list<double>::const_iterator it = l.begin(); it!=l.end(); ++it) {
        f << i++ << "\t" << *it << endl;
      }
      f.close();
    }
    
    /** convenience function for saving a vector to a file */
    inline void csave(const std::vector<double> &v, const std::string &name) {
      ofstream f;
      f.open(name);
      for (int i=0; i<v.size(); ++i) {
        f << i << "\t" << v[i] << endl;
      }
      f.close();
    }

    /** convenience function for saving a value to a file */
    inline void csave(const double &v, const std::string &name) {
      ofstream f;
      f.open(name);
      f << 0 << "\t" << v << endl;
      f.close();
    }

#endif

