
#include "sim_topology.h"


// returns the zone for an index of a node
int sim_topology::getzone(int nodeid) {
  
  #ifdef __debug_rangecheck
  if (nodeid>=nnodes) {
    cout << __PRETTY_FUNCTION__ << " nodeid is too large: " << nodeid << ">=" << nnodes << endl;
    exit(EXIT_FAILURE);
  }
  #endif

  return getzone(px[nodeid]*nsq, py[nodeid]*nsq);
}

// returns an index from an x, y coordinate
int sim_topology::getzone(int x, int y) {

  #ifdef __debug_rangecheck
  if ((x>nsq) || (y>nsq) || (x<-1) || (y<-1)) {
    cout << __PRETTY_FUNCTION__ << " wrong index: " << x << "," << y << endl;
    exit(EXIT_FAILURE);
  }
  #endif

  // shift all nodes with one position diagonally
  return (nsq+2) * (x+1) + (y+1);
}



sim_topology::sim_topology(int _nnodes, int _nsq) : nnodes(_nnodes), nsq(_nsq) {
    
    // sanity check
    if (nsq<=0) {
      cout << __PRETTY_FUNCTION__ << " number of squares is too small: " << nsq << endl;
      exit(EXIT_FAILURE);
    }
    
    // total number of zones is (nsq+2)*(nsq+2)
    nzones = (nsq+2)*(nsq+2);
    
    //resize data structures
    nbrs.resize(nnodes);
    ninzone.resize(nzones);
    px.resize(nnodes);
    py.resize(nnodes);
    
    // generate random positions for the nodes
    for (int i=0; i<nnodes; ++i) {    
      px[i] = ((double) generator() / generator.max());
      py[i] = ((double) generator() / generator.max());
    }
    
    // figure out the zones in which nodes are - O(n)
    unordered_map<int, list<int>> m; // temporary variable
    for (int i=0; i<nnodes; ++i) {
      int index = getzone(i);
      ninzone[index]++;      // increment counter for the nodes
      m[index].push_back(i); // record this node
    }
    
    // transform the map(list) into a map(vector)
    for (int i=0; i<nzones; ++i) {
      vector<int> v(make_move_iterator(begin(m[i])), make_move_iterator(end(m[i])));
      htn[i] = v;      
    }
    
    // fill in the nbrs data structure - take each square excluding the ones on the border
    for (int i=0; i<nsq; ++i) {
      for (int j=0; j<nsq; ++j) {
        
        int sum = -1;   // the node itself does not matter

        sum = sum + ninzone[getzone(i-1,j-1)];
        sum = sum + ninzone[getzone(i-1,j  )];
        sum = sum + ninzone[getzone(i-1,j+1)];
        sum = sum + ninzone[getzone(i  ,j-1)];
        sum = sum + ninzone[getzone(i  ,j  )];
        sum = sum + ninzone[getzone(i  ,j+1)];
        sum = sum + ninzone[getzone(i+1,j-1)];
        sum = sum + ninzone[getzone(i+1,j  )];
        sum = sum + ninzone[getzone(i+1,j+1)];
        
        // update all the nodes in this zone
        int z = getzone(i,j);
        for (int ind=0; ind<htn[z].size(); ++ind) {
          nbrs[htn[z][ind]] = sum;
        }
      }
    }
}

/* ugly duplication of the constructor code - do something better than this! */
sim_topology::sim_topology(int _nnodes, int _nsq, MatrixXf _pos) : nnodes(_nnodes), nsq(_nsq) {
    // sanity check
    if (nsq<=0) {
      cout << __PRETTY_FUNCTION__ << " number of squares is too small: " << nsq << endl;
      exit(EXIT_FAILURE);
    }
    
    // total number of zones is (nsq+2)*(nsq+2)
    nzones = (nsq+2)*(nsq+2);
    
    //resize data structures
    nbrs.resize(nnodes);
    ninzone.resize(nzones);
    px.resize(nnodes);
    py.resize(nnodes);
    
    // generate random positions for the nodes
    for (int i=0; i<nnodes; ++i) {   
      px[i] = _pos(i,0);
      py[i] = _pos(i,1);
    }
    
    // figure out the zones in which nodes are - O(n)
    unordered_map<int, list<int>> m; // temporary variable
    for (int i=0; i<nnodes; ++i) {
      int index = getzone(i);
      ninzone[index]++;      // increment counter for the nodes
      m[index].push_back(i); // record this node
    }
    
    // transform the map(list) into a map(vector)
    for (int i=0; i<nzones; ++i) {
      vector<int> v(make_move_iterator(begin(m[i])), make_move_iterator(end(m[i])));
      htn[i] = v;      
    }
    
    // fill in the nbrs data structure - take each square excluding the ones on the border
    for (int i=0; i<nsq; ++i) {
      for (int j=0; j<nsq; ++j) {
        
        int sum = -1;   // the node itself does not matter

        sum = sum + ninzone[getzone(i-1,j-1)];
        sum = sum + ninzone[getzone(i-1,j  )];
        sum = sum + ninzone[getzone(i-1,j+1)];
        sum = sum + ninzone[getzone(i  ,j-1)];
        sum = sum + ninzone[getzone(i  ,j  )];
        sum = sum + ninzone[getzone(i  ,j+1)];
        sum = sum + ninzone[getzone(i+1,j-1)];
        sum = sum + ninzone[getzone(i+1,j  )];
        sum = sum + ninzone[getzone(i+1,j+1)];
        
        // update all the nodes in this zone
        int z = getzone(i,j);
        for (int ind=0; ind<htn[z].size(); ++ind) {
          nbrs[htn[z][ind]] = sum;
        }
      }
    }
}



int sim_topology::getnnbr(int nodeid, int posid) {
  
  // go through the zones
  int x = getx(getzone(nodeid));
  int y = gety(getzone(nodeid));
    
  vector<int> v = {getzone(x-1,y-1), getzone(x-1,y), getzone(x-1,y+1), getzone(x,y-1), getzone(x,y), getzone(x,y+1), getzone(x+1,y-1), getzone(x+1,y), getzone(x+1,y+1)};
  
  int sum = 0;
  for (int i=0; i<v.size(); ++i) {
    if ((sum+ninzone[v[i]])>posid) {  
      return htn[v[i]][posid-sum];
    }
    sum = sum + ninzone[v[i]];
  }
  
  cout << __PRETTY_FUNCTION__ << " invalid position of neighbor requested (" << nodeid << "," << posid << ")" << endl;
  exit(EXIT_FAILURE);
}



int sim_topology::randnbr(int nodeid) {
  
  #ifdef __debug_rangecheck
    if (nodeid>=size()) {
      cout << __PRETTY_FUNCTION__ << " node " << nodeid << " requested out of " << size() << " nodes." << endl;
      exit(EXIT_FAILURE);
    }
  #endif
  
  // trivial case - single node
  if (nbrs[nodeid] == 0) {
    return nodeid;
  }
  
  //TODO: see why the generator is not initialized with a different value here!!!
  std::chrono::high_resolution_clock::time_point t = std::chrono::high_resolution_clock::now();
  std::chrono::high_resolution_clock::duration d = t - t.min();
  generator.seed(abs(d.count()%100000));
  
  // pick number of random neighbor
  int temp = generator();
  int pos = temp % nbrs[nodeid];

  int res = getnnbr(nodeid, pos);
  if (res==nodeid) {
    if (pos==0) {
      return getnnbr(nodeid, 1);
    } else {
      return getnnbr(nodeid, pos-1);
    }
  }
    
  return res;
}



// includes the node itself - actually it returns the neighbors for all the nodes in the zone
list<int> sim_topology::getnbrsnode(int nodeid) {

  #ifdef __debug_rangecheck
    if (nodeid>=size()) {
      cout << __PRETTY_FUNCTION__ << " node " << nodeid << " requested out of " << size() << " nodes." << endl;
      exit(EXIT_FAILURE);
    }
  #endif

  return getnbrszone(getzone(nodeid));
}
  
  
  
list<int> sim_topology::getnbrszone(int zoneid) {

  #ifdef __debug_rangecheck
    if (zoneid>=nzones) {
      cout << __PRETTY_FUNCTION__ << " zone " << zoneid << " requested out of " << nzones << " zones." << endl;
      exit(EXIT_FAILURE);
    }
  #endif

  list<int> res;

  int x = getx(zoneid);
  int y = gety(zoneid);

  // create a list of interesting zones
  list<int> zones;
  zones.push_back(getzone(x-1,y-1));
  zones.push_back(getzone(x-1,y  ));
  zones.push_back(getzone(x-1,y+1));
  zones.push_back(getzone(x  ,y-1));
  zones.push_back(getzone(x  ,y  ));
  zones.push_back(getzone(x  ,y+1));
  zones.push_back(getzone(x+1,y-1));
  zones.push_back(getzone(x+1,y  ));
  zones.push_back(getzone(x+1,y+1));

  // add neighbors
  for (list<int>::iterator it=zones.begin(); it!=zones.end(); ++it) {
    for (int jt=0; jt<htn[*it].size(); ++jt) {
        res.push_back(htn[*it][jt]);
    }
  }
  
  return res;
}



void sim_topology::savelinks(string name) {
  ofstream f;
  f.open(name);
  
  vector<int> flag(nnodes, 0);
  
  for (int i=0; i<nsq; ++i) {
    for (int j=0; j<nsq; ++j) {
    
      // get the nodes in the current square
      vector<int> tnodes = htn[getzone(i,j)];
      // get all the neighbors around
      list<int> tnbrs = getnbrszone(getzone(i,j));
      
      for (int ind=0; ind<tnodes.size(); ++ind) {
        flag[tnodes[ind]] = 1;
        for (list<int>::iterator it=tnbrs.begin(); it!=tnbrs.end(); ++it) {
          if (!flag[*it]) {
            f << px[tnodes[ind]] << "\t" << py[tnodes[ind]] << endl;
            f << px[*it] << "\t" << py[*it] << endl;
            f << endl;
          }
        }
      }
      
    }
  }
  
  f.close();
}



MatrixXf sim_topology::getnnbrs(void) {
  MatrixXf res = MatrixXf::Zero(nnodes,1);
  for (int i=0; i<nnodes; ++i) {
    res(i) = nbrs[i];
  }
  return res;
}




void sim_topology::savepos(string name) {

  ofstream f;
  f.open(name);

  for (int i=0; i<nnodes; ++i) {
    f << px[i] << "\t" << py[i] << endl;
  }

  f.close();
}


#ifdef __debug_topology_printall
void sim_topology::debug(void) {
  
  cout << "positions of nodes:" << endl;
  for (int i=0; i<nnodes; ++i) {
    cout << "(" << pos(i,0) << "," << pos(i,1) << ") " << endl;
  }
  cout << endl;
  
  cout << "connectivity matrix:" << endl;
  for (int i=0; i<nnodes; ++i) {
    for (int j=0; j<nnodes; ++j) {
      cout << cm(i,j);
    }
    cout << endl;
  }
  cout << endl;
  
  cout << "nr neighbors all nodes:" << endl;
  for (int i=0; i<nnodes; ++i) {
    cout << nbrs[i] << endl;
  }
  cout << endl;
}
#endif

