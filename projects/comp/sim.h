
#include "../../core/sim_header.h"

#ifndef __sim_h
#define __sim_h

struct res_struct {
  double err;
  int nruns;
};

enum inittypes_t {UNIFORM, NORMAL, THREE};

class sim_exp {

    sim_topology top;
    sim_gossip_pushpull<sim_data_exp> alg;
      
    double compsum;
        
  public: 
    
    sim_exp(int _nnodes, int _nsquares, int _nvals, int _maxrange) : top(_nnodes, _nsquares), alg(&top) {
      
      // initialize the data structures - random values between 0 and MAXRANGE
      alg.v.resize(top.size());
      std::vector<double> v(top.size(), 0);
      for (int i=0; i<v.size(); ++i) {
        v[i] = 1+(double)(generator() % 1000000) / 1000000.0 * _maxrange;
      }
            
      // record the mean value
      compsum = csum(v);
     
      // initialize local vectors
      for (int i=0; i<v.size(); ++i) {
        alg.v[i].init(v[i], _nvals);
      }      
    };
    

    // returns the number of runs till value converges
    struct res_struct run() {
    
      std::list<double> ressum;
       
      int count = 0;
      bool go_on = true;
      while (go_on) {
        go_on = alg.runround();
        count++;
      }
      
      compsum = alg.getvalue(0) / compsum;

      struct res_struct res;
      res.err = compsum;
      res.nruns = count;
      
      return res;
    };
};


class sim_double {

    sim_topology top;
    sim_gossip_pushpull<sim_data_double> alg;
      
    double compmean;
    double comperr; 
    double precision;
    
  public: 
    
    sim_double(int _nnodes, int _nsquares, inittypes_t _inittype, double _precision) : top(_nnodes, _nsquares), alg(&top), precision(_precision) {
      
      // initialize the data structures - random values between 0 and MAXRANGE
      alg.v.resize(top.size());
      std::vector<double> v(top.size(), 0);
      
      switch (_inittype) {
        case UNIFORM:
          for (int i=0; i<v.size(); ++i) {
            v[i] = 1 + (double)(generator() % 1000000) / 1000000.0;
          }
          break;
        case NORMAL:
          for (int i=0; i<v.size(); ++i) {
            const double s = .1;
            const double m1 = 0;
            const double m2 = 0;
            const double C = 10;
            v[i] = C * exp( -(top.getpx(i)-m1)*(top.getpx(i)-m1)/2/s/s -(top.getpy(i)-m2)*(top.getpy(i)-m2)/2/s/s); 
          }
          break;
        case THREE:
          v[0] = 0;
          v[1] = 0;
          v[2] = 10;
          break;
        default: ;
      };
      
      // record the positions and values of the nodes
      ofstream f;
      f.open("./projects/comp/values.txt");
      for (int i=0; i<v.size(); ++i) {
        f << top.getpx(i) << "\t" << top.getpy(i) << "\t" << v[i] << endl;
      }
      f.close();
      
      // record the mean and stddev values
      compmean = cmean(v);
      cout << "computed mean: " << compmean << endl;
      
      // initialize local vectors
      for (int i=0; i<v.size(); ++i) {
        alg.v[i].init(v[i]);
      }
      
      cout << "vectors initialized" << endl;
    }; 

    std::vector<double> getnorm(const std::vector<double> &v, const double &x) {
      int n = v.size();
      std::vector<double> res(n, 0);
      
      for (int i=0; i<n; ++i) {
        res[i] = v[i] / x;
      }
      
      return res;
    }  


    // returns the number of runs till value converges
    int run_time(const std::string &name) {

      std::list<double> resstd;
      double laststd;

      cout << "here 1" << endl;

      laststd = cstd(getnorm(alg.getvalue(), compmean));
      resstd.push_back(laststd);

      int cnt = 0;
      while (laststd > precision) {
        
        alg.runround();
        laststd = cstd(getnorm(alg.getvalue(), compmean));
        resstd.push_back(laststd);
     
        cout << "step " << cnt << " laststd " << laststd << endl;
     
        cnt++;
      }
      cout << "number of steps: " << cnt << endl;
      
      csave(resstd, name);
      
      return cnt;
    };
    
    
    // returns the number of runs till value converges
    int run_hops() {

      std::list<double> resstd;
      double laststd;

      laststd = cstd(getnorm(alg.getvalue(), compmean));
      resstd.push_back(laststd);

      int cnt = 0;
      while (laststd > precision) {
        
        alg.runround();
        laststd = cstd(getnorm(alg.getvalue(), compmean));
        resstd.push_back(laststd);
     
        cnt++;
      }

      return cnt;
    };

};


#endif

