
#include "sim_data_exp.h"  
  
void sim_data_exp::init(double _val, int _nvals) {
  
  nvals = _nvals;
  origlambda = _val;
  
  vals.resize(nvals);
  
  exponential_distribution<double> distribution(origlambda);
  
  for (int i=0; i<nvals; ++i) {
    vals[i] = distribution(generator);
  }

  initialize();
}


#ifdef __debug_expfunctionality
void sim_data_exp::mix(sim_data *d2, sim_data *res1) {
  
  #ifdef __debug_rangecheck
  if (!isinitialized()) {
    cout << "sim_data_exp::mix() - data structure not initialized. call init." << endl;
    exit(EXIT_FAILURE);
  }
  #endif
  
  sim_data_exp temp;
  temp.init(1, getnv());
  mix(d2, res1, &temp);
}
#endif


bool sim_data_exp::mix(sim_data *d2, sim_data *res1, sim_data *res2) {
  
  #ifdef __debug_rangecheck
  if (!isinitialized()) {
    cout << "sim_data_exp::mix() - data structure not initialized. call init." << endl;
    exit(EXIT_FAILURE);
  }
  #endif

  bool resval = false;
  
  for (int i=0; i<nvals; ++i) {
    double v1 = vals[i];
    double v2 = ((sim_data_exp *)d2)->getvector(i);
    double res = ((v1>v2) ? v2 : v1);
    if (v1!=v2) {
      resval = true;
    }
    ((sim_data_exp *)res1)->setvector(i, res);
    ((sim_data_exp *)res2)->setvector(i, res);
  }
  
  return resval;
}


double sim_data_exp::getvalue(void) {
  
  #ifdef __debug_rangecheck
  if (!isinitialized()) {
    cout << "sim_data_exp::mix() - data structure not initialized. call init." << endl;
    exit(EXIT_FAILURE);
  }
  #endif
  
  double res = 0;
  for (int i=0; i<nvals; ++i) {
    res = res + vals[i];
  }
  res = ((double)nvals) / res;
  return res;
}
  

double sim_data_exp::getlambda(void) {
  
  #ifdef __debug_rangecheck
  if (!isinitialized()) {
    cout << "sim_data_exp::mix() - data structure not initialized. call init." << endl;
    exit(EXIT_FAILURE);
  }
  #endif
  
  return origlambda;
}


vector<double> sim_data_exp::getvector(void) {
  
  #ifdef __debug_rangecheck
  if (!isinitialized()) {
    cout << "sim_data_exp::mix() - data structure not initialized. call init." << endl;
    exit(EXIT_FAILURE);
  }
  #endif
  
  return vals;
}


void sim_data_exp::setvalue(double _val) {
  
  #ifdef __debug_rangecheck
  if (!isinitialized()) {
    cout << "sim_data_exp::mix() - data structure not initialized. call init." << endl;
    exit(EXIT_FAILURE);
  }
  #endif
  
  init(_val, nvals); 
}


void sim_data_exp::setvector(vector<double> v) {
  
  #ifdef __debug_rangecheck
  if (!isinitialized()) {
    cout << "sim_data_exp::mix() - data structure not initialized. call init." << endl;
    exit(EXIT_FAILURE);
  }
  #endif
  
  vals = v;
}


int sim_data_exp::getnv(void) {
  
  #ifdef __debug_rangecheck
  if (!isinitialized()) {
    cout << "sim_data_exp::mix() - data structure not initialized. call init." << endl;
    exit(EXIT_FAILURE);
  }
  #endif
  
  return nvals;
}


double sim_data_exp::getvector(int pos) {
  
  #ifdef __debug_rangecheck
  if (!isinitialized()) {
    cout << "sim_data_exp::mix() - data structure not initialized. call init." << endl;
    exit(EXIT_FAILURE);
  }
  #endif

  #ifdef __debug_rangecheck
  if (pos>=nvals) {
    cout << "sim_data_exp::getvector() error: pos exceeds vector size " << pos << ">=" << nvals << endl;
    exit(EXIT_FAILURE);
  }
  #endif

  return vals[pos];
}

void sim_data_exp::setvector(int pos, double val) {
  
  #ifdef __debug_rangecheck
  if (!isinitialized()) {
    cout << "sim_data_exp::mix() - data structure not initialized. call init." << endl;
    exit(EXIT_FAILURE);
  }
  #endif
  
  #ifdef __debug_rangecheck
  if (pos>=nvals) {
    cout << "sim_data_exp::setvector() error: pos exceeds vector size " << pos << ">=" << nvals << endl;
    exit(EXIT_FAILURE);
  }
  #endif
  
  vals[pos] = val;
}

