
/** 
  
  \todo - resval might not be updated correctly in mix() - if a timer expires, it stays the same
  
 */

  

#include "sim_data_exp_fast.h"  
  
void sim_data_exp_fast::defaultinit(double _val, int _nvals, int _maxttl) {
  maxttl = _maxttl;
  // call original init
  sim_data_exp::init(_val, _nvals);
  // refresh ttl values
  vector<int> tempttl(_nvals, _maxttl);
  setttl(tempttl);
  // store the default origvector
  origvals = getvector();
}
  
  
void sim_data_exp_fast::init(double _val, int _nvals, int _maxttl) {

  if (isinitialized() && (getnv()==_nvals) && (getmaxttl()==_maxttl)) {
  
    //cout << "--old vector: " << getvector(0) << endl;
  
    // store original configuration
    sim_data_exp_fast temp = *this; 
  
    // default init
    defaultinit(_val, _nvals, _maxttl);
    
    //cout << "--old vector: " << getvector(0) << endl;
    
    // mix the values
    mix(&temp, this);

    //cout << "--old vector: " << getvector(0) << endl;
    
    // mark old values as negative
    for (int i=0; i<getnv(); ++i) {
      if ((temp.getorigvector(i)==getvector(i))) {
        setvector(i, -fabs(getvector(i)));
      }
    }

    //cout << "--old vector: " << getvector(0) << endl;
    //cout << "---origval: " << temp.getorigvector(0) << endl;
    //cout << "---newval:  " << getorigvector(0) << endl;
    
    
    return;
  }
  
  // not initialized - do the regular stuff
  defaultinit(_val, _nvals, _maxttl);
}


bool sim_data_exp_fast::mix(sim_data *d2, sim_data *res1, sim_data *res2) {
  
  #ifdef __debug_rangecheck
  if (!isinitialized()) {
    cout << "sim_data_exp_fast::mix() - data structure not initialized. call init." << endl;
    exit(EXIT_FAILURE);
  }
  #endif
  
  const double C = 2;
  bool resval = false;
  
  for (int i=0; i<getnv(); ++i) {
    
    // make sure we do not work with expired data
    if (checkexpired(i)) {
      resval = true;
    }
    if (((sim_data_exp_fast *)d2)->checkexpired(i)) {
      resval = true;
    }
    
    double v1 = getvector(i);
    double v2 = ((sim_data_exp_fast *)d2)->getvector(i);
    
    // mark change if values are different
    if (v1!=v2) {
      resval = true;
    } else {
      // also if one of them is negative
      if (v1<0) {
        resval = true;
      }
    }
    
    int t1 = getttl(i);
    int t2 = ((sim_data_exp_fast *)d2)->getttl(i);    
    
    // positive values - simple
    if ((v1>0) && (v2>0)) {
      if (v1>v2) {
        // update the values and ttl
        ((sim_data_exp_fast *)res1)->updatevalttl(i, v2, t2 - 1);
        if (res2!=NULL) ((sim_data_exp_fast *)res2)->updatevalttl(i, v2, t2);
      } else {
        if (v1<v2) {
          // update the values & ttl
          ((sim_data_exp_fast *)res1)->updatevalttl(i, v1, t1);
          if (res2!=NULL) ((sim_data_exp_fast *)res2)->updatevalttl(i, v1, t1 - 1);    
        } else {
          // values are equal - pick largest ttl
          if (t1 < t2) {
            ((sim_data_exp_fast *)res1)->updatevalttl(i, v1, t2 - 1);
            if (res2!=NULL) ((sim_data_exp_fast *)res2)->updatevalttl(i, v1, t2);
          } else {
            if (t1 > t2) {
              ((sim_data_exp_fast *)res1)->updatevalttl(i, v1, t1);
              if (res2!=NULL) ((sim_data_exp_fast *)res2)->updatevalttl(i, v1, t1 - 1);
            } else {
              // ttls are equal - do not change
              ((sim_data_exp_fast *)res1)->updatevalttl(i, v1, t1);
              if (res2!=NULL) ((sim_data_exp_fast *)res2)->updatevalttl(i, v1, t2);
            }
          }        
        }
      }
      continue;
    }
    
    // one negative, one positive
    if ((v1<0) && (v2>0)) {
      if ((-v1)==v2) {
        ((sim_data_exp_fast *)res1)->updatevalttl(i, v1, getmaxttl());
        if (res2!=NULL) ((sim_data_exp_fast *)res2)->updatevalttl(i, v1, getmaxttl());   // initiate new flood seed
      } else {
        ((sim_data_exp_fast *)res1)->updatevalttl(i, v1, t1/C);          // halve both values
        if (res2!=NULL) ((sim_data_exp_fast *)res2)->updatevalttl(i, v1, t1/C);
      }
      continue;
    }
    
    // symmetric
    if ((v1>0) && (v2<0)) {
      if (v1==(-v2)) {
        ((sim_data_exp_fast *)res1)->updatevalttl(i, v2, getmaxttl());
        if (res2!=NULL) ((sim_data_exp_fast *)res2)->updatevalttl(i, v2, getmaxttl());
      } else {
        ((sim_data_exp_fast *)res1)->updatevalttl(i, v2, t2/C);
        if (res2!=NULL) ((sim_data_exp_fast *)res2)->updatevalttl(i, v2, t2/C);
      }
      continue;
    }
    
    // both values are negative
    if (v1==v2) {
      ((sim_data_exp_fast *)res1)->updatevalttl(i, v1, t1/C); // make either of them half
      if (res2!=NULL) ((sim_data_exp_fast *)res2)->updatevalttl(i, v2, t2/C);
      continue;
    }
    
    if (v1<v2) {
      ((sim_data_exp_fast *)res1)->updatevalttl(i, v1, t1/C);
      if (res2!=NULL) ((sim_data_exp_fast *)res2)->updatevalttl(i, v1, t1/C);
    } else {
      ((sim_data_exp_fast *)res1)->updatevalttl(i, v2, t2/C);
      if (res2!=NULL) ((sim_data_exp_fast *)res2)->updatevalttl(i, v2, t2/C);
    }
  }
  
  return resval;
}


bool sim_data_exp_fast::checkexpired(int index) {

    bool resval = false;
  
    // value expired, replace
    if (getttl(index) <= 0) {
      if (getvector(index)!=getorigvector(index)) {
        resval = true;
      }
      updatevalttl(index, getorigvector(index), getmaxttl());
    }
    
    return resval;
}


void sim_data_exp_fast::step(void) {
  
  #ifdef __debug_rangecheck
  if (!isinitialized()) {
    cout << "sim_data_exp_fast::step() - data structure not initialized. call init." << endl;
    exit(EXIT_FAILURE);
  }
  #endif
  
  for (int i=0; i<getnv(); ++i) {
    
    // is it our own value?
    //if (fabs(getvector(i))==getorigvector(i)) {
    if (getvector(i)==getorigvector(i)) {
      setttl(i, getmaxttl());
      continue;
    }
    
    // decrease ttl
    setttl(i, getttl(i) - 1);
    
    // check if it expired
    checkexpired(i);
  }
}



double sim_data_exp_fast::getvalue(void) {
  #ifdef __debug_rangecheck
  if (!isinitialized()) {
    cout << "sim_data_exp::mix() - data structure not initialized. call init." << endl;
    exit(EXIT_FAILURE);
  }
  #endif
  
  double res = 0;
  for (int i=0; i<getnv(); ++i) {
    res = res + fabs(getvector(i));
  }
  res = ((double)getnv()) / res;
  return res;
}


int sim_data_exp_fast::getttl(int pos) {
  #ifdef __debug_rangecheck
  if (!isinitialized()) {
    cout << "sim_data_exp_ttl::getttl() - data structure not initialized. call init." << endl;
    exit(EXIT_FAILURE);
  }
  #endif
  
  #ifdef __debug_rangecheck
  if (pos>=getnv()) {
    cout << "sim_data_exp_ttl::getttl() error - position too large " << pos << ">=" << getnv() << endl;
    exit(EXIT_FAILURE);
  }
  #endif
  
  return ttl[pos];  
}

vector<int> sim_data_exp_fast::getttl(void) {
  #ifdef __debug_rangecheck
  if (!isinitialized()) {
    cout << "sim_data_exp_ttl::getttl() - data structure not initialized. call init." << endl;
    exit(EXIT_FAILURE);
  }
  #endif
  
  return ttl;  
}

void sim_data_exp_fast::setttl(vector<int> newttl) {
  #ifdef __debug_rangecheck
  if (!isinitialized()) {
    cout << "sim_data_exp_ttl::setttl() - data structure not initialized. call init." << endl;
    exit(EXIT_FAILURE);
  }
  #endif
  
  ttl = newttl;  
}

void sim_data_exp_fast::setttl(int pos, int val) {
  #ifdef __debug_rangecheck
  if (!isinitialized()) {
    cout << "sim_data_exp_ttl::setttl() - data structure not initialized. call init." << endl;
    exit(EXIT_FAILURE);
  }
  #endif
  
  #ifdef __debug_rangecheck
  if (pos>=getnv()) {
    cout << "sim_data_exp_ttl::setttl() error - position too large " << pos << ">=" << getnv() << endl;
    exit(EXIT_FAILURE);
  }
  #endif
  
  ttl[pos] = val;
}

void sim_data_exp_fast::updatevalttl(int pos, double val, int ttl) {
  
  #ifdef __debug_rangecheck
  if (!isinitialized()) {
    cout << "sim_data_exp_ttl::updatevalttl() - data structure not initialized. call init." << endl;
    exit(EXIT_FAILURE);
  }
  #endif
  
  #ifdef __debug_rangecheck
  if (pos>=getnv()) {
    cout << "sim_data_exp_ttl::updatevalttl() error - position too large " << pos << ">=" << getnv() << endl;
    exit(EXIT_FAILURE);
  }
  #endif
  
  setvector(pos, val);
  setttl(pos, ttl);
}

double sim_data_exp_fast::getorigvector(int pos) {
  
  #ifdef __debug_rangecheck
  if (!isinitialized()) {
    cout << "sim_data_exp_ttl::getorigvector() - data structure not initialized. call init." << endl;
    exit(EXIT_FAILURE);
  }
  #endif
  
  #ifdef __debug_rangecheck
  if (pos>=getnv()) {
    cout << "sim_data_exp_ttl::getorigvector() error: pos exceeds vector size " << pos << ">=" << getnv() << endl;
    exit(EXIT_FAILURE);
  }
  #endif
  
  return origvals[pos];
}

vector<double> sim_data_exp_fast::getorigvector(void) {
  
  #ifdef __debug_rangecheck
  if (!isinitialized()) {
    cout << "sim_data_exp_ttl::getorigvector() - data structure not initialized. call init." << endl;
    exit(EXIT_FAILURE);
  }
  #endif
  
  return origvals;
}
