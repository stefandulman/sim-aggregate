
#include "sim_data_double.h"

sim_data_double::sim_data_double() {
}


void sim_data_double::init(double _val) {
  val = _val;
  
  #ifdef __debug_rangecheck
  initialize();
  #endif
}


bool sim_data_double::mix(sim_data *d2, sim_data *res1, sim_data *res2) {
  
  #ifdef __debug_rangecheck
  if (!isinitialized()) {
    cout << "sim_data_double::mix() - data structure not initialized. call init." << endl;
    exit(EXIT_FAILURE);
  }
  #endif
  
  double temp = (this->getvalue() + d2->getvalue()) / 2.0;
  
  bool res = false;
  if (this->getvalue() != d2->getvalue()) {
    res = true;
  }
  
  res1->setvalue(temp);
  res2->setvalue(temp);
  
  return res;
}


double sim_data_double::getvalue(void) {
  
  #ifdef __debug_rangecheck
  if (!isinitialized()) {
    cout << "sim_data_double::mix() - data structure not initialized. call init." << endl;
    exit(EXIT_FAILURE);
  }
  #endif
  
  return val;
}


void sim_data_double::setvalue(double _val) {
  
  #ifdef __debug_rangecheck
  if (!isinitialized()) {
    cout << "sim_data_double::mix() - data structure not initialized. call init." << endl;
    exit(EXIT_FAILURE);
  }
  #endif
  
  val = _val;
}

