
#include "../../core/sim_header.h"

int main(void) {

  sim_networkdata x, temp;
  
  temp = x.sum() / x.count();
  double avg = temp.value();
  
  cout << "average: " << avg << endl;
  
  x = x.sum(2) / x.sum();
  cout << "min val: " << x.min().value() << endl;
  
  cout << "done!" << endl;

  return 0;
}

