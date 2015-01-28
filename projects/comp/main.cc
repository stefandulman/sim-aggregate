
#include "../../core/sim_header.h"
#include "sim.h"

int main(void) {

  ofstream f;

  f.open("./build/runhops_1000.txt");
  { sim_double s(1000, 1,  NORMAL, 0.05);  f << "1\t" << s.run_hops() << endl;  };
  { sim_double s(1000, 2,  NORMAL, 0.05);  f << "2\t" << s.run_hops() << endl;  };
  { sim_double s(1000, 5,  NORMAL, 0.05);  f << "5\t" << s.run_hops() << endl;  };
  { sim_double s(1000, 10, NORMAL, 0.05);  f << "10\t" << s.run_hops() << endl; };
  f.close();  

  return 0;
}

