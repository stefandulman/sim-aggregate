
/** 
  
  \file
  
 */

#ifndef __sim_algorithms_h
#define __sim_algorithms_h

// comp(T,T) - returns true if the first element is smaller than the second

template <class T>
bool alg_defaultcompare(const class& T x,  const class& T y) {
  if (x<y) return true;
  return false;
}

template <class C, class T, class Compare>
T alg_getmin(const C& container, Compare comp = alg_defaultcompare) {
  
  if (container.size()==0) {
    cout << "error: alg_getmin called with empty container." << endl;
    exit(EXIT_FAILURE);
  }
  
  C::iterator minval = container.begin(); 
  while (C::iterator it=container.begin(); it!=container.end(); ++it) {
    if (comp(*it, *minval) {
      minval = it;
    }
  }
  
  return *minval;
}


template <class C, class T, class Compare>
int alg_getminpos(const C& container, Compare comp = alg_defaultcompare) {
  
  if (container.size()==0) {
    cout << "error: alg_getminpos called with empty container." << endl;
    exit(EXIT_FAILURE);
  }
  
  int minpos = 0;
  int counter = 0;
  C::iterator minval = container.begin(); 
  while (C::iterator it=container.begin(); it!=container.end(); ++it) {
    if (comp(*it, *minval) {
      minval = it;
      minpos = counter;
    }
    counter++;
  }
  
  return counter;
}





#endif

