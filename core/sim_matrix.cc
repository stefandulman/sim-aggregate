
#include "sim_matrix.h"

void sim_matrix::init_random(double start, double end, int col) {  
  uniform_real_distribution<double> distribution(start, end);
  for (int i=0; i<_nrows; ++i) {
    if (col==-1) {
      for (int j=0; j<_ncols; ++j) {
        matrix[i][j] = distribution(generator);
      }
    } else {
      matrix[i][col] = distribution(generator);
    }
  }
}

void sim_matrix::init_const(double val, int col) {
  for (int i=0; i<_nrows; ++i) {
    if (col==-1) {
      for (int j=0; j<_ncols; ++j) {
        matrix[i][j] = val;
      }
    } else {
      matrix[i][col] = val;
    }
  }
}

void sim_matrix::init_inc(double start, double step, int col) {
  for (int i=0; i<_nrows; ++i) {
    double s = start;
    if (col==-1) {
      for (int j=0; j<_ncols; ++j) {
        matrix[i][j] = s;
      }
    } else {
      matrix[i][col] = s;
    }
    s += step;
  }
}

void sim_matrix::init_mat(sim_matrix &m) {
  // transfer elements
  matrix.clear();
  matrix = m.data();
  // adjust sizes
  _nrows = m.nrows();
  _ncols = m.ncols();
  // topology pointer copy
  _top = m._top;
  // other adjustements...
  
}



sim_matrix sim_matrix::average(void) {
  // compute average - numericaly unstable!
  vector<double> avg (_ncols, 0);
  for (int i=0; i<_nrows; ++i) {
    for (int j=0; j<_ncols; ++j) {
      avg[j] += matrix[i][j];
    }
  }
  
  sim_matrix m(*this);
  for (int j=0; j<_ncols; ++j) {
    m.init_const(avg[j]/_nrows,j);
  }
  
  return m;
}

sim_matrix sim_matrix::sum(void) {
  // compute sum - numericaly unstable!
  vector<double> sum (_ncols, 0);
  for (int i=0; i<_nrows; ++i) {
    for (int j=0; j<_ncols; ++j) {
      sum[j] += matrix[i][j];
    }
  }
  
  sim_matrix m(*this);
  for (int j=0; j<_ncols; ++j) {
    m.init_const(sum[j], j);
  }
  
  return m;  
}

sim_matrix sim_matrix::min(void) {
  vector<double> min (_ncols, DBL_MAX);
  for (int i=0; i<_nrows; ++i) {
    for (int j=0; j<_ncols; ++j) {
      if (min[j]>matrix[i][j]) {
        min[j] = matrix[i][j];
      }
    }
  }
  sim_matrix m(*this);
  for (int j=0; j<_ncols; ++j) {
    m.init_const(min[j], j);
  }
  return m;
}

sim_matrix sim_matrix::max(void) {
  // compute average - numericaly unstable!
  vector<double> max (_ncols, -DBL_MAX);
  for (int i=0; i<_nrows; ++i) {
    for (int j=0; j<_ncols; ++j) {
      if (max[j]<matrix[i][j]) {
        max[j] = matrix[i][j];
      }
    }
  }
  sim_matrix m(*this);
  for (int j=0; j<_ncols; ++j) {
    m.init_const(max[j], j);
  }
  return m;
}



sim_matrix sim_matrix::operator/(sim_matrix& rhs) {
  if ((_ncols!=rhs.ncols()) || (_nrows!=rhs.nrows())) {
    cout << "operator /= size mismatch." << endl;
    exit(1);
  }
  sim_matrix m(*this);
  for (int i=0; i<_nrows; ++i) {
    for (int j=0; j<_ncols; ++j) {
      m[i][j] = (*this)[i][j] / rhs[i][j];
    }
  }
  return m;
}

sim_matrix sim_matrix::operator+(sim_matrix& rhs) {
  if ((_ncols!=rhs.ncols()) || (_nrows!=rhs.nrows())) {
    cout << "operator += size mismatch." << endl;
    exit(1);
  }
  sim_matrix m(*this);
  for (int i=0; i<_nrows; ++i) {
    for (int j=0; j<_ncols; ++j) {
      m[i][j] = (*this)[i][j] + rhs[i][j];
    }
  }
  return m;
}

sim_matrix sim_matrix::operator*(sim_matrix& rhs) {
  if ((_ncols!=rhs.ncols()) || (_nrows!=rhs.nrows())) {
    cout << "operator *= size mismatch." << endl;
    exit(1);
  }
  sim_matrix m(*this);
  for (int i=0; i<_nrows; ++i) {
    for (int j=0; j<_ncols; ++j) {
      m[i][j] = (*this)[i][j] * rhs[i][j];
    }
  }
  return m;
}

sim_matrix sim_matrix::operator-(sim_matrix& rhs) {
  if ((_ncols!=rhs.ncols()) || (_nrows!=rhs.nrows())) {
    cout << "operator -= size mismatch." << endl;
    exit(1);
  }
  sim_matrix m(*this);
  for (int i=0; i<_nrows; ++i) {
    for (int j=0; j<_ncols; ++j) {
      m[i][j] = (*this)[i][j] - rhs[i][j];
    }
  }
  return m;
}

