/**
  
  \file
  
 */

#include "sim_header.h"

#ifndef __SIM_MATRX_H
#define __SIM_MATRX_H

// matrix types definition
typedef vector<vector <double>> Matrix;
typedef vector<double> Row;

/** matrix data storage file and basic matrix operations.
 * if a valid pointer to a network topology is provided, then simulation of networked behavior is possible.
 * if a null pointer is provided, then the needed operation is faked - results on complexity are bogus in this case.
 */ 
class sim_matrix {

    Matrix matrix; 

    int _nrows;
    int _ncols;

  public:

      sim_topology *_top;  // move it to private

      sim_matrix(int nrows, int ncols=1, sim_topology *top=NULL) : _nrows(nrows), _ncols(ncols), _top(top) {
          for (int i=0; i<_nrows; ++i) {
              Row row(_ncols);
              matrix.push_back(row);
          }
          // non-null topologies not supported. remove only after implementing all functionality
          //if (top!=NULL) {
          //  exit(1);
          //}
      };
      sim_matrix(sim_matrix &m) {
          init_mat(m);
      }
      ~sim_matrix() {};

      // accessors
      int nrows(void)      { return _nrows; };
      int ncols(void)      { return _ncols; };
      Matrix& data(void)   { return matrix; };

      // row and element access - to-do: not safe!
      vector<double>& operator[] (int i) {
          return matrix[i];
      }
      
      // print operator
      friend std::ostream& operator<<(std::ostream& str, const sim_matrix& m) {
          for (int i=0; i<m._nrows; ++i) {
              for (int j=0; j<m._ncols; ++j) {
                  str << m.matrix[i][j] << " ";
              }
              str << endl;
          }
          return str;
      };

      // math operators
      //sim_matrix& operator/=(sim_matrix& rhs);
      //sim_matrix& operator+=(sim_matrix& rhs);
      //sim_matrix& operator*=(sim_matrix& rhs);
      //sim_matrix& operator-=(sim_matrix& rhs);

      sim_matrix operator/(sim_matrix& rhs);
      sim_matrix operator+(sim_matrix& rhs);
      sim_matrix operator*(sim_matrix& rhs);
      sim_matrix operator-(sim_matrix& rhs);

      // initialization col=-1 means initialize all columns, otherwise only selected column
      void init_random(double start=0, double end=1, int col=-1);
      void init_const(double val, int col=-1);
      void init_inc(double start=0, double step=1, int col=-1);
      void init_mat(sim_matrix &m);

      // basic operations - they are always column-wise
      sim_matrix average(void);
      sim_matrix sum(void);
      sim_matrix min(void);
      sim_matrix max(void);
      
};

#endif
