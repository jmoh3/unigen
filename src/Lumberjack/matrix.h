/*
 * matrix.h
 *
 */

#ifndef MATRIX_H
#define MATRIX_H

#include "utils.h"
#include <cstring>
#include <fstream>
#include <list>
#include <random>

/// This class models a (k-Dollo) completion matrix
class Matrix
{
public:
  
  /// Constructor
  ///
  /// @param m Number of clones
  /// @param n Number of mutations
  Matrix(int m, int n);
  
  /// Default constructor
  Matrix();
  
  /// Construct matrix from file. Returns NULL if construction fails.
  ///
  /// @param filename Filename
  static Matrix* parse(const std::string& filename);
  
    
  /// Return number of clones
  int getNrClones() const
  {
    return _m;
  }
  
  /// Return number of mutations
  int getNrMutations() const
  {
    return _n;
  }
  
  /// Return maximum number of losses per character
  int getMaxNrLosses() const
  {
    return _k;
  }
  
  
  /// Return entry in matrix
  int getEntry(int p, int c) const
  {
    assert(0 <= p && p < _m);
    assert(0 <= c && c < _n);
    
    return _D[p][c];
  }
  
  /// Set entry in matrix
  void setEntry(int p, int c, int i)
  {
    assert(0 <= p && p < _m);
    assert(0 <= c && c < _n);

    _D[p][c] = i;
    
    if (i - 1 > _k)
    {
      _k = i - 1;
    }
  }
  

protected:
  
  /// Number of taxa
  int _m;
  /// Number of characters
  int _n;
  /// Input matrix
  StlIntMatrix _D;
  /// Number of losses
  int _k;

  friend std::ostream& operator<<(std::ostream& out, const Matrix& D);
  friend std::istream& operator>>(std::istream& in, Matrix& D);
};

/// Write matrix to output stream
///
/// @param out Output stream
/// @param D Matrix
std::ostream& operator<<(std::ostream& out, const Matrix& D);

/// Read matrix from input stream
///
/// @param in Input stream
/// @param D Matrix
std::istream& operator>>(std::istream& in, Matrix& D);

#endif // MATRIX_H
