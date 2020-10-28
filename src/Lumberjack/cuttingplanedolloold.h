/*
 * cuttingplane.h
 *
 */

#ifndef CUTTINGPLANEDOLLO_H
#define CUTTINGPLANEDOLLO_H

#include <cryptominisat5/cryptominisat.h>
#include "matrix.h"
#include "utils.h"
#include <approxmc/cuttingplane.h>

using ApproxMC::CuttingPlane;

/// This class provides a cutting plane wrapper for CryptoMiniSAT
/// This can be used to solve the the k-DP problem .
class CuttingPlaneDolloOld: public CuttingPlane
{
protected:
  /// Triple (p,c,i)
  struct Triple
  {
  public:
    Triple(int p, int c, int i)
      : _p(p)
      , _c(c)
      , _i(i)
    {
    }
    
    Triple()
      : _p(-1)
      , _c(-1)
      , _i(-1)
    {
    }
    
    int _p;
    int _c;
    int _i;
  };

public:
    
  /// Constructor
  CuttingPlaneDolloOld(SATSolver* solver, const Matrix& B, const int m, const int n, const int k, StlIntMatrix& B2Var, StlBoolMatrix& activeEntries, Matrix& solA);
  
  /// Return solution matrix
  const Matrix& getSolA() const
  {
    return _solA;
  }
  
protected:
  
  /// Get current assignment from solver and input
  int getEntryAssignment(int p, int c);
  
  /// Identify violated constraints
  int separate();
  
  /// Extract solution from SAT solver
  void processSolution();
      
  /// Forbidden submatrix
  typedef std::array<Triple, 6> ViolatedConstraint;
  
  /// List of forbidden submatrices
  typedef std::list<ViolatedConstraint> ViolatedConstraintList;
  
protected:
  /// Input matrix
  const Matrix& _B;
  /// Number of taxa
  const int _m;
  /// Number of characters
  const int _n;
  /// Maximum number of losses
  const int _k;
  /// _B2Var[p][c] maps matrix entries to active variable index
  StlIntMatrix& _B2Var;
  /// Indicates which matrix entries are active
  StlBoolMatrix& _activeEntries;
  /// Solution matrix
  Matrix& _solA;
};

#endif // COLUMNGEN_H
