/*
 * cuttingplane.h
 *
 */

#ifndef SAMPLERDOLLO_H
#define SAMPLERDOLLO_H

#include <cryptominisat5/cryptominisat.h>
#include "matrix.h"
#include "utils.h"
#include <approxmc/approxmc.h>
#include "cuttingplanedollo.h"
#include "unigen/unigen.h"

using namespace UniGen;
using ApproxMC::AppMC;

/// This class provides a cutting plane wrapper for CryptoMiniSAT
/// This can be used to solve the the k-DP problem .
class SamplerDollo
{
public:
    
  /// Constructor
  ///
  /// @param B Input matrix
  /// @param k Maximum number of losses per character
  SamplerDollo(const Matrix& B, int k, AppMC* appmc, UniG* unigen);
  
  /// Initialize solver
  virtual void init();
  
  /// Return solution matrix
  const Matrix& getSolA() const
  {
    return _solA;
  }
  
  /// Solve
  int solve(const ApproxMC::SolCount *sol_count, uint32_t num_samples);
  
protected:

  /// Get current assignment from solver and input
  int getEntryAssignment(int p, int c);

  /// Get current assignment from solver and input
  lbool getAssignment(int var);

  /// Extract solution from SAT solver
  void processSolution();
  
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
  StlIntMatrix _B2Var;
  /// _var2B maps active variable index to triple indexing marix B entry
  std::vector<Triple> _var2B;
  /// Indicates which matrix entries are active
  StlBoolMatrix _activeEntries;
  /// Number of active variables
  int _nrActiveVariables;
  /// Number of constraints
  int _nrConstraints;
  /// Solution matrix
  Matrix _solA;
  /// Approx MC solver
  AppMC* _approxmc;
  /// UniGen
  UniG* _unigen;
  /// Cutting plane oracle
  CuttingPlaneDollo* _cuttingPlane;
};

#endif // COLUMNGEN_H
