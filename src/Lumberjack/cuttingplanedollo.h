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
class CuttingPlaneDollo: public CuttingPlane
{
public:

  /// Constructor
  CuttingPlaneDollo(SATSolver* solver,
                    const Matrix& B,
                    StlIntMatrix& loss_vars,
                    StlIntMatrix& false_neg_vars,
                    StlIntMatrix& false_pos_vars);
  
protected:
  
  /// Get current assignment from solver and input
  int getEntryAssignment(int p, int c);
  
  /// Identify violated constraints
  int separate();
  
protected:

  /// Input matrix
  const Matrix& B_;
  /// Number of taxa
  const int m_;
  /// Number of characters
  const int n_;

  /// loss_vars maps matrix entries to loss variables
  StlIntMatrix& loss_vars_;
  /// false_neg_vars maps matrix entries to false neg variables
  StlIntMatrix& false_neg_vars_;
  /// false_pos_vars maps matrix entries to false pos variables
  StlIntMatrix& false_pos_vars_;
  
};

#endif // COLUMNGEN_H
