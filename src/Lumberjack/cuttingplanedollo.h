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
#include <utility>
#include <vector>
#include <unordered_set>

using ApproxMC::CuttingPlane;
using std::vector;
using std::string;
using std::unordered_set;
using std::pair;

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

  /// Get literals for loss, false negative, and false positive variables corresponding to entry p, c
  vector<Lit> getLits(int p, int c);
  
  /// Identify violated constraints
  int separate();

  string getSubmatrixAsString(vector<pair<size_t, size_t>> positions);
  
protected:

  /// Input matrix
  const Matrix& B_;
  /// Number of taxa
  const int m_;
  /// Number of characters
  const int n_;
  /// Number of allowed losses
  const int k_ = 1;

  const unordered_set<string> forbidden_submatrices_ {"100111", "100112", "100211", "100212", "100121",
                                            "100122", "100221", "100222", "200111", "200112",
                                            "200211", "200212", "200121", "200122", "200221",
                                            "200222", "110212", "110222", "210212", "210222",
                                            "201121", "201122", "201221", "201222", "211222"};

  /// loss_vars maps matrix entries to loss variables
  StlIntMatrix& loss_vars_;
  /// false_neg_vars maps matrix entries to false neg variables
  StlIntMatrix& false_neg_vars_;
  /// false_pos_vars maps matrix entries to false pos variables
  StlIntMatrix& false_pos_vars_;
  
};

#endif // COLUMNGEN_H
