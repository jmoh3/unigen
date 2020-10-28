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
using ApproxMC::SolCount;

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
  virtual void Init();
  
  /// Solve
  int Sample(const SolCount *sol_count, uint32_t num_samples);
  
protected:
  /// Initializes variable matrices that define entries of corrected matrix
  void InitializeVariableMatrices();

  /// Get clauses that prevent conflicting values
  std::vector<std::vector<Lit>> GetConflictingValuesClauses();

  /// Get current assignment from solver and input
  int GetEntryAssignment(size_t clone, size_t mutation);

  /// Get current assignment from solver and input
  lbool GetAssignment(size_t var);
  
protected:

  /// Input matrix
  const Matrix& B_;
  /// Number of taxa
  const int m_;
  /// Number of characters
  const int n_;
  /// Maximum number of losses
  const int k_;
  
  /// loss_vars_[p][c] maps matrix entries to their loss variables
  StlIntMatrix loss_vars_;
  /// false_pos_vars_ maps matrix entries to false positive variables
  StlIntMatrix false_pos_vars_;
  /// false_neg_vars_ maps matrix entries to false positive variables
  StlIntMatrix false_neg_vars_;
  /// Number of variables
  int num_vars_;
  /// Number of constraints
  int num_constraints_;
  
  /// Approx MC solver
  AppMC* approxmc_;
  /// UniGen
  UniG* unigen_;
  /// Cutting plane oracle
  CuttingPlaneDollo* cutting_plane_;
};

#endif // COLUMNGEN_H
