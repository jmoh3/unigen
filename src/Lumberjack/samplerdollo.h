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

  /*
    Constructor
    @param B Input matrix
    @param k Maximum number of losses per character
  */
  SamplerDollo(const Matrix& B, size_t k, AppMC* appmc, UniG* unigen);
  
  /*
    Initializes solver
  */
  virtual void Init();
  
  /*
    Samples solutions from current 1-Dollo instance
    @param sol_count
    @param num_samples desired number of samples
  */
  void Sample(const SolCount *sol_count, uint32_t num_samples);
  
protected:
  /*
    Initializes variable matrices that define entries of corrected matrix
  */
  void InitializeVariableMatrices();

  /*
    Get clauses that prevent conflicting values
    @return vector of clauses that prevent conflicting values
  */
  void AddConflictingValuesClauses();

  /*
    Get current assignment for a mutation in a clone from solver and input
    @param clone index of clone to get assignment for
    @param mutation index of mutation to get assignment for
    @return 0 (mutation not present), 1 (mutation present), or 2 (mutation lost)
  */
  int GetEntryAssignment(size_t clone, size_t mutation);

  /*
    Get current assignment of a variable from solver and input
    @param var label for variable to get assignment for
    @return true or false
  */
  lbool GetAssignment(size_t var);
  
  void PrintVariableMatrices();

protected:

  /// Input matrix
  const Matrix& B_;
  /// Number of taxa
  const size_t m_;
  /// Number of characters
  const size_t n_;
  /// Maximum number of losses
  const size_t k_;
  
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