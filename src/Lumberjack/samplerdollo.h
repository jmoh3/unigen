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
#include "adder.h"
#include <map>

using namespace UniGen;
using ApproxMC::AppMC;
using ApproxMC::SolCount;
using std::map;

/// This class provides a cutting plane wrapper for CryptoMiniSAT
/// This can be used to solve the the k-DP problem .
class SamplerDollo
{
public:
  
  /// Constructor
  /// @param B Input matrix
  /// @param k Maximum number of losses per character
  /// @param appmc pointer to ApproxMC object
  /// @param unigen pointer to UniGen object
  /// @param false_pos_rate rate at which false positives occur in SCS data
  /// @param false_neg_rate rate at which false negatives occur in SCS data
  SamplerDollo(const Matrix& B, size_t k, AppMC* appmc, UniG* unigen, double false_pos_rate=0.01, double false_neg_rate=0.5);
  
  /// Initializes solver
  virtual void Init();
  
  /// Samples solutions from current 1-Dollo instance
  /// @param sol_count
  /// @param num_samples desired number of samples
  void Sample(const SolCount *sol_count, uint32_t num_samples);
  
protected:

  /// Initializes variable matrices that define entries of corrected matrix
  void InitializeVariableMatrices();

  /// Get clauses that prevent conflicting values
  /// @return vector of clauses that prevent conflicting values
  void AddConflictingValuesClauses();

  /// Get current assignment of a variable from solver and input
  /// @param var label for variable to get assignment for
  /// @return true or false
  lbool GetAssignment(size_t var);

  /// Gets a map representing truth assignments for a solution
  /// @param solution a vector of ints each entry is a variable, which is assigned 
  /// true if positive and false o/w
  /// @return a map of variable label to truth assignment
  map<int, bool> GetSolutionMap(const vector<int>& solution);

  /// Get current assignment of a variable from a solution map
  /// @param solution a map of variable label to truth assignment
  /// @param clone
  /// @param mutation
  /// @return 0, 1, or 2
  int GetAssignmentFromSolution(map<int, bool>& solution, size_t clone, size_t mutation);

  /// Helper method that gets an adder object for the current instance
  Adder GetAdder();

  void UpdateSamplingSet();
  
  /// Helper method that prints out all the variable matrices for an instance
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

  /// False negative rate
  const double fn_rate_;
  /// False positive rate
  const double fp_rate_;

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
