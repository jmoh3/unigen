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
#include <vector>

using namespace UniGen;
using ApproxMC::AppMC;
using ApproxMC::SolCount;
using std::map;
using std::vector;

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
  SamplerDollo(const Matrix& B, size_t k, AppMC* appmc, UniG* unigen, 
     size_t cell_clusters, size_t mutation_clusters,
     double false_pos_rate=0.01, double false_neg_rate=0.5);
  
  /// Initializes solver
  virtual void Init();
  
  /// Samples solutions from current 1-Dollo instance
  /// @param sol_count
  /// @param num_samples desired number of samples
  void Sample(const SolCount *sol_count, uint32_t num_samples);
  
protected:

  /// Initializes variable matrices that define entries of corrected matrix
  void InitializeVariableMatrices();

  /// Add clauses that prevent conflicting values to CNF formula
  void AddConflictingValuesClauses();

  /// Adds clauses that enforce the values of the pair in column equal variables
  /// i.e. B[row1][col] == B[row2][col] => pair_in_col_equal[col][row1][row2]
  void AddColPairsEqualClauses();

  /// Adds clauses that enforce the values of the pair in row equal variables
  /// i.e. B[row][col1] == B[row][col2] => pair_in_row_equal[row][col1][col2]
  void AddRowPairsEqualClauses();

  void AddRowDuplicateClauses();

  void AddColDuplicateClauses();

  /// Get current assignment of a variable from solver and input
  /// @param var label for variable to get assignment for
  /// @return true or false
  lbool GetAssignment(size_t var);

  /// Gets a map representing truth assignments for a solution
  /// @param solution a vector of ints each entry is a variable, which is assigned 
  /// true if positive and false o/w
  /// @return a map of variable label to truth assignment
  map<int, bool> GetSolutionMap(const vector<int>& solution);

  /// Gets resulting solution matrix
  vector<vector<int>> GetSolMatrix(map<int, bool> sol_map);

  /// Checks to make sure all clustering variables are all correctly set.
  void ValidateSolution(map<int, bool> sol_map, vector<vector<int>> sol_matrix);

  /// Get current assignment of a variable from a solution map
  /// @param solution a map of variable label to truth assignment
  /// @param clone
  /// @param mutation
  /// @return 0, 1, or 2
  int GetAssignmentFromSolution(map<int, bool>& solution, size_t clone, size_t mutation);

  /// Helper method that gets an adder object for the current instance
  Adder GetAdder();

  /// Updates approxmc_'s sampling set with new variables
  void UpdateSamplingSet();
  
  /// Helper method that prints out all the variable matrices for an instance
  void PrintVariableMatrices();

  /// Gets a new unused variable label
  /// @return new variable
  int GetNewVar();

  /// Gets the variable label corresponding to entry at row, col being one
  /// @return variable label
  int GetEntryIsOneVar(size_t row, size_t col) const;

  /// Gets variable labels corresponding to entry at row, col being zero
  /// @return variable label
  vector<int> GetEntryIsZeroVars(size_t row, size_t col) const;

  /// Adds clauses to current formula
  /// @param clauses clauses to add
  void AddClauses(vector<vector<int>> clauses);

  /// Adds clause to current formula
  /// @param clause clause to add
  void AddClause(vector<int> clause);

  /// Adds a clause to imply lhs => rhs in formula
  /// @param lhs the literals on the left hand side of implication
  /// @param rhs the literal on the right hand side of the implication
  void AddImplyClause(vector<int> lhs, int rhs);
  
  /// Adds clauses to imply lhs => rhs in formula
  /// rhs.size() clauses will be added
  /// @param lhs the literals on the left hand side of implication
  /// @param rhs the literals on the right hand side of the implication
  void AddImplyClauses(vector<int> lhs, vector<int> rhs);

  /// Adds clauses to imply entry1 == entry2 => pair_equal_var
  void SetPairOfVarsEqual(int entry1, int entry2, int pair_equal_var);

  /// Adds clauses to imply entry1 == entry2 => pair_equal_var
  /// Overloaded function for when multiple variables correspond to one entry
  void SetPairOfVarsEqual(vector<int> entry1, vector<int> entry2, int pair_equal_var);

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

  const size_t num_cell_clusters_;
  const size_t num_mutation_clusters_;

  /// loss_vars_ maps matrix entries to their loss variables
  StlIntMatrix loss_vars_;
  /// false_pos_vars_ maps matrix entries to false positive variables
  StlIntMatrix false_pos_vars_;
  /// false_neg_vars_ maps matrix entries to false positive variables
  StlIntMatrix false_neg_vars_;

  /// pair_in_row_equal_[i][j][k] is true if the ith element of row j and row k are equal
  vector<vector<vector<int>>> pair_in_row_equal_;
  /// pair_in_col_equal_[i][j][k] is true if the kth element of col i and col j are equal
  vector<vector<vector<int>>> pair_in_col_equal_;

  /// row_is_duplicate_of_[i][j] is true if row i is equal to row j
  vector<vector<int>> row_is_duplicate_of_;
  /// col_is_duplicate_of_[i][j] is true if col i is equal to col j
  vector<vector<int>> col_is_duplicate_of_;

  /// row_is_duplicate_[i] is true if row i is a duplicate of some previous row 0,1,..,i-1
  vector<int> row_is_duplicate_;
  /// col_is_duplicate_[i] is true if col i is a duplicate of some previous col 0,1,..,i-1
  vector<int> col_is_duplicate_;
  
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
