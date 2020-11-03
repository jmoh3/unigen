/*
 * cuttingplanedollo.cpp
 *
 */

#include "cuttingplanedollo.h"
#include <vector>
#include <iostream>
#include <string>

using namespace CMSat;
using std::vector;
using std::pair;

CuttingPlaneDollo::CuttingPlaneDollo(SATSolver* solver,
                    const Matrix& B,
                    StlIntMatrix& loss_vars,
                    StlIntMatrix& false_neg_vars,
                    StlIntMatrix& false_pos_vars)
  : CuttingPlane(solver)
  , B_(B)
  , m_(B.getNrClones())
  , n_(B.getNrMutations())
  , loss_vars_(loss_vars)
  , false_neg_vars_(false_neg_vars)
  , false_pos_vars_(false_pos_vars)
{
}

int CuttingPlaneDollo::getEntryAssignment(int p, int c) {
  int loss_var = loss_vars_[p][c];

  if (getAssignment(loss_var) == l_True) {
    return 2;
  }
  
  int original_val = B_.getEntry(p, c);

  if (original_val == 0) {
    int false_neg_var = false_neg_vars_[p][c];

    if (getAssignment(false_neg_var) == l_True) {
      return 1;
    } else {
      return 0;
    }
  } else {
    int false_pos_var = false_pos_vars_[p][c];

    if (getAssignment(false_pos_var) == l_True) {
      return 0;
    } else {
      return 1;
    }
  }

  throw std::runtime_error("Error: Solver did not assign truth value to variable.");
}

string CuttingPlaneDollo::getSubmatrixAsString(vector<pair<size_t, size_t>> positions) {
  vector<int> submatrix_assignments;

  for (auto position : positions) {
    submatrix_assignments.push_back(getEntryAssignment(position.first, position.second));
  }

  string submatrix_str;

  for (auto entry : submatrix_assignments) {
    submatrix_str.append(std::to_string(entry));
  }

  return submatrix_str;
}

vector<Lit> CuttingPlaneDollo::getLits(int p, int c) {
  int loss_var = loss_vars_[p][c];
  int false_neg_var = false_neg_vars_[p][c];
  int false_pos_var = false_pos_vars_[p][c];

  Lit loss_lit(loss_var, getAssignment(loss_var) == l_True);
  Lit false_neg_lit(false_neg_var, getAssignment(false_neg_var) == l_True);
  Lit false_pos_lit(false_pos_var, getAssignment(false_pos_var) == l_True);

  vector<Lit> lits { loss_lit, false_neg_lit, false_pos_lit };

  return lits;
}

int CuttingPlaneDollo::separate() {
  for (size_t row1 = 0; row1 < m_; row1++) {
    for (size_t row2 = 0; row2 < m_; row2++) {
      if (row1 == row2) {
        continue;
      }
      for (size_t row3 = 0; row3 < n_; row3++) {
        if (row3 == row2 || row3 == row1) {
          continue;
        }
        for (size_t col1 = 0; col1 < n_; col1++) {
          for (size_t col2 = 0; col2 < n_; col2++) {
            if (col1 == col2) {
              continue;
            }

            pair<size_t, size_t> b_11_pos(row1, col1);
            pair<size_t, size_t> b_12_pos(row1, col2);

            pair<size_t, size_t> b_21_pos(row2, col1);
            pair<size_t, size_t> b_22_pos(row2, col2);

            pair<size_t, size_t> b_31_pos(row3, col1);
            pair<size_t, size_t> b_32_pos(row3, col2);

            vector<pair<size_t, size_t>> positions {b_11_pos, b_12_pos, b_21_pos, b_22_pos, b_31_pos, b_32_pos};
            string submatrix_str = getSubmatrixAsString(positions);

            if (forbidden_submatrices_.find(submatrix_str) != forbidden_submatrices_.end()) {
              // submatrix is forbidden

              // get literals corresponding to each entry
              vector<Lit> clause;
              for (auto position : positions) {
                vector<Lit> entry_lits = getLits(position.first, position.second);
                for (size_t i = 0; i < entry_lits.size(); i++) {
                  // negate each one
                  entry_lits[i] = ~entry_lits[i];
                  clause.push_back(entry_lits[i]);
                }
              }
              
              addClause(clause);
            }
          }
        }
      }
    }
  }

  return 0;
}
