#include "samplerdollo.h"
#include "cuttingplanedollo.h"
#include "adder.h"
#include <map>

using std::map;
using std::vector;

SamplerDollo::SamplerDollo(const Matrix &B, size_t k, AppMC *appmc, UniG *unigen,
                           size_t cell_clusters, size_t mutation_clusters, double false_pos_rate, double false_neg_rate)
    : B_(B), m_(B.getNrClones()), n_(B.getNrMutations()), k_(k), num_vars_(1), loss_vars_(), false_pos_vars_(), false_neg_vars_(), num_constraints_(0), approxmc_(appmc), unigen_(unigen), fp_rate_(false_pos_rate), fn_rate_(false_neg_rate), num_cell_clusters_(cell_clusters), num_mutation_clusters_(mutation_clusters)
{
}

void SamplerDollo::Init()
{
  InitializeVariableMatrices();
  PrintVariableMatrices();
  cutting_plane_ = new CuttingPlaneDollo(approxmc_->get_solver(), B_, loss_vars_, false_neg_vars_, false_pos_vars_);
  unigen_->set_cutting_plane(cutting_plane_);
  approxmc_->setCuttingPlane(cutting_plane_);

  Adder adder = GetAdder();
  num_vars_ += adder.GetNumVarsAdded();

  std::cout << "Updating sampling set\n";
  UpdateSamplingSet();

  vector<Lit> tmp{Lit(0, false)};
  approxmc_->add_clause(tmp);

  std::cout << "Adding conflicting clauses\n";
  AddConflictingValuesClauses();

  std::cout << "Adding clustering clauses\n";
  AddColPairsEqualClauses();
  AddRowPairsEqualClauses();

  std::cout << "Adding row duplicate clauses\n";
  AddRowDuplicateClauses();
  AddColDuplicateClauses();

  std::cout << "Adding fp and fn constraint clauses\n";
  vector<vector<Lit>> adder_clauses = adder.GetClauses();

  for (auto clause : adder_clauses)
  {
    approxmc_->add_clause(clause);
  }
}

void SamplerDollo::InitializeVariableMatrices()
{
  loss_vars_.resize(m_);
  false_pos_vars_.resize(m_);
  false_neg_vars_.resize(m_);

  for (size_t i = 0; i < m_; i++)
  {
    loss_vars_[i].resize(n_);
    false_pos_vars_[i].resize(n_);
    false_neg_vars_[i].resize(n_);

    for (size_t j = 0; j < n_; j++)
    {
      loss_vars_[i][j] = GetNewVar();

      if (B_.getEntry(i, j) == 1)
      {
        false_pos_vars_[i][j] = GetNewVar();
      }
      else
      {
        false_neg_vars_[i][j] = GetNewVar();
      }
    }
  }

  // begin clustering vars
  pair_in_row_equal_.resize(m_);
  for (size_t i = 0; i < m_; i++)
  {
    pair_in_row_equal_[i].resize(n_);

    for (size_t j = 0; j < n_; j++)
    {
      pair_in_row_equal_[i][j].resize(n_);

      for (size_t k = j + 1; k < n_; k++)
      {
        pair_in_row_equal_[i][j][k] = GetNewVar();
      }
    }
  }

  pair_in_col_equal_.resize(n_);
  for (size_t i = 0; i < n_; i++)
  {
    pair_in_col_equal_[i].resize(m_);

    for (size_t j = 0; j < m_; j++)
    {
      pair_in_col_equal_[i][j].resize(m_);

      for (size_t k = j + 1; k < m_; k++)
      {
        pair_in_col_equal_[i][j][k] = GetNewVar();
      }
    }
  }

  row_is_duplicate_of_.resize(m_);
  for (size_t i = 0; i < m_; i++)
  {
    row_is_duplicate_of_[i].resize(m_);
    for (size_t j = i + 1; j < m_; j++)
    {
      row_is_duplicate_of_[i][j] = GetNewVar();
    }
  }

  col_is_duplicate_of_.resize(n_);
  for (size_t i = 0; i < n_; i++)
  {
    col_is_duplicate_of_[i].resize(n_);
    for (size_t j = i + 1; j < n_; j++)
    {
      col_is_duplicate_of_[i][j] = GetNewVar();
    }
  }

  row_is_duplicate_.resize(m_);
  for (size_t i = 0; i < m_; i++)
  {
    row_is_duplicate_[i] = GetNewVar();
  }

  col_is_duplicate_.resize(n_);
  for (size_t i = 0; i < n_; i++)
  {
    col_is_duplicate_[i] = GetNewVar();
  }
}

Adder SamplerDollo::GetAdder()
{
  Adder adder(num_vars_);

  // false negative constraints
  vector<int> false_neg_flattened;
  for (auto row : false_neg_vars_)
  {
    for (auto var : row)
    {
      if (var != 0)
      {
        false_neg_flattened.push_back(var);
      }
    }
  }
  size_t max_fn = fn_rate_ * false_neg_flattened.size();
  std::cout << "Max num false negatives: " << max_fn << std::endl;
  adder.EncodeEqualToK(false_neg_flattened, max_fn);

  // false positive constraints
  vector<int> false_pos_flattened;
  for (auto row : false_pos_vars_)
  {
    for (auto var : row)
    {
      if (var != 0)
      {
        false_pos_flattened.push_back(var);
      }
    }
  }
  size_t max_fp = fp_rate_ * false_pos_flattened.size();
  std::cout << "Max num false positives: " << max_fp << std::endl;
  adder.EncodeEqualToK(false_pos_flattened, max_fp);

  // Row clustering constraints
  size_t num_row_duplicates = m_ - num_cell_clusters_;
  std::cout << "Num row duplicates: " << num_row_duplicates << std::endl;
  adder.EncodeEqualToK(row_is_duplicate_, num_row_duplicates);

  // Column clustering constraints
  size_t num_col_duplicates = n_ - num_mutation_clusters_;
  std::cout << "Num col duplicates: " << num_col_duplicates << std::endl;
  adder.EncodeEqualToK(col_is_duplicate_, num_col_duplicates);

  return adder;
}

void SamplerDollo::AddConflictingValuesClauses()
{
  for (size_t i = 0; i < m_; i++)
  {
    for (size_t j = 0; j < n_; j++)
    {
      vector<int> clause;

      if (B_.getEntry(i, j) == 0)
      {
        clause.push_back(-loss_vars_[i][j]);
        clause.push_back(-false_neg_vars_[i][j]);
      }
      else
      {
        clause.push_back(-loss_vars_[i][j]);
        clause.push_back(false_pos_vars_[i][j]);
      }

      AddClause(clause);
    }
  }
}

void SamplerDollo::AddRowDuplicateClauses()
{
  for (size_t row2 = 1; row2 < m_; row2++)
  {
    std::cout << "LARGER ROW: " << row2 << std::endl;
    // row_is_duplicate[row] =>
    // row_is_duplicate_of[0][row] or ... row_is_duplicate_of[row-1][row]
    vector<int> clause_only_if{-row_is_duplicate_[row2]};
    for (size_t row1 = 0; row1 < row2; row1++)
    {
      std::cout << "SMALLER ROW: " << row1 << std::endl;
      // if all pairs in column are equal, then row is duplicate of prev row
      // pair_in_col_equal[0][smaller][row] and ... pair_in_col_equal[n][smaller][row]
      // => row_is_duplicate_of[smaller][row]
      vector<int> clause_if;

      for (size_t col = 0; col < n_; col++)
      {
        std::cout << "COL: " << col << std::endl;
        clause_if.push_back(-pair_in_col_equal_[col][row1][row2]);

        // row_is_duplicate_of[smaller][row] => pair_in_col_equal[col][smaller][row]
        vector<int> pair_in_col_clause{-row_is_duplicate_of_[row1][row2], pair_in_col_equal_[col][row1][row2]};
        AddClause(pair_in_col_clause);
      }

      clause_if.push_back(row_is_duplicate_of_[row1][row2]);
      AddClause(clause_if);

      // row_is_duplicate_of[smaller][row] => row_is_duplicate[row]
      vector<int> row_is_duplicate_clause{-row_is_duplicate_of_[row1][row2], row_is_duplicate_[row2]};
      AddClause(row_is_duplicate_clause);

      clause_only_if.push_back(row_is_duplicate_of_[row1][row2]);
    }

    AddClause(clause_only_if);
  }

  // first row cannot be a duplicate
  vector<int> first_row_clause{-row_is_duplicate_[0]};
  AddClause(first_row_clause);
}

void SamplerDollo::AddColDuplicateClauses()
{
  for (size_t col2 = 1; col2 < n_; col2++)
  {
    // col_is_duplicate[col] =>
    // col_is_duplicate_of[0][col] or ... col_is_duplicate_of[col-1][col]
    vector<int> clause_only_if{-col_is_duplicate_[col2]};
    for (size_t col1 = 0; col1 < col2; col1++)
    {
      // pair_in_row_equal[0][smaller_col][col] and ... pair_in_row_equal[n][smaller_col][col]
      // => col_is_duplicate_of[smaller_col][col]
      vector<int> clause_if;

      for (size_t row = 0; row < m_; row++)
      {
        clause_if.push_back(-pair_in_row_equal_[row][col1][col2]);

        // col_is_duplicate_of[smaller_col][col] => pair_in_row_equal[row][smaller_col][col]
        vector<int> pair_in_col_clause{-col_is_duplicate_of_[col1][col2], pair_in_row_equal_[row][col1][col2]};
        AddClause(pair_in_col_clause);
      }

      clause_if.push_back(col_is_duplicate_of_[col1][col2]);
      AddClause(clause_if);

      // col_is_duplicate_of[smaller][col] => col_is_duplicate[col]
      vector<int> row_is_duplicate_clause{-col_is_duplicate_of_[col1][col2], col_is_duplicate_[col2]};
      AddClause(row_is_duplicate_clause);

      clause_only_if.push_back(col_is_duplicate_of_[col1][col2]);
    }

    AddClause(clause_only_if);
  }

  // first col cannot be a duplicate
  vector<int> first_col_clause{-col_is_duplicate_[0]};
  AddClause(first_col_clause);
}

void SamplerDollo::AddRowPairsEqualClauses()
{
  for (size_t row = 0; row < m_; row++)
  {
    for (size_t col1 = 0; col1 < n_; col1++)
    {
      for (size_t col2 = col1 + 1; col2 < n_; col2++)
      {
        int pair_in_row_equal_var = pair_in_row_equal_[row][col1][col2];

        /// BOTH ENTRIES ARE 1
        int rowcol1_is_one = GetEntryIsOneVar(row, col1);
        int rowcol2_is_one = GetEntryIsOneVar(row, col2);

        SetPairOfVarsEqual(rowcol1_is_one, rowcol2_is_one, pair_in_row_equal_var);

        /// BOTH ENTRIES ARE 2
        int rowcol1_is_two = loss_vars_[row][col1];
        int rowcol2_is_two = loss_vars_[row][col2];

        SetPairOfVarsEqual(rowcol1_is_two, rowcol2_is_two, pair_in_row_equal_var);

        // BOTH ENTRIES ARE 0
        vector<int> rowcol1_is_zero = GetEntryIsZeroVars(row, col1);
        vector<int> rowcol2_is_zero = GetEntryIsZeroVars(row, col2);

        SetPairOfVarsEqual(rowcol1_is_zero, rowcol2_is_zero, pair_in_row_equal_var);
      }
    }
  }
}

void SamplerDollo::AddColPairsEqualClauses()
{
  for (size_t col = 0; col < n_; col++)
  {
    for (size_t row1 = 0; row1 < m_; row1++)
    {
      for (size_t row2 = row1 + 1; row2 < m_; row2++)
      {
        int pair_in_col_equal_var = pair_in_col_equal_[col][row1][row2];

        /// BOTH ENTRIES ARE 1
        int row1col_is_one = GetEntryIsOneVar(row1, col);
        int row2col_is_one = GetEntryIsOneVar(row2, col);

        SetPairOfVarsEqual(row1col_is_one, row2col_is_one, pair_in_col_equal_var);

        /// BOTH ENTRIES ARE 2
        int row1col_is_two = loss_vars_[row1][col];
        int row2col_is_two = loss_vars_[row2][col];

        SetPairOfVarsEqual(row1col_is_two, row2col_is_two, pair_in_col_equal_var);

        // BOTH ENTRIES ARE 0
        vector<int> row1col_is_zero = GetEntryIsZeroVars(row1, col);
        vector<int> row2col_is_zero = GetEntryIsZeroVars(row2, col);

        SetPairOfVarsEqual(row1col_is_zero, row2col_is_zero, pair_in_col_equal_var);
      }
    }
  }
}

void SamplerDollo::UpdateSamplingSet()
{
  // Update sampling set
  std::cout << num_vars_ << " vars created total\n";
  SATSolver *solver = approxmc_->get_solver();
  solver->new_vars(num_vars_);

  vector<uint32_t> sampling_set;
  for (int i = 0; i < num_vars_; ++i)
  {
    sampling_set.push_back(i);
  }
  approxmc_->set_sampling_set(sampling_set);
}

void SamplerDollo::Sample(const ApproxMC::SolCount *sol_count, uint32_t num_samples)
{
  vector<vector<int>> solutions = unigen_->sample(sol_count, num_samples);

  std::cout << solutions.size() << " solutions sampled\n";

  for (auto solution : solutions)
  {
    std::cout << "===================\n";
    map<int, bool> sol_map = GetSolutionMap(solution);
    vector<vector<int>> sol_matrix = GetSolMatrix(sol_map);

    ValidateSolution(sol_map, sol_matrix);

    size_t num_fn = 0;
    size_t num_fp = 0;

    for (size_t i = 0; i < m_; i++)
    {
      for (size_t j = 0; j < n_; j++)
      {
        std::cout << sol_matrix[i][j] << " ";

        if (sol_matrix[i][j] == 1 && B_.getEntry(i, j) == 0)
        {
          num_fn++;
        }
        if ((sol_matrix[i][j] == 2 || sol_matrix[i][j] == 0) && B_.getEntry(i, j) == 1)
        {
          num_fp++;
        }
      }
      std::cout << "\n";
    }

    std::cout << "# false negatives: " << num_fn << "\n# false positives: " << num_fp << "\n";
  }
}

void SamplerDollo::ValidateSolution(map<int, bool> sol_map, vector<vector<int>> sol_matrix)
{
  // Verifies values of pair in row equal/pair in column equal variables
  for (size_t i = 0; i < m_; i++)
  {
    for (size_t j = 0; j < n_; j++)
    {
      for (size_t k = j + 1; k < n_; k++)
      {
        bool pair_equal = sol_map[pair_in_row_equal_[i][j][k]];
        if (pair_equal)
        {
          assert(sol_matrix[i][j] == sol_matrix[i][k]);
        }
        else
        {
          assert(sol_matrix[i][j] != sol_matrix[i][k]);
        }
      }
    }
  }

  for (size_t i = 0; i < n_; i++)
  {
    for (size_t j = 0; j < m_; j++)
    {
      for (size_t k = j + 1; k < m_; k++)
      {
        bool pair_equal = sol_map[pair_in_col_equal_[i][j][k]];
        if (pair_equal)
        {
          assert(sol_matrix[j][i] == sol_matrix[k][i]);
        }
        else
        {
          assert(sol_matrix[j][i] != sol_matrix[k][i]);
        }
      }
    }
  }

  // Verifies row is duplicate of and col is duplicate of variables

  for (size_t j = 1; j < m_; j++)
  {
    bool row_j_is_duplicate_actual = false;

    for (size_t i = 0; i < j; i++)
    {
      bool row_i_duplicate_of_j_expected = sol_map[row_is_duplicate_of_[i][j]];

      bool row_i_duplicate_of_j_actual = true;
      for (size_t k = 0; k < n_; k++)
      {
        if (sol_matrix[i][k] != sol_matrix[j][k])
        {
          row_i_duplicate_of_j_actual = false;
          break;
        }
      }
      if (row_i_duplicate_of_j_actual)
      {
        row_j_is_duplicate_actual = true;
      }

      assert(row_i_duplicate_of_j_expected == row_i_duplicate_of_j_actual);
    }
    bool row_j_is_duplicate_expected = sol_map[row_is_duplicate_[j]];
    assert(row_j_is_duplicate_actual == row_j_is_duplicate_expected);
  }

  for (size_t j = 1; j < n_; j++)
  {
    bool col_j_is_duplicate_actual = false;

    for (size_t i = 0; i < j; i++)
    {
      bool col_i_duplicate_of_j_expected = sol_map[col_is_duplicate_of_[i][j]];

      bool col_i_duplicate_of_j_actual = true;
      for (size_t k = 0; k < m_; k++)
      {
        if (sol_matrix[k][i] != sol_matrix[k][j])
        {
          col_i_duplicate_of_j_actual = false;
          break;
        }
      }

      if (col_i_duplicate_of_j_actual)
      {
        col_j_is_duplicate_actual = true;
      }

      assert(col_i_duplicate_of_j_expected == col_i_duplicate_of_j_actual);
    }

    bool col_j_is_duplicate_expected = sol_map[col_is_duplicate_[j]];
    assert(col_j_is_duplicate_actual == col_j_is_duplicate_expected);
  }
}

lbool SamplerDollo::GetAssignment(size_t var)
{
  SATSolver *solver = approxmc_->get_solver();
  return solver->get_model()[var];
}

vector<vector<int>> SamplerDollo::GetSolMatrix(map<int, bool> sol_map)
{
  vector<vector<int>> sol_matrix;

  for (size_t i = 0; i < m_; i++)
  {
    vector<int> matrix_row;
    for (size_t j = 0; j < n_; j++)
    {
      int entry = GetAssignmentFromSolution(sol_map, i, j);
      matrix_row.push_back(entry);
    }
    sol_matrix.push_back(matrix_row);
  }

  return sol_matrix;
}

void SamplerDollo::PrintVariableMatrices()
{
  std::cout << "Loss variables\n";

  for (size_t i = 0; i < m_; i++)
  {
    for (size_t j = 0; j < n_; j++)
    {
      std::cout << loss_vars_[i][j] + 1 << " ";
    }

    std::cout << "\n";
  }

  std::cout << "False negative variables\n";

  for (size_t i = 0; i < m_; i++)
  {
    for (size_t j = 0; j < n_; j++)
    {
      if (B_.getEntry(i, j) == 0)
      {
        std::cout << false_neg_vars_[i][j] + 1 << " ";
      }
      else
      {
        std::cout << "0 ";
      }
    }

    std::cout << "\n";
  }

  std::cout << "False positive variables\n";

  for (size_t i = 0; i < m_; i++)
  {
    for (size_t j = 0; j < n_; j++)
    {
      if (B_.getEntry(i, j) == 1)
      {
        std::cout << false_pos_vars_[i][j] + 1 << " ";
      }
      else
      {
        std::cout << "0 ";
      }
    }

    std::cout << "\n";
  }

  std::cout << "Pair in row equal" << std::endl;

  for (size_t i = 0; i < m_; i++)
  {
    std::cout << "Row " << i << std::endl;
    for (size_t j = 0; j < n_; j++)
    {
      for (size_t k = j + 1; k < n_; k++)
      {
        std::cout << "B_[" << i << "][" << j << "] = B_[" << i << "][" << k << "] <=> var " << pair_in_row_equal_[i][j][k] + 1 << std::endl;
      }
    }
  }

  std::cout << "Pair in column equal" << std::endl;

  for (size_t i = 0; i < n_; i++)
  {
    std::cout << "Column " << i << std::endl;
    for (size_t j = 0; j < m_; j++)
    {
      for (size_t k = j + 1; k < m_; k++)
      {
        std::cout << "B_[" << j << "][" << i << "] = B_[" << k << "][" << i << "] <=> var " << pair_in_col_equal_[i][j][k] + 1 << std::endl;
      }
    }
  }

  std::cout << "Row is duplicate of " << std::endl;

  for (size_t i = 0; i < m_; i++)
  {
    for (size_t j = i + 1; j < m_; j++)
    {
      std::cout << "B[" << i << "] == B[" << j << "] <=> var " << row_is_duplicate_of_[i][j] + 1 << std::endl;
    }
  }

  std::cout << "Col is duplicate of " << std::endl;

  for (size_t i = 0; i < n_; i++)
  {
    for (size_t j = i + 1; j < n_; j++)
    {
      std::cout << "B[:," << i << "] == B[:," << j << "] <=> var " << col_is_duplicate_of_[i][j] + 1 << std::endl;
    }
  }

  std::cout << "Row is duplicate" << std::endl;

  for (size_t i = 0; i < m_; i++)
  {
    std::cout << "Row " << i << " is duplicate <=> var " << row_is_duplicate_[i] + 1 << std::endl;
  }

  std::cout << "Col is duplicate" << std::endl;

  for (size_t i = 0; i < n_; i++)
  {
    std::cout << "Col " << i << " is duplicate <=> var " << col_is_duplicate_[i] + 1 << std::endl;
  }
}

map<int, bool> SamplerDollo::GetSolutionMap(const vector<int> &solution)
{
  map<int, bool> sol_map;
  for (auto num : solution)
  {
    sol_map[abs(num) - 1] = num > 0;
  }
  return sol_map;
}

int SamplerDollo::GetAssignmentFromSolution(map<int, bool> &solution, size_t clone, size_t mutation)
{
  int loss_var = loss_vars_[clone][mutation];
  int false_neg_var = false_neg_vars_[clone][mutation];
  int false_pos_var = false_pos_vars_[clone][mutation];

  int entry = B_.getEntry(clone, mutation);

  if (solution[loss_var])
  {
    assert((entry == 0 && !solution[false_neg_var]) || (entry == 1 && solution[false_pos_var]));
    return 2;
  }

  if ((entry == 0 && !solution[false_neg_var]) || (entry == 1 && solution[false_pos_var]))
  {
    return 0;
  }

  return 1;
}

int SamplerDollo::GetNewVar()
{
  int new_var = num_vars_;
  num_vars_++;
  return new_var;
}

int SamplerDollo::GetEntryIsOneVar(size_t row, size_t col) const
{
  if (B_.getEntry(row, col) == 0)
  {
    return false_neg_vars_[row][col];
  }
  return -false_pos_vars_[row][col];
}

vector<int> SamplerDollo::GetEntryIsZeroVars(size_t row, size_t col) const
{
  vector<int> zero_vars;
  zero_vars.push_back(-loss_vars_[row][col]);

  if (B_.getEntry(row, col) == 0)
  {
    zero_vars.push_back(-false_neg_vars_[row][col]);
  }
  else
  {
    zero_vars.push_back(false_pos_vars_[row][col]);
  }

  return zero_vars;
}

void SamplerDollo::AddClause(vector<int> clause)
{
  vector<Lit> lits;

  for (auto var : clause)
  {
    int label = abs(var);
    bool is_inverted = var < 0;

    Lit lit(label, is_inverted);
    lits.push_back(lit);
  }

  std::cout << "Adding clause " << lits << std::endl;
  approxmc_->add_clause(lits);
}

void SamplerDollo::AddClauses(vector<vector<int>> clauses)
{
  for (auto clause : clauses)
  {
    AddClause(clause);
  }
}

void SamplerDollo::AddImplyClause(vector<int> lhs, int rhs)
{
  vector<int> clause;

  for (auto var : lhs)
  {
    clause.push_back(-var);
  }

  clause.push_back(rhs);

  AddClause(clause);
}

void SamplerDollo::AddImplyClauses(vector<int> lhs, vector<int> rhs)
{
  for (auto var : rhs)
  {
    AddImplyClause(lhs, var);
  }
}

void SamplerDollo::SetPairOfVarsEqual(int entry1, int entry2, int pair_equal_var)
{
  // (entry1 == value) and (entry2 == value) => pair_equal
  vector<int> lhs{entry1, entry2};
  AddImplyClause(lhs, pair_equal_var);

  // pair_equal and (entry1 == value) => (entry2 == value)
  lhs = vector<int>{entry1, pair_equal_var};
  AddImplyClause(lhs, entry2);

  // pair_equal and (entry2 == value) => (entry1 == value)
  lhs = vector<int>{entry2, pair_equal_var};
  AddImplyClause(lhs, entry1);
}

void SamplerDollo::SetPairOfVarsEqual(vector<int> entry1, vector<int> entry2, int pair_equal_var)
{
  // entry1 == 1 and entry2 == 1 => pair_equal
  vector<int> lhs;
  lhs.insert(lhs.begin(), entry1.begin(), entry1.end());
  lhs.insert(lhs.begin(), entry2.begin(), entry2.end());
  AddImplyClause(lhs, pair_equal_var);

  // pair_equal and entry1 == 1 => entry2 == 1
  lhs = vector<int>{pair_equal_var};
  lhs.insert(lhs.begin(), entry1.begin(), entry1.end());
  AddImplyClauses(lhs, entry2);

  // pair_equal and entry2 == 1 => entry1 == 1
  lhs = vector<int>{pair_equal_var};
  lhs.insert(lhs.begin(), entry2.begin(), entry2.end());
  AddImplyClauses(lhs, entry1);
}