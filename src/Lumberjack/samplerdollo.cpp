#include "samplerdollo.h"
#include "cuttingplanedollo.h"
#include <map>

using std::vector;
using std::map;

SamplerDollo::SamplerDollo(const Matrix& B, size_t k, AppMC* appmc, UniG* unigen)
  : B_(B)
  , m_(B.getNrClones())
  , n_(B.getNrMutations())
  , k_(k)
  , num_vars_(0)
  , loss_vars_()
  , false_pos_vars_()
  , false_neg_vars_()
  , num_constraints_(0)
  , approxmc_(appmc)
  , unigen_(unigen)
{
}

void SamplerDollo::Init() {
  InitializeVariableMatrices();
  PrintVariableMatrices();
  cutting_plane_ = new CuttingPlaneDollo(approxmc_->get_solver(), B_, loss_vars_, false_neg_vars_, false_pos_vars_);

  // Update sampling set
  SATSolver* solver = approxmc_->get_solver();
  solver->new_vars(num_vars_);

  vector<uint32_t> sampling_set;
  for (int i = 0; i < num_vars_; ++i)
  {
    sampling_set.push_back(i);
  }
  approxmc_->set_sampling_set(sampling_set);

  std::cout << "Adding conflicting clauses\n";
  AddConflictingValuesClauses();
  approxmc_->get_solver()->dump_irred_clauses(&std::cout);
}

void SamplerDollo::InitializeVariableMatrices() {
  loss_vars_.resize(m_);
  false_pos_vars_.resize(m_);
  false_neg_vars_.resize(m_);

  for (size_t i = 0; i < m_; i++) {
    loss_vars_[i].resize(n_);
    false_pos_vars_[i].resize(n_);
    false_neg_vars_[i].resize(n_);

    for (size_t j = 0; j < n_; j++) {
      num_vars_++;
      loss_vars_[i][j] = num_vars_;
      
      if (B_.getEntry(i, j) == 1) {
        num_vars_++;
        false_pos_vars_[i][j] = num_vars_;
      } else {
        num_vars_++;
        false_neg_vars_[i][j] = num_vars_;
      }
    }
  }
}

void SamplerDollo::AddConflictingValuesClauses() {
  for (size_t i = 0; i < m_; i++) {
    for (size_t j = 0; j < n_; j++) {
      vector<Lit> clause;

      if (B_.getEntry(i, j) == 0) {
        clause.push_back(Lit(loss_vars_[i][j], false));
        std::cout << "Var " << loss_vars_[i][j] << " false\n";
        std::cout << "clause1 " << clause << std::endl;
        clause.push_back(Lit(false_neg_vars_[i][j], false));
        std::cout << "Var " << false_neg_vars_[i][j] << " false\n";
        std::cout << "clause2 " << clause << std::endl;
      } else {
        clause.push_back(Lit(loss_vars_[i][j], false));
        std::cout << "Var " << loss_vars_[i][j] << " false\n";
        std::cout << "clause1 " << clause << std::endl;
        clause.push_back(Lit(false_pos_vars_[i][j], true));
        std::cout << "Var " << false_pos_vars_[i][j] << " true\n";
        std::cout << "clause2 " << clause << std::endl;
      }

      approxmc_->add_clause(clause);
    }
  }
}

void SamplerDollo::Sample(const ApproxMC::SolCount *sol_count, uint32_t num_samples) {
  vector<vector<int>> solutions = unigen_->sample(sol_count, num_samples);

  std::cout << solutions.size() << " solutions sampled\n";

  for (auto solution : solutions) {
    std::cout << "Solution:\n";
    map<int, bool> sol_map = GetSolutionMap(solution);
    for (size_t i = 0; i < m_; i++) {
      for (size_t j = 0; j < n_; j++) {
        int entry = GetAssignmentFromSolution(sol_map, i, j);
        std::cout << entry << " ";
      }
      std::cout << "\n";
    }
  }
}

lbool SamplerDollo::GetAssignment(size_t var) {
  SATSolver* solver = approxmc_->get_solver();
  return solver->get_model()[var];
}

void SamplerDollo::PrintVariableMatrices() {

  std::cout << "Loss variables\n";

  for (size_t i = 0; i < m_; i++) {
    for (size_t j = 0; j < n_; j++) {
      std::cout << loss_vars_[i][j] << " ";
    }

    std::cout << "\n";
  }

  std::cout << "False negative variables\n";

  for (size_t i = 0; i < m_; i++) {
    for (size_t j = 0; j < n_; j++) {
      std::cout << false_neg_vars_[i][j] << " ";
    }

    std::cout << "\n";
  }

  std::cout << "False positive variables\n";

  for (size_t i = 0; i < m_; i++) {
    for (size_t j = 0; j < n_; j++) {
      std::cout << false_pos_vars_[i][j] << " ";
    }

    std::cout << "\n";
  }
}

map<int, bool> SamplerDollo::GetSolutionMap(const vector<int>& solution) {
  map<int, bool> sol_map;
  for (auto num : solution) {
    sol_map[abs(num)] = num > 0;
  }
  return sol_map;
}

int SamplerDollo::GetAssignmentFromSolution(map<int, bool>& solution, size_t clone, size_t mutation) {
  int loss_var = loss_vars_[clone][mutation];
  int false_neg_var = false_neg_vars_[clone][mutation];
  int false_pos_var = false_pos_vars_[clone][mutation];

  int entry = B_.getEntry(clone, mutation);

  if (solution[loss_var]) {
    if( (entry == 0 && !solution[false_neg_var]) || (entry == 1 && solution[false_pos_var]) ) {
      std::cout << "Something went wrong ";
      std::cout << (entry == 0 && !solution[false_neg_var]) << " ";
      std::cout << (entry == 1 && solution[false_pos_var]) << std::endl;
    }
    return 2;
  }

  if ((entry == 0 && !solution[false_neg_var]) || (entry == 1 && solution[false_pos_var])) {
    return 0;
  }

  return 1;
}