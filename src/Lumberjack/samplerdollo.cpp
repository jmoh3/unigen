#include "samplerdollo.h"
#include "cuttingplanedollo.h"

using std::vector;

SamplerDollo::SamplerDollo(const Matrix& B, int k, AppMC* appmc, UniG* unigen)
  : B_(B)
  , m_(B.getNrClones())
  , n_(B.getNrMutations())
  , k_(k)
  , num_vars_(1)
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

  // Update sampling set
  vector<uint32_t> sampling_set;
  for (int i = 0; i < num_vars_; ++i)
  {
    sampling_set.push_back(i);
  }
  approxmc_->set_sampling_set(sampling_set);
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

vector<vector<Lit>> SamplerDollo::GetConflictingValuesClauses() {
  vector<vector<Lit>> clauses;

  for (size_t i = 0; i < m_; i++) {
    for (size_t j = 0; j < n_; j++) {
      vector<Lit> clause;

      if (B_.getEntry(i, j) == 1) {
        clause.push_back(Lit(loss_vars_[i][j], false));
        clause.push_back(Lit(false_neg_vars_[i][j], false));
      } else {
        clause.push_back(Lit(loss_vars_[i][j], false));
        clause.push_back(Lit(false_pos_vars_[i][j], true));
      }
    }
  }

  return clauses;
}

int SamplerDollo::Sample(const ApproxMC::SolCount *sol_count, uint32_t num_samples) {
  unigen_->sample(sol_count, num_samples);
  return 0;
}

lbool SamplerDollo::GetAssignment(size_t var) {
  SATSolver* solver = approxmc_->get_solver();
  return solver->get_model()[var];
}

int SamplerDollo::GetEntryAssignment(size_t clone, size_t mutation){
  return 0;
}
