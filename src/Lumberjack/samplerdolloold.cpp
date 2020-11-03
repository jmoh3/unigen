#include "samplerdolloold.h"
#include "cuttingplanedolloold.h"

SamplerDolloOld::SamplerDolloOld(const Matrix& B, int k, AppMC* appmc, UniG* unigen)
  : _B(B)
  , _m(B.getNrClones())
  , _n(B.getNrMutations())
  , _k(k)
  , _B2Var()
  , _var2B()
  , _activeEntries()
  , _nrActiveVariables(0)
  , _nrConstraints(0)
  , _solA(_B.getNrClones(), _B.getNrMutations())
  , _approxmc(appmc)
  , _unigen(unigen)
{
}

void SamplerDolloOld::init() {
  
  /// Resize containers for each clone
  _activeEntries.resize(_m);
  _B2Var.resize(_m);
  
  for (int p = 0; p < _m; p++)
  {
    /// Resize container for each mutation
    _activeEntries[p].resize(_n);
    _B2Var[p].resize(_n);

    for (int c = 0; c < _n; ++c)
    {
      _activeEntries[p][c] = false;
      
      /// Check if this entry is elligible for being active
      if (_B.getEntry(p, c) == 0){
        /// Note that is it active and generate variables
        _activeEntries[p][c] = true;
        _B2Var[p][c] = _nrActiveVariables;
        _nrActiveVariables++;
      }
    }
  }

  // _falsePositiveVars.resize(_m);
  // _falseNegativeVars.resize(_m);

  // std::vector<std::vector<Lit>> clauses;

  // for (int p = 0; p < _m; p++)
  // {
  //   /// Resize container for each mutation
  //   _falsePositiveVars[p].resize(_n);
  //   _falseNegativeVars[p].resize(_n);

  //   for (int c = 0; c < _n; ++c)
  //   {
  //     if (_B.getEntry(p, c) == 0){
  //       // add a false negative variable
  //       _falseNegativeVars[p][c] = _nrActiveVariables;
  //       _nrActiveVariables++;

  //       _falsePositiveVars[p][c] = 0;

  //       // prevent "false negative and entry is 2"
  //       std::vector<Lit> clause;
  //       clause.push_back(Lit(_B2Var[p][c], false));
  //       clause.push_back(Lit(_falseNegativeVars[p][c], false));

  //       clauses.push_back(clause);
  //     } else {
  //       _falsePositiveVars[p][c] = _nrActiveVariables;
  //       _nrActiveVariables++;

  //       _falseNegativeVars[p][c] = 0;

  //       // prevent "not false positive and entry is 2"
  //       std::vector<Lit> clause;
  //       clause.push_back(Lit(_B2Var[p][c], false));
  //       clause.push_back(Lit(_falsePositiveVars[p][c], true));

  //       clauses.push_back(clause);
  //     }
  //   }
  // }

  _cuttingPlane = new CuttingPlaneDolloOld(_approxmc->get_solver(), _B, _m, _n, _k, _B2Var, _activeEntries, _solA);
  _unigen->set_cutting_plane(_cuttingPlane);
  _approxmc->setCuttingPlane(_cuttingPlane);

  /// Create variables
  SATSolver* solver = _approxmc->get_solver();
  solver->new_vars(_nrActiveVariables);

  std::vector<uint32_t> sampling_set;
  
  // Update sampling set
  for (int i = 0; i < _nrActiveVariables; ++i)
  {
    sampling_set.push_back(i);
  }
  _approxmc->set_sampling_set(sampling_set);
}

int SamplerDolloOld::solve(const ApproxMC::SolCount *sol_count, uint32_t num_samples) {
  _unigen->sample(sol_count, num_samples); // _approxmc->solve(_conf);
  // if (num_solutions > 0) {
  //   processSolution();
  // }
  processSolution();

  return 0;
}

/// Get current assignment from solver and input
lbool SamplerDolloOld::getAssignment(int var) {
  SATSolver* solver = _approxmc->get_solver();
  return solver->get_model()[var]; // see if this is const
}

int SamplerDolloOld::getEntryAssignment(int p, int c){
  
  if (!_activeEntries[p][c]){
    return _B.getEntry(p, c);
  }
  
  int var = _B2Var[p][c];

  if (getAssignment(var) == l_True) {
    return 2;
  }
  
  if (getAssignment(var) == l_False) {
    return 0;
  }

  throw std::runtime_error("Error: Solver did not assign truth value to variable.");
}

void SamplerDolloOld::processSolution() {
  for (int p = 0; p < _m; p++){
    for (int c = 0; c < _n; ++c){
      _solA.setEntry(p, c, getEntryAssignment(p, c));
    }
  }
}
