/*
 * cuttingplane.h
 *
 */

#ifndef CUTTINGPLANE_H
#define CUTTINGPLANE_H

#include <fstream>
#include <random>
#include <map>
#include <cstdint>
#include <cryptominisat5/cryptominisat.h>

using namespace CMSat;

namespace UniGen
{
    /// This class provides a cutting plane wrapper for ApproxMC
  class CuttingPlane
  {
  public:

    /// Constructor
    CuttingPlane(SATSolver* solver): _solver(solver) {}

    virtual ~CuttingPlane() {}

    /// Identify violated constraints
    virtual int separate() = 0;

    int getNumClausesAdded() const {
      return _numClausesAdded;
    }
    
  protected:
    
    /// Get current assignment from solver and input
    lbool getAssignment(int var) {
      return _solver->get_model()[var]; // see if this is const
    }

    void addClause(const std::vector<Lit>& clause) {
      std::cout << "Added clause " << clause << std::endl;
      _solver->add_clause(clause);
      _numClausesAdded++;
    }

  private:
    SATSolver* _solver;
    int _numClausesAdded;
  };
}

#endif
