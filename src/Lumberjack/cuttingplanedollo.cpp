/*
 * cuttingplanedollo.cpp
 *
 */

#include "cuttingplanedollo.h"
using namespace CMSat;

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
  // TODO: implement
  return 0;
}

int CuttingPlaneDollo::separate() {
  // TODO: implement
  return 0;
}
