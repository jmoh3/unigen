/*
 * cuttingplanedollo.cpp
 *
 */

#include "cuttingplanedollo.h"
using namespace CMSat;

CuttingPlaneDollo::CuttingPlaneDollo(SATSolver* solver, const Matrix& B, const int m, const int n, const int k, StlIntMatrix& B2Var, StlBoolMatrix& activeEntries, Matrix& solA)
  : CuttingPlane(solver)
  , _B(B)
  , _m(m)
  , _n(n)
  , _k(k)
  , _B2Var(B2Var)
  , _activeEntries(activeEntries)
  , _solA(solA)
{
}

int CuttingPlaneDollo::getEntryAssignment(int p, int c){
  
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

int CuttingPlaneDollo::separate()
{
  int nrNewClauses = 0;
  std::vector<Lit> clause;
  
  for (int c = 0; c < _n; ++c)
  {
    for (int d = c + 1; d < _n; ++d)
    {
      // condition 1
      for (int i = 1; i <= _k + 1; ++i)
      {
        for (int i_prime = 1; i_prime <= _k + 1; ++i_prime)
        {
          for (int j = 1; j <= _k + 1; ++j)
          {
            for (int j_prime = 1; j_prime <= _k + 1; ++j_prime)
            {
              StlIntSet pSet;
              for (int p = 0; p < _m; p++)
              {
                if ((getEntryAssignment(p,c) == i) && (getEntryAssignment(p,d) == 0))
                {
                  pSet.insert(p);
                }
              }
              
              StlIntSet qSet;
              for (int q = 0; q < _m; q++)
              {
                if ((getEntryAssignment(q,c) == 0) && (getEntryAssignment(q,d) == j))
                {
                  qSet.insert(q);
                }
              }
              
              StlIntSet rSet;
              for (int r = 0; r < _m; r++)
              {
                if ((getEntryAssignment(r,c) == i_prime) && (getEntryAssignment(r,d) == j_prime))
                {
                  rSet.insert(r);
                }
              }
              
              for (int p : pSet)
              {
                for (int q : qSet)
                {
                  for (int r : rSet)
                  {
                    clause.clear();
                                        
                    if (_activeEntries[p][c]){
                      /// constraint[0] = Triple(p, c, i);
                      /// This entry is currently set to TRUE
                      /// This constraint will not be violated if it is FALSE next time
                      clause.push_back(Lit(_B2Var[p][c], true));
                    }
                    
                    if (_activeEntries[p][d]){
                      /// constraint[1] = Triple(p, d, 0);
                      /// This entry is currently set to FALSE
                      /// This constraint will not be violated if it is TRUE next time
                      clause.push_back(Lit(_B2Var[p][d], false));
                    }
                    
                    if (_activeEntries[q][c]){
                      /// constraint[2] = Triple(q, c, 0);
                      /// This entry is currently set to FALSE
                      /// This constraint will not be violated if it is TRUE next time
                      clause.push_back(Lit(_B2Var[q][c], false));
                    }
                    
                    if (_activeEntries[q][d]){
                      /// constraint[3] = Triple(q, d, j);
                      /// This entry is currently set to TRUE
                      /// This constraint will not be violated if it is FALSE next time
                      clause.push_back(Lit(_B2Var[q][d], true));
                    }
                    
                    if (_activeEntries[r][c]){
                      /// constraint[4] = Triple(r, c, i_prime);
                      /// This entry is currently set to TRUE
                      /// This constraint will not be violated if it is FALSE next time
                      clause.push_back(Lit(_B2Var[r][c], true));
                    }
                    
                    if (_activeEntries[r][d]){
                      /// constraint[5] = Triple(r, d, j_prime);
                      /// This entry is currently set to TRUE
                      /// This constraint will not be violated if it is FALSE next time
                      clause.push_back(Lit(_B2Var[r][d], true));
                    }
//                    sortClause(clause);
                    addClause(clause);
                    nrNewClauses = nrNewClauses+1;
                  }
                }
              }
            }
          }
        }
      }
      
      // condition 2
      for (int i = 1; i <= _k + 1; ++i)
      {
        for (int i_prime = 1; i_prime <= _k + 1; ++i_prime)
        {
          for (int j = 2; j <= _k + 1; ++j)
          {
            for (int j_prime = 1; j_prime <= _k + 1; ++j_prime)
            {
              if (j_prime == j) continue;
          
              StlIntSet pSet;
              for (int p = 0; p < _m; p++)
              {
                if ((getEntryAssignment(p,c) == i) && (getEntryAssignment(p,d) == j_prime))
                {
                  pSet.insert(p);
                }
              }
              
              StlIntSet qSet;
              for (int q = 0; q < _m; q++)
              {
                if ((getEntryAssignment(q,c) == 0) && (getEntryAssignment(q,d) == j))
                {
                  qSet.insert(q);
                }
              }
              
              StlIntSet rSet;
              for (int r = 0; r < _m; r++)
              {
                if ((getEntryAssignment(r,c) == i_prime) && (getEntryAssignment(r,d) == j))
                {
                  rSet.insert(r);
                }
              }
              
              for (int p : pSet)
              {
                for (int q : qSet)
                {
                  for (int r : rSet)
                  {
                    clause.clear();
                    
                    if (_activeEntries[p][c]){
                      /// constraint[0] = Triple(p, c, i);
                      clause.push_back(Lit(_B2Var[p][c], true));
                    }
                    
                    if (_activeEntries[p][d]){
                      /// constraint[1] = Triple(p, d, j_prime);
                      clause.push_back(Lit(_B2Var[p][d], true));
                    }
                    
                    if (_activeEntries[q][c]){
                      /// constraint[2] = Triple(q, c, 0);
                      clause.push_back(Lit(_B2Var[q][c], false));
                    }
                    
                    if (_activeEntries[q][d]){
                      /// constraint[3] = Triple(q, d, j);
                      clause.push_back(Lit(_B2Var[q][d], true));
                    }
                    
                    if (_activeEntries[r][c]){
                      /// constraint[4] = Triple(r, c, i_prime);
                      clause.push_back(Lit(_B2Var[r][c], true));
                    }
                    
                    if (_activeEntries[r][d]){
                      /// constraint[5] = Triple(r, d, j);
                      clause.push_back(Lit(_B2Var[r][d], true));
                    }
//                    sortClause(clause);
                    addClause(clause);
                    nrNewClauses = nrNewClauses+1;
                  }
                }
              }
            }
          }
        }
      }
      
      // condition 3
      for (int i = 2; i <= _k + 1; ++i)
      {
        for (int i_prime = 1; i_prime <= _k + 1; ++i_prime)
        {
          if (i_prime == i) continue;
          for (int j = 1; j <= _k + 1; ++j)
          {
            for (int j_prime = 1; j_prime <= _k + 1; ++j_prime)
            {
              StlIntSet pSet;
              for (int p = 0; p < _m; p++)
              {
                if ((getEntryAssignment(p,c) == i) && (getEntryAssignment(p,d) == 0))
                {
                  pSet.insert(p);
                }
              }
              
              StlIntSet qSet;
              for (int q = 0; q < _m; q++)
              {
                if ((getEntryAssignment(q,c) == i_prime) && (getEntryAssignment(q,d) == j))
                {
                  qSet.insert(q);
                }
              }
              
              StlIntSet rSet;
              for (int r = 0; r < _m; r++)
              {
                if ((getEntryAssignment(r,c) == i) && (getEntryAssignment(r,d) == j_prime))
                {
                  rSet.insert(r);
                }
              }
              
              for (int p : pSet)
              {
                for (int q : qSet)
                {
                  for (int r : rSet)
                  {
                    clause.clear();
                    
                    if (_activeEntries[p][c]){
                      /// constraint[0] = Triple(p, c, i);
                      clause.push_back(Lit(_B2Var[p][c], true));
                    }
                    
                    if (_activeEntries[p][d]){
                      /// constraint[1] = Triple(p, d, 0);
                      clause.push_back(Lit(_B2Var[p][d], false));
                    }
                    
                    if (_activeEntries[q][c]){
                      /// constraint[2] = Triple(q, c, i_prime);
                      clause.push_back(Lit(_B2Var[q][c], true));
                    }
                    
                    if (_activeEntries[q][d]){
                      /// constraint[3] = Triple(q, d, j);
                      clause.push_back(Lit(_B2Var[q][d], true));
                    }
                    
                    if (_activeEntries[r][c]){
                      /// constraint[4] = Triple(r, c, i);
                      clause.push_back(Lit(_B2Var[r][c], true));
                    }
                    
                    if (_activeEntries[r][d]){
                      /// constraint[5] = Triple(r, d, j_prime);
                      clause.push_back(Lit(_B2Var[r][d], true));
                    }
//                    sortClause(clause);
                    addClause(clause);
                    nrNewClauses = nrNewClauses+1;
                  }
                }
              }
            }
          }
        }
      }
      
      // condition 4
      for (int i = 2; i <= _k + 1; ++i)
      {
        for (int i_prime = 1; i_prime <= _k + 1; ++i_prime)
        {
          if (i_prime == i) continue;
          for (int j = 2; j <= _k + 1; ++j)
          {
            for (int j_prime = 1; j_prime <= _k + 1; ++j_prime)
            {
              if (j_prime == j) continue;
              
              StlIntSet pSet;
              for (int p = 0; p < _m; p++)
              {
                if ((getEntryAssignment(p,c) == i) && (getEntryAssignment(p,d) == j_prime))
                {
                  pSet.insert(p);
                }
              }
              
              StlIntSet qSet;
              for (int q = 0; q < _m; q++)
              {
                if ((getEntryAssignment(q,c) == i_prime) && (getEntryAssignment(q,d) == j))
                {
                  qSet.insert(q);
                }
              }
              
              StlIntSet rSet;
              for (int r = 0; r < _m; r++)
              {
                if ((getEntryAssignment(r,c) == i) && (getEntryAssignment(r,d) == j))
                {
                  rSet.insert(r);
                }
              }
              
              for (int p : pSet)
              {
                for (int q : qSet)
                {
                  for (int r : rSet)
                  {
                    clause.clear();
                    
                    if (_activeEntries[p][c]){
                      /// constraint[0] = Triple(p, c, i);
                      clause.push_back(Lit(_B2Var[p][c], true));
                    }
                    
                    if (_activeEntries[p][d]){
                      /// constraint[1] = Triple(p, d, j_prime);
                      clause.push_back(Lit(_B2Var[p][d], true));
                    }
                    
                    if (_activeEntries[q][c]){
                      /// constraint[2] = Triple(q, c, i_prime);
                      clause.push_back(Lit(_B2Var[q][c], true));
                    }
                    
                    if (_activeEntries[q][d]){
                      /// constraint[3] = Triple(q, d, j);
                      clause.push_back(Lit(_B2Var[q][d], true));
                    }
                    
                    if (_activeEntries[r][c]){
                      /// constraint[4] = Triple(r, c, i);
                      clause.push_back(Lit(_B2Var[r][c], true));
                    }
                    
                    if (_activeEntries[r][d]){
                      /// constraint[5] = Triple(r, d, j);
                      clause.push_back(Lit(_B2Var[r][d], true));
                    }
//                    sortClause(clause);
                    addClause(clause);
                    nrNewClauses = nrNewClauses+1;
                  }
                }
              }
            }
          }
        }
      }
    }
  }
//  std::cout << "Added " << nrNewClauses << " clauses\n";
  return nrNewClauses;
}

void CuttingPlaneDollo::processSolution(){
  for (int p = 0; p < _m; p++){
    for (int c = 0; c < _n; ++c){
      _solA.setEntry(p, c, getEntryAssignment(p, c));
    }
  }
}
