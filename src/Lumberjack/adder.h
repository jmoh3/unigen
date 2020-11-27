/**
 * adder.h 
 * 
 */
#ifndef ADDER_H
#define ADDER_H

#include <cryptominisat5/cryptominisat.h>
#include "matrix.h"
#include "utils.h"
#include <approxmc/cuttingplane.h>
#include <utility>
#include <vector>

using std::vector;

class Adder {
  public:
    Adder(int start_var);

    /**
     * Enforces that the sum of boolean variables in vars_to_sum must equal k
     * @param vars_to_sum 
     * @param k
     */
    void EncodeEqualToK(const vector<int>& vars_to_sum, size_t k);

    /**
     * Enforces that the sum of boolean variables in vars_to_sum must be less than or equal to k
     * @param vars_to_sum 
     * @param k
     */
    void EncodeLeqToK(const vector<int>& vars_to_sum, size_t k);
    
    /**
     * Encodes addition of two numbers, each encoded as vector of binary variables.
     * @param a first number, encoded as a vector of binary variables
     * @param b second number, encoded as a vector of binary variables
     * @return a vector of binary variables representing the addition of these two numbers
     */
    vector<int> Add(const vector<int>& a, const vector<int>& b);

    /**
     * Encodes a less-than-or-equal-to constraint: a <= b
     * @param a the first number, encoded as a vector of binary variables
     * @param b second number, which must be greater than a, encoded as a vector of binary variables
     * @return a variable that will represent whether this constraint is satisfied or not
     */
    int Leq(const vector<int>& a, const vector<int>& b);

    /**
     * Encodes an equality constraint: a = b
     * @param a the first number, encoded as a vector of binary variables
     * @param b the first number, encoded as a vector of binary variables
     */
    void Equal(const vector<int>& a, const vector<int>& b);

    const vector<vector<Lit>>& GetClauses() const;
    const vector<int>& GetIndependentSet() const;
  
  protected:
    int GetTrueVar() const;
    int GetFalseVar() const;

    int GetNewVar(bool independent=false);

    int AND(int a, int b, int r);
    int OR(int a, int b, int r);
    int XOR(int a, int b, int r);

    void SetTrue(int var);

    void HalfAdder(int a, int b, int result, int carry);
    void FullAdder(int a, int b, int c, int result, int carry);
  
  private:
    vector<vector<Lit>> clauses_;
    vector<int> independent_set_;
    int current_var_;
    /** Reserved variable that is set to true in CNF formula */
    int true_var_;
};

#endif // ADDER_H