/**
 * adder.h 
 * 
 */
#ifndef ADDER_H
#define ADDER_H

#include <cryptominisat5/cryptominisat.h>
#include <utility>
#include <vector>

using std::vector;
using CMSat::Lit;

class Adder {
  public:
    /**
     * Constructor
     * @param start_var the first unused variable label for the formula
     */
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

    const vector<vector<Lit>>& GetClauses() const {
      return clauses_;
    }

    const vector<int>& GetIndependentSet() const {
      return independent_set_;
    }

    size_t GetFirstUnusedVar() const {
      return current_var_;
    }

    size_t GetNumVarsAdded() const {
      return current_var_ - true_var_;
    }
  
  protected:
    int GetTrueVar() const {
      return true_var_;
    }

    int GetFalseVar() const {
      return -true_var_;
    }

    int GetNewVar(bool independent=false);
    vector<int> GetNewVarVector(size_t num_vars);

    int AND(int a, int b, int* r = nullptr);
    int OR(int a, int b, int* r = nullptr);
    int XOR(int a, int b, int* r = nullptr);

    void SetTrue(int var);

    void HalfAdder(int a, int b, int result, int carry);
    void FullAdder(int a, int b, int c, int result, int carry);

    vector<Lit> ConvertIntVecToClause(const vector<int>& clause_ints) const;

    vector<int> ConvertIntToBinList(int val, size_t num_bits) const;
  
  private:
    vector<vector<Lit>> clauses_;
    vector<int> independent_set_;
    int current_var_;
    /** Reserved variable that is set to true in CNF formula */
    int true_var_;
};

#endif // ADDER_H