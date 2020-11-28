#include "adder.h"
#include <cryptominisat5/cryptominisat.h>
#include <bitset>

using std::vector;
using CMSat::Lit;

Adder::Adder(int start_var) {
    true_var_ = start_var;
    current_var_ = start_var + 1;

    SetTrue(true_var_);
}

void Adder::EncodeEqualToK(const vector<int>& vars_to_sum, size_t k) {
    vector<vector<int>> to_sum;
    for (auto var : vars_to_sum) {
        vector<int> var_vec { var };
        to_sum.push_back(var_vec);
    }

    while (to_sum.size() > 1) {
        vector<vector<int>> result;
        size_t end = to_sum.size() - (to_sum.size() % 2);

        for (size_t i = 0; i < end; i += 2) {
            to_sum[i].push_back(GetFalseVar()); // pad by 1
            to_sum[i+1].push_back(GetFalseVar()); // pad by 1
            result.push_back(Add(to_sum[i], to_sum[i+1]));
        }
        
        if (to_sum.size() % 2 != 0) {
            to_sum[to_sum.size() - 1].push_back(GetFalseVar());
            result.push_back(to_sum[to_sum.size() - 1]);
        }
        
        to_sum = result;
    }

    vector<int> k_binary = ConvertIntToBinList(k, to_sum[0].size());
    Equal(to_sum[0], k_binary);
}

void Adder::EncodeLeqToK(const vector<int>& vars_to_sum, size_t k) {
    vector<vector<int>> to_sum;
    for (auto var : vars_to_sum) {
        vector<int> var_vec { var };
        to_sum.push_back(var_vec);
    }

    while (to_sum.size() > 1) {
        vector<vector<int>> result;
        size_t end = to_sum.size() - (to_sum.size() % 2);

        for (size_t i = 0; i < end; i += 2) {
            to_sum[i].push_back(GetFalseVar()); // pad by 1
            to_sum[i+1].push_back(GetFalseVar()); // pad by 1
            result.push_back(Add(to_sum[i], to_sum[i+1]));
        }
        
        if (to_sum.size() % 2 != 0) {
            to_sum[to_sum.size() - 1].push_back(GetFalseVar());
            result.push_back(to_sum[to_sum.size() - 1]);
        }
        
        to_sum = result;
    }

    vector<int> k_binary = ConvertIntToBinList(k, to_sum[0].size());
    Leq(to_sum[0], k_binary);
}

int Adder::Leq(const vector<int>& a, const vector<int>& b) {
    size_t num_bits = a.size();
    vector<int> comp_a;
    for (auto var : a) {
        comp_a.push_back(-var);
    }

    vector<int> r = GetNewVarVector(num_bits);
    vector<int> c = GetNewVarVector(num_bits);

    FullAdder(comp_a[0], b[0], GetTrueVar(), r[0], c[0]);
    for (size_t i = 1; i < num_bits; i++) {
        FullAdder(comp_a[i], b[i], c[i-1], r[i], c[i]);
    }

    vector<int> last_clause { c[num_bits-1] };
    clauses_.push_back(ConvertIntVecToClause(last_clause));

    return c[num_bits-1];
}

void Adder::Equal(const vector<int>& a, const vector<int>& b) {
    size_t num_bits = a.size();

    for (size_t i = 0; i < num_bits; i++) {
        vector<int> clause_1 { -a[i], b[i] };
        vector<int> clause_2 { a[i], -b[i] };

        clauses_.push_back(ConvertIntVecToClause(clause_1));
        clauses_.push_back(ConvertIntVecToClause(clause_2));
    }
}

vector<int> Adder::Add(const vector<int>& a, const vector<int>& b) {
    size_t num_bits = a.size();

    vector<int> r = GetNewVarVector(num_bits);
    vector<int> c = GetNewVarVector(num_bits);

    HalfAdder(a[0], b[0], r[0], c[0]);

    for (size_t i = 1; i < num_bits; i++) {
        FullAdder(a[i], b[i], c[i-1], r[i], c[i]);
    }

    vector<int> last_clause { -c[num_bits-1] };
    clauses_.push_back(ConvertIntVecToClause(last_clause));

    return r;
}

void Adder::HalfAdder(int a, int b, int result, int carry) {
    XOR(a, b, &result);
    AND(a, b, &carry);
}

void Adder::FullAdder(int a, int b, int c, int result, int carry) {
    int r1 = GetNewVar();
    int c1 = GetNewVar();
    HalfAdder(a, b, r1, c1);

    int c2 = GetNewVar();
    HalfAdder(r1, c, result, c2);
    
    OR(c1, c2, &carry);
}

int Adder::AND(int a, int b, int* r = nullptr) {
    if (r == nullptr) {
        int new_var = GetNewVar();
        r = &new_var;
    }
    int r_var = *r;

    vector<int> clause_1 { r_var, -a, -b };
    vector<int> clause_2 { -r_var, a };
    vector<int> clause_3 { -r_var, b };

    clauses_.push_back(ConvertIntVecToClause(clause_1));
    clauses_.push_back(ConvertIntVecToClause(clause_2));
    clauses_.push_back(ConvertIntVecToClause(clause_3));

    return r_var;
}

int Adder::OR(int a, int b, int* r = nullptr) {
    if (r == nullptr) {
        int new_var = GetNewVar();
        r = &new_var;
    }
    int r_var = *r;

    vector<int> clause_1 { -r_var, a, b };
    vector<int> clause_2 { r_var, -a };
    vector<int> clause_3 { r_var, -b };

    clauses_.push_back(ConvertIntVecToClause(clause_1));
    clauses_.push_back(ConvertIntVecToClause(clause_2));
    clauses_.push_back(ConvertIntVecToClause(clause_3));

    return r_var;
}

int Adder::XOR(int a, int b, int* r = nullptr) {
    if (r == nullptr) {
        int new_var = GetNewVar();
        r = &new_var;
    }
    int r_var = *r;

    vector<int> clause_1 { -r_var, a, b };
    vector<int> clause_2 { -r_var, -a, -b };
    vector<int> clause_3 { r_var, a, -b };
    vector<int> clause_4 { r_var, -a, b };

    clauses_.push_back(ConvertIntVecToClause(clause_1));
    clauses_.push_back(ConvertIntVecToClause(clause_2));
    clauses_.push_back(ConvertIntVecToClause(clause_3));
    clauses_.push_back(ConvertIntVecToClause(clause_4));

    return r_var;
}

vector<int> Adder::GetNewVarVector(size_t num_vars) {
    vector<int> new_vars;

    for (size_t i = 0; i < num_vars; i++) {
        new_vars.push_back(GetNewVar());
    }

    return new_vars;
}

int Adder::GetNewVar(bool independent=false) {
    int new_var = current_var_;
    if (independent) {
        independent_set_.push_back(new_var);
    }
    current_var_++;
    return new_var;
}

void Adder::SetTrue(int var) {
    Lit true_lit(var, false);
    vector<Lit> true_clause { true_lit };
    clauses_.push_back(true_clause);
}

vector<Lit> Adder::ConvertIntVecToClause(const vector<int>& clause_ints) const {
    vector<Lit> clause;

    for (auto var : clause_ints) {
        if (var > 0) {
            Lit lit(var, false);
            clause.push_back(lit);
        } else {
            Lit lit(-var, true);
            clause.push_back(lit);
        }
    }

    return clause;
}

vector<int> Adder::ConvertIntToBinList(int val, size_t num_bits) const {
    vector<int> binary_vec;

    const size_t N = 32;
    assert(num_bits <= N);
    
    // credit to https://stackoverflow.com/questions/22746429/c-decimal-to-binary-converting
    std::string binary_str = std::bitset<N>(val).to_string();

    for (size_t i = N-1; i >= N - num_bits; i--) {
        char curr_char = binary_str[i];
        if (curr_char == '1') {
            binary_vec.push_back(GetTrueVar());
        } else {
            binary_vec.push_back(GetFalseVar());
        }
    }

    return binary_vec;
}