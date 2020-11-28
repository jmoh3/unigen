#include "adder.h"
#include <cryptominisat5/cryptominisat.h>

using std::vector;
using CMSat::Lit;

Adder::Adder(int start_var) {
    true_var_ = start_var;
    current_var_ = start_var + 1;

    SetTrue(true_var_);
}

void Adder::EncodeEqualToK(const vector<int>& vars_to_sum, size_t k) {

}

void Adder::EncodeLeqToK(const vector<int>& vars_to_sum, size_t k) {

}

int Adder::Leq(const vector<int>& a, const vector<int>& b) {

}

void Adder::Equal(const vector<int>& a, const vector<int>& b) {

}

vector<int> Adder::Add(const vector<int>& a, const vector<int>& b) {

}

void Adder::HalfAdder(int a, int b, int result, int carry) {

}

void Adder::FullAdder(int a, int b, int c, int result, int carry) {

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