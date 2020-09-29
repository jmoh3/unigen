/*
 * utils.h
 *
 */

#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <stdexcept>
#include <set>
#include <map>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>


typedef std::vector<bool> StlBoolVector;
typedef std::vector<double> StlDoubleVector;
typedef std::vector<int> StlIntVector;
typedef std::vector<StlIntVector> StlIntMatrix;
typedef std::vector<StlIntMatrix> StlInt3Matrix;
typedef std::set<int> StlIntSet;
typedef std::vector<std::string> StringVector;
typedef std::pair<int, int> IntPair;
typedef std::set<IntPair> IntPairSet;
typedef std::vector<IntPair> IntPairVector;
typedef std::vector<StlBoolVector> StlBoolMatrix;
typedef std::vector<StlBoolMatrix> StlBool3Matrix;

/// Platform independent extraction of line from input string
///
/// @param is Input stream
/// @param t Output string
std::istream& getline(std::istream& is, std::string& t);

/// Return a string containing the current line number
std::string getLineNumber();

/// Current line number
extern int g_lineNumber;

#endif // UTILS_H
