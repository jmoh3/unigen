/*
 * matrix.cpp
 *
 */

#include "matrix.h"

Matrix::Matrix()
  : _m(0)
  , _n(0)
  , _D()
  , _k(0)
{
}

Matrix::Matrix(int m, int n)
  : _m(m)
  , _n(n)
  , _D(m, StlIntVector(n, 0))
{
}

Matrix* Matrix::parse(const std::string& filename)
{
  Matrix* pMatrix = NULL;
  
  if (filename == "-")
  {
    pMatrix = new Matrix();
    std::cin >> *pMatrix;
  }
  else
  {
    std::ifstream inA(filename.c_str());
    if (!inA.good())
    {
      std::cerr << "Error: could not open '" << filename << "' for reading" << std::endl;
      return NULL;
    }
    pMatrix = new Matrix();
    inA >> *pMatrix;
    inA.close();
  }
  
  return pMatrix;
}


std::ostream& operator<<(std::ostream& out, const Matrix& D)
{
  out << D._m << " #taxa" << std::endl;
  out << D._n << " #characters" << std::endl;
  for (int p = 0; p < D._m; ++p)
  {
    bool first = true;
    for (int c = 0; c < D._n; ++c)
    {
      if (first)
        first = false;
      else
        out << " ";
      
      out << D._D[p][c];
    }
    out << std::endl;
  }
  
  return out;
}

std::istream& operator>>(std::istream& in, Matrix& D)
{
  std::string line;
  getline(in, line);
  
  std::stringstream ss(line);

  int m = -1;
  ss >> m;
  if (m < 0)
  {
    throw std::runtime_error(getLineNumber()
                             + "Error: number of clones should be positive.");
  }
  D._m = m;

  int n = -1;
  getline(in, line);
  ss.clear();
  ss.str(line);
  ss >> n;
  if (n < 0)
  {
    throw std::runtime_error(getLineNumber()
                             + "Error: number of characters should be positive.");
  }
  D._n = n;
  
  D._k = 0;
  D._D = StlIntMatrix(m, StlIntVector(n, 0));
  for (int p = 0; p < m; ++p)
  {
    StringVector s;
    getline(in, line);
    boost::split(s, line, boost::is_any_of(" \t"));
    if (s.size() < (size_t)n)
    {
      throw std::runtime_error(getLineNumber()
                               + "Error: insufficient number of characters.");
    }
    
    for (int c = 0; c < n; ++c)
    {
      int i = atoi(s[c].c_str());
      D._D[p][c] = i;
      if (i - 1 > D._k)
      {
        D._k = i - 1;
      }
    }
  }
  
  return in;
}
