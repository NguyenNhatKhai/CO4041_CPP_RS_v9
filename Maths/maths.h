////////////////////////////////////////////////////////////////////////////////////////////////////
//
// File: maths.h
// Author: Nhat Khai Nguyen
//
////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef _MATHS_H_
#define _MATHS_H_

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

class Field;
class Element;
class Polynomial;
class Matrix;

////////////////////////////////////////////////////////////////////////////////////////////////////

class Polynomial {
    public:
    Field* field;
    vector<Element> coefficients;

    public:
    Polynomial();
    Polynomial(Field* field, vector<Element> coefficients);
    ~Polynomial() = default;

    public:
    bool operator==(const Polynomial& polynomial) const;
    bool operator!=(const Polynomial& polynomial) const;
    Polynomial operator-() const;
    Polynomial operator*(Element scalar) const;
    Polynomial operator+(const Polynomial& polynomial) const;
    Polynomial operator-(const Polynomial& polynomial) const;
    Polynomial operator*(const Polynomial& polynomial) const;
    Polynomial operator/(const Polynomial& polynomial) const;
    Polynomial operator%(const Polynomial& polynomial) const;

    public:
    int degree() const;
    Polynomial redegree(int degree) const;
    Polynomial align() const;
    Element evaluate(Element argument) const;
    Polynomial derivative() const;
};

ostream& operator<<(ostream& output, const Polynomial& polynomial);

namespace polynomials {
    static Polynomial default_polynomial;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// class Matrix {
//     public:
//     Field* field;
//     vector<vector<Element>> elements;

//     public:
//     Matrix();
//     Matrix(Field* field, vector<vector<Element>> elements);
//     ~Matrix() = default;

//     public:
//     bool operator==(const Matrix& matrix) const;
//     bool operator!=(const Matrix& matrix) const;
//     Matrix operator-() const;
//     Matrix operator~() const;
//     Matrix operator*(Element scalar) const;
//     Matrix operator+(const Matrix& matrix) const;
//     Matrix operator-(const Matrix& matrix) const;
//     Matrix operator*(const Matrix& matrix) const;
//     Matrix operator/(const Matrix& matrix) const;

//     public:
//     int row() const;
//     int column() const;
//     int size() const;
//     Matrix resize(int row, int column) const;

//     public:
//     Matrix transpose() const;
//     Element determinant() const;
// };

// ostream& operator<<(ostream& output, const Matrix& matrix);

// namespace matrices {
//     static Matrix default_matrix;
// }

////////////////////////////////////////////////////////////////////////////////////////////////////

#include "../FFA/ffa.h"
// #include "matrix.cpp"
#include "polynomial.cpp"

#endif

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////