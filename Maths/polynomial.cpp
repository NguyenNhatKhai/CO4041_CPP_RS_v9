////////////////////////////////////////////////////////////////////////////////////////////////////
//
// File: polynomial.cpp
// Author: Nhat Khai Nguyen
//
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "maths.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

Polynomial::Polynomial() {
    this->field = &fields::default_field;
    this->coefficients = vector<Element>();
}

Polynomial::Polynomial(Field* field, vector<Element> coefficients) {
    if (coefficients.size() == 0) {
        throw "Maths\\Polynomial\\Polynomial(Field*, vector<Element>)\\coefficients\\size";
    }
    for (int i = 0; i < coefficients.size(); i ++) {
        if (*coefficients[i].field != *field) {
            throw "Maths\\Polynomial\\Polynomial(Field*, vector<Element>)\\coefficients\\field";
        }
    }
    this->field = field;
    this->coefficients = coefficients;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

bool Polynomial::operator==(const Polynomial& polynomial) const {
    if (this == &polynomials::default_polynomial && &polynomial != &polynomials::default_polynomial || this != &polynomials::default_polynomial && &polynomial == &polynomials::default_polynomial) return false;
    if (this->field == polynomial.field && this->align().coefficients == polynomial.align().coefficients) return true;
    return *this->field == *polynomial.field && this->align().coefficients == polynomial.align().coefficients;
}

bool Polynomial::operator!=(const Polynomial& polynomial) const {
    return !(*this == polynomial);
}

ostream& operator<<(ostream& output, const Polynomial& polynomial) {
    output << "{";
    for (int i = 0; i < polynomial.coefficients.size(); i ++) {
        output << polynomial.coefficients[i];
        if (i != polynomial.coefficients.size() - 1) {
            output << ", ";
        }
    }
    output << "}";
    return output;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

int Polynomial::degree() const {
    for (int i = 0; i < this->coefficients.size(); i ++) {
        if (this->coefficients[this->coefficients.size() - i - 1] != this->field->zero_element()) {
            return this->coefficients.size() - i - 1;
        }
    }
    return 0;
}

Polynomial Polynomial::redegree(int degree) const {
    if (*this == polynomials::default_polynomial) {
        throw "Maths\\Polynomial\\redegree(int)";
    }
    if (degree < 0) {
        throw "Maths\\Polynomial\\redegree(int)\\degree";
    }
    vector<Element> new_coefficients(degree + 1, this->field->zero_element());
    for (int i = 0; i <= min(degree, this->degree()); i ++) {
        new_coefficients[i] = this->coefficients[i];
    }
    return Polynomial(this->field, new_coefficients);
}

Polynomial Polynomial::align() const {
    if (*this == polynomials::default_polynomial) {
        throw "Maths\\Polynomial\\align()";
    }
    return this->redegree(this->degree());
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

Polynomial Polynomial::operator+(const Polynomial& polynomial) const {
    if (*this == polynomials::default_polynomial) {
        throw "Maths\\Polynomial\\operator+(const Polynomial&)";
    }
    if (polynomial == polynomials::default_polynomial) {
        throw "Maths\\Polynomial\\operator+(const Polynomial&)\\polynomial";
    }
    if (*polynomial.field != *this->field) {
        throw "Maths\\Polynomial\\operator+(const Polynomial&)\\polynomial\\field";
    }
    int new_degree = max(this->degree(), polynomial.degree());
    Polynomial temp_0_polynomial = this->redegree(new_degree);
    Polynomial temp_1_polynomial = polynomial.redegree(new_degree);
    vector<Element> new_coefficients(new_degree + 1, this->field->zero_element());
    for (int i = 0; i <= new_degree; i ++) {
        new_coefficients[i] = temp_0_polynomial.coefficients[i] + temp_1_polynomial.coefficients[i];
    }
    return Polynomial(this->field, new_coefficients);
}

Polynomial Polynomial::operator*(const Polynomial& polynomial) const {
    if (*this == polynomials::default_polynomial) {
        throw "Maths\\Polynomial\\operator*(const Polynomial&)";
    }
    if (polynomial == polynomials::default_polynomial) {
        throw "Maths\\Polynomial\\operator*(const Polynomial&)\\polynomial";
    }
    if (polynomial.field != this->field) {
        throw "Maths\\Polynomial\\operator*(const Polynomial&)\\polynomial\\field";
    }
    int new_degree = this->degree() + polynomial.degree();
    Polynomial temp_0_polynomial = this->align();
    Polynomial temp_1_polynomial = polynomial.align();
    vector<Element> new_coefficients(new_degree + 1, this->field->zero_element());
    for (int i = 0; i <= temp_0_polynomial.degree(); i ++) {
        for (int j = 0; j <= temp_1_polynomial.degree(); j ++) {
            new_coefficients[i + j] = new_coefficients[i + j] + (temp_0_polynomial.coefficients[i] * temp_1_polynomial.coefficients[j]);
        }
    }
    return Polynomial(this->field, new_coefficients);
}

Polynomial Polynomial::operator/(const Polynomial& polynomial) const {
    if (*this == polynomials::default_polynomial) {
        throw "Maths\\Polynomial\\operator/(const Polynomial&)";
    }
    if (polynomial == polynomials::default_polynomial) {
        throw "Maths\\Polynomial\\operator/(const Polynomial&)\\polynomial";
    }
    if (polynomial.field != this->field) {
        throw "Maths\\Polynomial\\operator/(const Polynomial&)\\polynomial\\field";
    }
    int new_degree = this->degree() - polynomial.degree();
    if (new_degree < 0) return Polynomial(this->field, {this->field->zero_element()});
    Polynomial temp_0_polynomial = this->align();
    Polynomial temp_1_polynomial = polynomial.align();
    vector<Element> new_coefficients(new_degree + 1, this->field->zero_element());
    for (int i = 0; i <= new_degree; i ++) {
        new_coefficients[new_degree - i] = temp_0_polynomial.coefficients[this->degree() - i] / temp_1_polynomial.coefficients[temp_1_polynomial.degree()];
        temp_0_polynomial = *this + (temp_1_polynomial * Polynomial(this->field, new_coefficients));
    }
    return Polynomial(this->field, new_coefficients);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

Polynomial Polynomial::operator-() const {
    if (*this == polynomials::default_polynomial) {
        throw "Maths\\Polynomial\\operator-()";
    }
    int new_degree = this->degree();
    vector<Element> new_coefficients(new_degree + 1, this->field->zero_element());
    for (int i = 0; i <= new_degree; i ++) {
        new_coefficients[i] = -this->coefficients[i];
    }
    return Polynomial(this->field, new_coefficients);
}

Polynomial Polynomial::operator-(const Polynomial& polynomial) const {
    if (*this == polynomials::default_polynomial) {
        throw "Maths\\Polynomial\\operator-(const Polynomial&)";
    }
    if (polynomial == polynomials::default_polynomial) {
        throw "Maths\\Polynomial\\operator-(const Polynomial&)\\polynomial";
    }
    if (polynomial.field != this->field) {
        throw "Maths\\Polynomial\\operator-(const Polynomial&)\\polynomial\\field";
    }
    return (*this) + (-polynomial);
}

Polynomial Polynomial::operator%(const Polynomial& polynomial) const {
    if (*this == polynomials::default_polynomial) {
        throw "Maths\\Polynomial\\operator%(const Polynomial&)";
    }
    if (polynomial == polynomials::default_polynomial) {
        throw "Maths\\Polynomial\\operator%(const Polynomial&)\\polynomial";
    }
    if (polynomial.field != this->field) {
        throw "Maths\\Polynomial\\operator%(const Polynomial&)\\polynomial\\field";
    }
    int new_degree = polynomial.degree() - 1;
    if (new_degree < 0) return Polynomial(this->field, {this->field->zero_element()});
    return (*this - (polynomial * (*this / polynomial))).redegree(new_degree);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////