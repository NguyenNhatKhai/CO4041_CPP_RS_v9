////////////////////////////////////////////////////////////////////////////////////////////////////
//
// File: rs.cpp
// Author: Nhat Khai Nguyen
//
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "fec.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

RS::RS(Field* symbol_field, int correction_capability) {
    if (correction_capability <= 0 || 2 * correction_capability >= symbol_field->size() - 1) {
        throw "FEC\\RS\\RS(Field*, int)\\correction_capability";
    }
    this->symbol_field = symbol_field;
    this->generator_polynomial = Polynomial(this->symbol_field, {this->symbol_field->unit_element()});
    for (int i = 0; i < 2 * correction_capability; i ++) {
        this->generator_polynomial = this->generator_polynomial * Polynomial(this->symbol_field, {this->symbol_field->general_elements[i + 1], this->symbol_field->unit_element()});
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

int RS::codeword_length() const {
    return this->symbol_field->size() - 1;
}

int RS::message_length() const {
    return this->symbol_field->size() - this->generator_polynomial.degree() - 1;
}

int RS::parity_length() const {
    return this->codeword_length() - this->message_length();
}

int RS::symbol_size() const {
    return this->symbol_field->primitive_polynomial->degree();
}

int RS::detection_capability() const {
    return this->parity_length();
}

int RS::correction_capability() const {
    return this->parity_length() / 2;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

Polynomial RS::systematic_encode(const Polynomial& message, const string& output_path) const {
    if (message.coefficients.size() != this->message_length()) {
        throw "FEC\\RS\\systematic_encode(const Polynomial&, const string&)\\message\\coefficients\\size";
    } else if (*message.field != *this->symbol_field) {
        throw "FEC\\RS\\systematic_encode(const Polynomial&, const string&)\\message\\field";
    }
    Polynomial temp_0_polynomial(this->symbol_field, vector<Element>(this->parity_length() + 1, this->symbol_field->zero_element()));
    temp_0_polynomial.coefficients[this->parity_length()] = this->symbol_field->unit_element();
    temp_0_polynomial = ((message * temp_0_polynomial) + ((message * temp_0_polynomial) % this->generator_polynomial)).redegree(this->codeword_length() - 1);
    if (output_path != "output") {
        ofstream output_file("Output/" + output_path + ".txt");
        output_file << left << setw(16) << "polMes" << message << endl;
        output_file << left << setw(16) << "polCod" << temp_0_polynomial << endl;
        output_file.close();
    }
    return temp_0_polynomial;
}

Polynomial RS::nonsystematic_encode(const Polynomial& message, const string& output_path) const {
    if (message.coefficients.size() != this->message_length()) {
        throw "FEC\\RS\\nonsystematic_encode(const Polynomial&, const string&)\\message\\coefficients\\size";
    } else if (*message.field != *this->symbol_field) {
        throw "FEC\\RS\\nonsystematic_encode(const Polynomial&, const string&)\\message\\field";
    }
    Polynomial temp_0_polynomial = (message * this->generator_polynomial).redegree(this->codeword_length() - 1);
    if (output_path != "output") {
        ofstream output_file("Output/" + output_path + ".txt");
        output_file << left << setw(16) << "polMes" << message << endl;
        output_file << left << setw(16) << "polCod" << temp_0_polynomial << endl;
        output_file.close();
    }
    return temp_0_polynomial;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

Polynomial RS::add_error(const Polynomial& codeword, const Polynomial& error, const string& output_path) const {
    if (codeword.coefficients.size() != this->codeword_length()) {
        throw "FEC\\RS\\add_error(const Polynomial&, const Polynomial&, const string&)\\codeword\\coefficients\\size";
    } else if (error.coefficients.size() != this->codeword_length()) {
        throw "FEC\\RS\\add_error(const Polynomial&, const Polynomial&, const string&)\\error\\coefficients\\size";
    } else if (*codeword.field != *this->symbol_field) {
        throw "FEC\\RS\\add_error(const Polynomial&, const Polynomial&, const string&)\\codeword\\field";
    } else if (*error.field != *this->symbol_field) {
        throw "FEC\\RS\\add_error(const Polynomial&, const Polynomial&, const string&)\\error\\field";
    }
    Polynomial temp_0_polynomial = (codeword + error).redegree(this->codeword_length() - 1);
    if (output_path != "output") {
        ofstream output_file("Output/" + output_path + ".txt");
        output_file << left << setw(16) << "polCod" << codeword << endl;
        output_file << left << setw(16) << "polErr" << error << endl;
        output_file << left << setw(16) << "polRec" << temp_0_polynomial << endl;
        output_file.close();
    }
    return temp_0_polynomial;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

Polynomial RS::pgz_decode(const Polynomial& received, const string& output_path) const {
    if (received.coefficients.size() != this->codeword_length()) {
        throw "FEC\\RS\\pgz_decode(const Polynomial&, const string&)\\received\\coefficients\\size";
    } else if (*received.field != *this->symbol_field) {
        throw "FEC\\RS\\pgz_decode(const Polynomial&, const string&)\\received\\field";
    }
    vector<Element> temp_0_syndrome = this->syndrome(received);
    Polynomial temp_1_syndrome = this->syndrome(temp_0_syndrome);
    Polynomial temp_2_error_locator = this->pgz_error_locator(temp_0_syndrome);
    Polynomial temp_3_error_evaluator = this->pgz_error_evaluator(temp_1_syndrome, temp_2_error_locator);
    vector<Element> temp_4_roots = this->chien_roots(temp_2_error_locator);
    vector<Element> temp_5_error_values = this->pgz_error_values(temp_2_error_locator, temp_3_error_evaluator, temp_4_roots);
    Polynomial temp_6_estimated_error = this->estimated_error(temp_4_roots, temp_5_error_values);
    Polynomial temp_7_estimated_codeword = this->estimated_codeword(received, temp_6_estimated_error);
    Polynomial temp_8_estimated_message = this->estimated_message(temp_7_estimated_codeword);
    if (output_path != "output") {
        ofstream output_file("Output/" + output_path + ".txt");
        output_file << left << setw(16) << "polRec" << received << endl;
        for (int i = 0; i < temp_0_syndrome.size(); i ++) {
            output_file << left << setw(16) << "vecSyn[" + to_string(i) + "]" << temp_0_syndrome[i] << endl;
        }
        output_file << left << setw(16) << "polSyn" << temp_1_syndrome << endl;
        output_file << left << setw(16) << "polErrLoc" << temp_2_error_locator << endl;
        output_file << left << setw(16) << "polErrEva" << temp_3_error_evaluator << endl;
        for (int i = 0; i < temp_4_roots.size(); i ++) {
            output_file << left << setw(16) << "vecRoo[" + to_string(i) + "]" << temp_4_roots[i] << endl;
        }
        for (int i = 0; i < temp_5_error_values.size(); i ++) {
            output_file << left << setw(16) << "vecErrVal[" + to_string(i) + "]" << temp_5_error_values[i] << endl;
        }
        output_file << left << setw(16) << "polEstErr" << temp_6_estimated_error << endl;
        output_file << left << setw(16) << "polEstCod" << temp_7_estimated_codeword << endl;
        output_file << left << setw(16) << "polEstMes" << temp_8_estimated_message << endl;
        output_file.close();
    }
    return temp_8_estimated_message;
}

Polynomial RS::bm_decode(const Polynomial& received, const string& output_path) const {
    if (received.coefficients.size() != this->codeword_length()) {
        throw "FEC\\RS\\bm_decode(const Polynomial&, const string&)\\received\\coefficients\\size";
    } else if (*received.field != *this->symbol_field) {
        throw "FEC\\RS\\bm_decode(const Polynomial&, const string&)\\received\\field";
    }
    vector<Element> temp_0_syndrome = this->syndrome(received);
    Polynomial temp_1_syndrome = this->syndrome(temp_0_syndrome);
    Polynomial temp_2_error_locator = this->bm_error_locator(temp_0_syndrome);
    Polynomial temp_3_error_evaluator = this->bm_error_evaluator(temp_1_syndrome, temp_2_error_locator);
    vector<Element> temp_4_roots = this->chien_roots(temp_2_error_locator);
    vector<Element> temp_5_error_values = this->bm_error_values(temp_2_error_locator, temp_3_error_evaluator, temp_4_roots);
    Polynomial temp_6_estimated_error = this->estimated_error(temp_4_roots, temp_5_error_values);
    Polynomial temp_7_estimated_codeword = this->estimated_codeword(received, temp_6_estimated_error);
    Polynomial temp_8_estimated_message = this->estimated_message(temp_7_estimated_codeword);
    if (output_path != "output") {
        ofstream output_file("Output/" + output_path + ".txt");
        output_file << left << setw(16) << "polRec" << received << endl;
        for (int i = 0; i < temp_0_syndrome.size(); i ++) {
            output_file << left << setw(16) << "vecSyn[" + to_string(i) + "]" << temp_0_syndrome[i] << endl;
        }
        output_file << left << setw(16) << "polSyn" << temp_1_syndrome << endl;
        output_file << left << setw(16) << "polErrLoc" << temp_2_error_locator << endl;
        output_file << left << setw(16) << "polErrEva" << temp_3_error_evaluator << endl;
        for (int i = 0; i < temp_4_roots.size(); i ++) {
            output_file << left << setw(16) << "vecRoo[" + to_string(i) + "]" << temp_4_roots[i] << endl;
        }
        for (int i = 0; i < temp_5_error_values.size(); i ++) {
            output_file << left << setw(16) << "vecErrVal[" + to_string(i) + "]" << temp_5_error_values[i] << endl;
        }
        output_file << left << setw(16) << "polEstErr" << temp_6_estimated_error << endl;
        output_file << left << setw(16) << "polEstCod" << temp_7_estimated_codeword << endl;
        output_file << left << setw(16) << "polEstMes" << temp_8_estimated_message << endl;
        output_file.close();
    }
    return temp_8_estimated_message;
}

Polynomial RS::euclidean_decode(const Polynomial& received, const string& output_path) const {
    if (received.coefficients.size() != this->codeword_length()) {
        throw "FEC\\RS\\euclidean_decode(const Polynomial&, const string&)\\received\\coefficients\\size";
    } else if (*received.field != *this->symbol_field) {
        throw "FEC\\RS\\euclidean_decode(const Polynomial&, const string&)\\received\\field";
    }
    vector<Element> temp_0_syndrome = this->syndrome(received);
    Polynomial temp_1_syndrome = this->syndrome(temp_0_syndrome);
    Polynomial temp_2_error_locator = this->euclidean_error_locator(temp_1_syndrome);
    Polynomial temp_3_error_evaluator = this->euclidean_error_evaluator(temp_1_syndrome);
    vector<Element> temp_4_roots = this->chien_roots(temp_2_error_locator);
    vector<Element> temp_5_error_values = this->euclidean_error_values(temp_2_error_locator, temp_3_error_evaluator, temp_4_roots);
    Polynomial temp_6_estimated_error = this->estimated_error(temp_4_roots, temp_5_error_values);
    Polynomial temp_7_estimated_codeword = this->estimated_codeword(received, temp_6_estimated_error);
    Polynomial temp_8_estimated_message = this->estimated_message(temp_7_estimated_codeword);
    if (output_path != "output") {
        ofstream output_file("Output/" + output_path + ".txt");
        output_file << left << setw(16) << "polRec" << received << endl;
        for (int i = 0; i < temp_0_syndrome.size(); i ++) {
            output_file << left << setw(16) << "vecSyn[" + to_string(i) + "]" << temp_0_syndrome[i] << endl;
        }
        output_file << left << setw(16) << "polSyn" << temp_1_syndrome << endl;
        output_file << left << setw(16) << "polErrLoc" << temp_2_error_locator << endl;
        output_file << left << setw(16) << "polErrEva" << temp_3_error_evaluator << endl;
        for (int i = 0; i < temp_4_roots.size(); i ++) {
            output_file << left << setw(16) << "vecRoo[" + to_string(i) + "]" << temp_4_roots[i] << endl;
        }
        for (int i = 0; i < temp_5_error_values.size(); i ++) {
            output_file << left << setw(16) << "vecErrVal[" + to_string(i) + "]" << temp_5_error_values[i] << endl;
        }
        output_file << left << setw(16) << "polEstErr" << temp_6_estimated_error << endl;
        output_file << left << setw(16) << "polEstCod" << temp_7_estimated_codeword << endl;
        output_file << left << setw(16) << "polEstMes" << temp_8_estimated_message << endl;
        output_file.close();
    }
    return temp_8_estimated_message;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

vector<Element> RS::syndrome(const Polynomial& received) const {
    if (received.coefficients.size() != this->codeword_length()) {
        throw "FEC\\RS\\syndrome(const Polynomial&)\\received\\coefficients\\size";
    } else if (*received.field != *this->symbol_field) {
        throw "FEC\\RS\\syndrome(const Polynomial&)\\received\\field";
    }
    vector<Element> temp_0_elements(this->parity_length(), this->symbol_field->zero_element());
    for (int i = 0; i < temp_0_elements.size(); i ++) {
        temp_0_elements[i] = received.evaluate(this->symbol_field->general_elements[i + 1]);
    }
    return temp_0_elements;
}

Polynomial RS::syndrome(const vector<Element>& syndrome) const {
    if (syndrome.size() != 2 * this->correction_capability()) {
        throw "FEC\\RS\\syndrome(const vector<Element>&)\\syndrome\\size";
    } else if (true) {
        for (int i = 0; i < syndrome.size(); i ++) {
            if (*syndrome[i].field != *this->symbol_field) {
                throw "FEC\\RS\\syndrome(const vector<Element>&)\\syndrome\\field";
            }
        }
    }
    return Polynomial(this->symbol_field, syndrome);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

Polynomial RS::pgz_error_locator(const vector<Element>& syndrome) const {
    if (syndrome.size() != 2 * this->correction_capability()) {
        throw "FEC\\RS\\pgz_error_locator(const vector<Element>&)\\syndrome\\size";
    } else if (true) {
        for (int i = 0; i < syndrome.size(); i ++) {
            if (*syndrome[i].field != *this->symbol_field) {
                throw "FEC\\RS\\pgz_error_locator(const vector<Element>&)\\syndrome\\field";
            }
        }
    }
    vector<Element> temp_0_coefficients(1, this->symbol_field->unit_element());
    for (int i = this->correction_capability(); i > 0; i --) {
        vector<vector<Element>> temp_1_elements(i, vector<Element>(i, this->symbol_field->zero_element()));
        for (int j = 0; j < i; j ++) {
            for (int k = 0; k < i; k ++) {
                temp_1_elements[j][k] = syndrome[j + k];
            }
        }
        Matrix temp_2_matrix(this->symbol_field, temp_1_elements);
        if (temp_2_matrix.determinant() != this->symbol_field->zero_element()) {
            vector<vector<Element>> temp_3_elements(i, vector<Element>(1, this->symbol_field->zero_element()));
            for (int j = 0; j < i; j ++) {
                temp_3_elements[j][0] = -syndrome[i + j];
            }
            Matrix temp_4_matrix = (~temp_2_matrix) * Matrix(this->symbol_field, temp_3_elements);
            for (int j = 1; j <= i; j ++) {
                temp_0_coefficients.push_back(temp_4_matrix.elements[i - j][0]);
            }
            break;
        }
    }
    return Polynomial(this->symbol_field, temp_0_coefficients).redegree(this->correction_capability());
}

Polynomial RS::bm_error_locator(const vector<Element>& syndrome) const {
    if (syndrome.size() != 2 * this->correction_capability()) {
        throw "FEC\\RS\\bm_error_locator(const vector<Element>&)\\syndrome\\size";
    } else if (true) {
        for (int i = 0; i < syndrome.size(); i ++) {
            if (*syndrome[i].field != *this->symbol_field) {
                throw "FEC\\RS\\bm_error_locator(const vector<Element>&)\\syndrome\\field";
            }
        }
    }
    Polynomial temp_0_polynomial(this->symbol_field, {this->symbol_field->unit_element()});
    Polynomial temp_1_polynomial(this->symbol_field, {this->symbol_field->unit_element()});
    for (int i = 0; i < syndrome.size(); i ++) {
        Element temp_2_discrepancy = this->symbol_field->zero_element();
        for (int j = 0; j <= temp_1_polynomial.degree(); j ++) {
            temp_2_discrepancy = temp_2_discrepancy + temp_1_polynomial.coefficients[j] * syndrome[i - j];
        }
        Polynomial temp_3_polynomial = Polynomial(this->symbol_field, {this->symbol_field->zero_element(), this->symbol_field->unit_element()}) * temp_0_polynomial;
        if (temp_2_discrepancy == this->symbol_field->zero_element()) {
            temp_0_polynomial = temp_3_polynomial;
        } else {
            Polynomial temp_4_polynomial = temp_1_polynomial + Polynomial(this->symbol_field, {this->symbol_field->zero_element(), temp_2_discrepancy}) * temp_0_polynomial;
            if (2 * temp_1_polynomial.degree() <= i - 1) {
                temp_0_polynomial = temp_3_polynomial;
                temp_1_polynomial = temp_4_polynomial;
            } else {
                temp_0_polynomial = temp_1_polynomial * (~temp_2_discrepancy);
                temp_1_polynomial = temp_4_polynomial;
            }
        }
    }
    return temp_1_polynomial.redegree(this->correction_capability());
}

Polynomial RS::euclidean_error_locator(const Polynomial& syndrome) const {
    if (syndrome.degree() >= 2 * this->correction_capability()) {
        throw "FEC\\RS\\euclidean_error_locator(const Polynomial&)\\syndrome\\degree";
    } else if (*syndrome.field != *this->symbol_field) {
        throw "FEC\\RS\\euclidean_error_locator(const Polynomial&)\\syndrome\\field";
    }
    Polynomial temp_0_polynomial(this->symbol_field, vector<Element>(2 * this->correction_capability() + 1, this->symbol_field->zero_element()));
    temp_0_polynomial.coefficients[2 * this->correction_capability()] = this->symbol_field->unit_element();
    Polynomial temp_1_polynomial = syndrome;
    Polynomial temp_2_polynomial(this->symbol_field, {this->symbol_field->zero_element()});
    Polynomial temp_3_polynomial(this->symbol_field, {this->symbol_field->unit_element()});
    while (temp_1_polynomial.degree() >= this->correction_capability()) {
        temp_2_polynomial = temp_2_polynomial - (temp_0_polynomial / temp_1_polynomial) * temp_3_polynomial;
        swap(temp_2_polynomial, temp_3_polynomial);
        temp_0_polynomial = temp_0_polynomial % temp_1_polynomial;
        swap(temp_0_polynomial, temp_1_polynomial);
    }
    return temp_3_polynomial.redegree(this->correction_capability());
}

////////////////////////////////////////////////////////////////////////////////////////////////////

Polynomial RS::pgz_error_evaluator(const Polynomial& syndrome, const Polynomial& error_locator) const {
    if (syndrome.degree() >= 2 * this->correction_capability()) {
        throw "FEC\\RS\\pgz_error_evaluator(const Polynomial&, const Polynomial&)\\syndrome\\degree";
    } else if (error_locator.degree() > this->correction_capability()) {
        throw "FEC\\RS\\pgz_error_evaluator(const Polynomial&, const Polynomial&)\\error_locator\\degree";
    } else if (*syndrome.field != *this->symbol_field) {
        throw "FEC\\RS\\pgz_error_evaluator(const Polynomial&, const Polynomial&)\\syndrome\\field";
    } else if (*error_locator.field != *this->symbol_field) {
        throw "FEC\\RS\\pgz_error_evaluator(const Polynomial&, const Polynomial&)\\error_locator\\field";
    }
    vector<Element> temp_0_coefficients(2 * this->correction_capability() + 1, this->symbol_field->zero_element());
    temp_0_coefficients[2 * this->correction_capability()] = this->symbol_field->unit_element();
    Polynomial temp_1_polynomial(this->symbol_field, temp_0_coefficients);
    return ((syndrome * error_locator) % temp_1_polynomial).redegree(this->correction_capability() - 1);
}

Polynomial RS::bm_error_evaluator(const Polynomial& syndrome, const Polynomial& error_locator) const {
    if (syndrome.degree() >= 2 * this->correction_capability()) {
        throw "FEC\\RS\\bm_error_evaluator(const Polynomial&, const Polynomial&)\\syndrome\\degree";
    } else if (error_locator.degree() > this->correction_capability()) {
        throw "FEC\\RS\\bm_error_evaluator(const Polynomial&, const Polynomial&)\\error_locator\\degree";
    } else if (*syndrome.field != *this->symbol_field) {
        throw "FEC\\RS\\bm_error_evaluator(const Polynomial&, const Polynomial&)\\syndrome\\field";
    } else if (*error_locator.field != *this->symbol_field) {
        throw "FEC\\RS\\bm_error_evaluator(const Polynomial&, const Polynomial&)\\error_locator\\field";
    }
    return this->pgz_error_evaluator(syndrome, error_locator);
}

Polynomial RS::euclidean_error_evaluator(const Polynomial& syndrome) const {
    if (syndrome.degree() >= 2 * this->correction_capability()) {
        throw "FEC\\RS\\euclidean_error_evaluator(const Polynomial&)\\syndrome\\degree";
    } else if (*syndrome.field != *this->symbol_field) {
        throw "FEC\\RS\\euclidean_error_evaluator(const Polynomial&)\\syndrome\\field";
    }
    Polynomial temp_0_polynomial(this->symbol_field, vector<Element>(2 * this->correction_capability() + 1, this->symbol_field->zero_element()));
    temp_0_polynomial.coefficients[2 * this->correction_capability()] = this->symbol_field->unit_element();
    Polynomial temp_1_polynomial = syndrome;
    while (temp_1_polynomial.degree() >= this->correction_capability()) {
        temp_0_polynomial = temp_0_polynomial % temp_1_polynomial;
        swap(temp_0_polynomial, temp_1_polynomial);
    }
    return temp_1_polynomial.redegree(this->correction_capability() - 1);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

vector<Element> RS::chien_roots(const Polynomial& error_locator) const {
    if (error_locator.degree() > this->correction_capability()) {
        throw "FEC\\RS\\chien_roots(const Polynomial&)\\error_locator\\degree";
    } else if (*error_locator.field != *this->symbol_field) {
        throw "FEC\\RS\\chien_roots(const Polynomial&)\\error_locator\\field";
    }
    vector<Element> temp_0_elements;
    for (int i = 0; i < this->symbol_field->size() - 1; i ++) {
        if (error_locator.evaluate(this->symbol_field->general_elements[i]) == this->symbol_field->zero_element()) {
            temp_0_elements.push_back(this->symbol_field->general_elements[i]);
        }
    }
    return temp_0_elements;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

vector<Element> RS::pgz_error_values(const Polynomial& error_locator, const Polynomial& error_evaluator, const vector<Element>& roots) const {
    if (error_locator.degree() > this->correction_capability()) {
        throw "FEC\\RS\\pgz_error_values(const Polynomial&, const Polynomial&, const vector<Element>&)\\error_locator\\degree";
    } else if (error_evaluator.degree() >= 2 * this->correction_capability()) {
        throw "FEC\\RS\\pgz_error_values(const Polynomial&, const Polynomial&, const vector<Element>&)\\error_evaluator\\degree";
    } else if (roots.size() > this->correction_capability()) {
        throw "FEC\\RS\\pgz_error_values(const Polynomial&, const Polynomial&, const vector<Element>&)\\roots\\size";
    } else if (*error_locator.field != *this->symbol_field) {
        throw "FEC\\RS\\pgz_error_values(const Polynomial&, const Polynomial&, const vector<Element>&)\\error_locator\\field";
    } else if (*error_evaluator.field != *this->symbol_field) {
        throw "FEC\\RS\\pgz_error_values(const Polynomial&, const Polynomial&, const vector<Element>&)\\error_evaluator\\field";
    } else if (true) {
        for (int i = 0; i < roots.size(); i ++) {
            if (*roots[i].field != *this->symbol_field) {
                throw "FEC\\RS\\pgz_error_values(const Polynomial&, const Polynomial&, const vector<Element>&)\\roots\\field";
            }
        }
    }
    vector<Element> temp_0_elements;
    for (int i = 0; i < roots.size(); i ++) {
        temp_0_elements.push_back(error_evaluator.evaluate(roots[i]) / error_locator.derivative().evaluate(roots[i]));
    }
    return temp_0_elements;
}

vector<Element> RS::bm_error_values(const Polynomial& error_locator, const Polynomial& error_evaluator, const vector<Element>& roots) const {
    if (error_locator.degree() > this->correction_capability()) {
        throw "FEC\\RS\\bm_error_values(const Polynomial&, const Polynomial&, const vector<Element>&)\\error_locator\\degree";
    } else if (error_evaluator.degree() >= 2 * this->correction_capability()) {
        throw "FEC\\RS\\bm_error_values(const Polynomial&, const Polynomial&, const vector<Element>&)\\error_evaluator\\degree";
    } else if (roots.size() > this->correction_capability()) {
        throw "FEC\\RS\\bm_error_values(const Polynomial&, const Polynomial&, const vector<Element>&)\\roots\\size";
    } else if (*error_locator.field != *this->symbol_field) {
        throw "FEC\\RS\\bm_error_values(const Polynomial&, const Polynomial&, const vector<Element>&)\\error_locator\\field";
    } else if (*error_evaluator.field != *this->symbol_field) {
        throw "FEC\\RS\\bm_error_values(const Polynomial&, const Polynomial&, const vector<Element>&)\\error_evaluator\\field";
    } else if (true) {
        for (int i = 0; i < roots.size(); i ++) {
            if (*roots[i].field != *this->symbol_field) {
                throw "FEC\\RS\\bm_error_values(const Polynomial&, const Polynomial&, const vector<Element>&)\\roots\\field";
            }
        }
    }
    return this->pgz_error_values(error_locator, error_evaluator, roots);
}

vector<Element> RS::euclidean_error_values(const Polynomial& error_locator, const Polynomial& error_evaluator, const vector<Element>& roots) const {
    if (error_locator.degree() > this->correction_capability()) {
        throw "FEC\\RS\\euclidean_error_values(const Polynomial&, const Polynomial&, const vector<Element>&)\\error_locator\\degree";
    } else if (error_evaluator.degree() >= 2 * this->correction_capability()) {
        throw "FEC\\RS\\euclidean_error_values(const Polynomial&, const Polynomial&, const vector<Element>&)\\error_evaluator\\degree";
    } else if (roots.size() > this->correction_capability()) {
        throw "FEC\\RS\\euclidean_error_values(const Polynomial&, const Polynomial&, const vector<Element>&)\\roots\\size";
    } else if (*error_locator.field != *this->symbol_field) {
        throw "FEC\\RS\\euclidean_error_values(const Polynomial&, const Polynomial&, const vector<Element>&)\\error_locator\\field";
    } else if (*error_evaluator.field != *this->symbol_field) {
        throw "FEC\\RS\\euclidean_error_values(const Polynomial&, const Polynomial&, const vector<Element>&)\\error_evaluator\\field";
    } else if (true) {
        for (int i = 0; i < roots.size(); i ++) {
            if (*roots[i].field != *this->symbol_field) {
                throw "FEC\\RS\\euclidean_error_values(const Polynomial&, const Polynomial&, const vector<Element>&)\\roots\\field";
            }
        }
    }
    vector<Element> temp_0_elements;
    for (int i = 0; i < roots.size(); i ++) {
        temp_0_elements.push_back(-error_evaluator.evaluate(roots[i]) / error_locator.derivative().evaluate(roots[i]));
    }
    return temp_0_elements;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

Polynomial RS::estimated_error(const vector<Element>& roots, const vector<Element>& error_values) const {
    if (roots.size() > this->correction_capability() || roots.size() != error_values.size()) {
        throw "FEC\\RS\\estimated_error(const vector<Element>&, const vector<Element>&)\\roots\\size";
    } else if (error_values.size() > this->correction_capability() || error_values.size() != roots.size()) {
        throw "FEC\\RS\\estimated_error(const vector<Element>&, const vector<Element>&)\\error_values\\size";
    } else if (true) {
        for (int i = 0; i < roots.size(); i ++) {
            if (*roots[i].field != *this->symbol_field) {
                throw "FEC\\RS\\estimated_error(const vector<Element>&, const vector<Element>&)\\roots\\field";
            }
        }
    } else if (true) {
        for (int i = 0; i < error_values.size(); i ++) {
            if (*error_values[i].field != *this->symbol_field) {
                throw "FEC\\RS\\estimated_error(const vector<Element>&, const vector<Element>&)\\error_values\\field";
            }
        }
    }
    vector<Element> temp_0_coefficients(this->codeword_length(), this->symbol_field->zero_element());
    for (int i = 0; i < roots.size(); i ++) {
        Element temp_1_element = ~roots[i];
        for (int j = 0; j < this->symbol_field->size() - 1; j ++) {
            if (temp_1_element == this->symbol_field->general_elements[j]) {
                temp_0_coefficients[j] = error_values[i];
                break;
            }
            if (j == this->symbol_field->size() - 2) {
                throw "FEC\\RS\\estimated_error(const vector<Element>&, const vector<Element>&)";
            }
        }
    }
    return Polynomial(this->symbol_field, temp_0_coefficients);
}

Polynomial RS::estimated_codeword(const Polynomial& received, const Polynomial& estimated_error) const {
    if (received.coefficients.size() != this->codeword_length()) {
        throw "FEC\\RS\\estimated_codeword(const Polynomial&, const Polynomial&)\\received\\coefficients\\size";
    } else if (estimated_error.coefficients.size() != this->codeword_length()) {
        throw "FEC\\RS\\estimated_codeword(const Polynomial&, const Polynomial&)\\estimated_error\\coefficients\\size";
    } else if (*received.field != *this->symbol_field) {
        throw "FEC\\RS\\estimated_codeword(const Polynomial&, const Polynomial&)\\received\\field";
    } else if (*estimated_error.field != *this->symbol_field) {
        throw "FEC\\RS\\estimated_codeword(const Polynomial&, const Polynomial&)\\estimated_error\\field";
    }
    return (received - estimated_error).redegree(this->codeword_length() - 1);
}

Polynomial RS::estimated_message(const Polynomial& estimated_codeword) const {
    if (estimated_codeword.coefficients.size() != this->codeword_length()) {
        throw "FEC\\RS\\estimated_message(const Polynomial&)\\estimated_codeword\\coefficients\\size";
    } else if (*estimated_codeword.field != *this->symbol_field) {
        throw "FEC\\RS\\estimated_message(const Polynomial&)\\estimated_codeword\\field";
    }
    return estimated_codeword.redegree(this->message_length() - 1);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////