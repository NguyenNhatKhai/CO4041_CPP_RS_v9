////////////////////////////////////////////////////////////////////////////////////////////////////
//
// File: fec.h
// Author: Nhat Khai Nguyen
//
////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef _FEC_H_
#define _FEC_H_

#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "../FFA/ffa.h"
#include "../Maths/maths.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

class RS;

////////////////////////////////////////////////////////////////////////////////////////////////////

class RS {
    public:
    Field* symbol_field;
    Polynomial generator_polynomial;

    public:
    RS() = delete;
    RS(Field* symbol_field, int correction_capability);
    ~RS() = default;

    public:
    int codeword_length() const;
    int message_length() const;
    int parity_length() const;
    int symbol_size() const;
    int detection_capability() const;
    int correction_capability() const;

    public:
    Polynomial systematic_encode(const Polynomial& message, const string& output_path = "output") const;
    Polynomial nonsystematic_encode(const Polynomial& message, const string& output_path = "output") const;

    public:
    Polynomial add_error(const Polynomial& codeword, const Polynomial& error, const string& output_path = "output") const;

    public:
    Polynomial pgz_decode(const Polynomial& received, const string& output_path = "output") const;
    Polynomial bm_decode(const Polynomial& received, const string& output_path = "output") const;
    Polynomial euclidean_decode(const Polynomial& received, const string& output_path = "output") const;

    private:
    vector<Element> syndrome(const Polynomial& received) const;
    Polynomial syndrome(const vector<Element>& syndrome) const;

    private:
    Polynomial pgz_error_locator(const vector<Element>& syndrome) const;
    Polynomial bm_error_locator(const vector<Element>& syndrome) const;
    Polynomial euclidean_error_locator(const Polynomial& syndrome) const;

    private:
    Polynomial pgz_error_evaluator(const Polynomial& syndrome, const Polynomial& error_locator) const;
    Polynomial bm_error_evaluator(const Polynomial& syndrome, const Polynomial& error_locator) const;
    Polynomial euclidean_error_evaluator(const Polynomial& syndrome) const;

    private:
    vector<Element> chien_roots(const Polynomial& error_locator) const;

    private:
    vector<Element> pgz_error_values(const Polynomial& error_locator, const Polynomial& error_evaluator, const vector<Element>& roots) const;
    vector<Element> bm_error_values(const Polynomial& error_locator, const Polynomial& error_evaluator, const vector<Element>& roots) const;
    vector<Element> euclidean_error_values(const Polynomial& error_locator, const Polynomial& error_evaluator, const vector<Element>& roots) const;

    private:
    Polynomial estimated_error(const vector<Element>& roots, const vector<Element>& error_values) const;
    Polynomial estimated_codeword(const Polynomial& received, const Polynomial& estimated_error) const;
    Polynomial estimated_message(const Polynomial& estimated_codeword) const;
};

////////////////////////////////////////////////////////////////////////////////////////////////////

#include "rs.cpp"

#endif

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////