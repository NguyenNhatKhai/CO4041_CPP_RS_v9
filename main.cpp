////////////////////////////////////////////////////////////////////////////////////////////////////
//
// File: main.cpp
// Author: Nhat Khai Nguyen
//
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "main.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

int main() {
    cout << "Harw" << endl;
    return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Example 6.7 on page 239

// #define BIN0 Element(false)
// #define BIN1 Element(true)
// int INT0 = 3;

// Polynomial POL0(&fields::default_field, {BIN1, BIN1, BIN0, BIN0, BIN1});
// Field FIE0(&POL0);
// RS RS0(&FIE0, INT0);

// #define ELE(i) FIE0.general_elements[i]

// int main() {
//     try {        
//         Polynomial polMes(&FIE0, vector<Element>(RS0.message_length(), FIE0.zero_element()));
//         Polynomial polErr(&FIE0, vector<Element>(RS0.codeword_length(), FIE0.zero_element()));
//         polErr.coefficients[3] = ELE(7); polErr.coefficients[6] = ELE(3); polErr.coefficients[12] = ELE(4);

//         Polynomial polCod = RS0.systematic_encode(polMes, "rs_encode");
//         Polynomial polRec = RS0.add_error(polCod, polErr, "rs_add_error");
//         Polynomial polEstMes = RS0.bm_decode(polRec, "rs_decode");

//         cout << "n" << RS0.codeword_length() << ", k" << RS0.message_length() << ", p" << RS0.parity_length() << ", s" << RS0.symbol_size() << ", d" << RS0.detection_capability() << ", c" << RS0.correction_capability() << ", g" << RS0.generator_polynomial << endl;
//         cout << boolalpha << (polMes == polEstMes) << endl;
//     } catch (const char* error_message) {
//         cout << error_message << endl;
//     }
//     return 0;
// }

////////////////////////////////////////////////////////////////////////////////////////////////////
// Confirmed that encoding can be performed as matrix multiplication

// #define BIN0 Element(false)
// #define BIN1 Element(true)
// int INT0 = 1;

// Polynomial POL0(&fields::default_field, {BIN1, BIN1, BIN0, BIN1});
// Field FIE0(&POL0);
// RS RS0(&FIE0, INT0);

// #define ELE(i) FIE0.general_elements[i]

// int main() {
//     try {
//         // ofstream output_file("Output/output.txt");
//         cout << left << setw(16) << "polGen" << RS0.generator_polynomial << endl;

//         vector<vector<Element>> temp_0_elements(RS0.message_length(), vector<Element>(RS0.codeword_length(), FIE0.zero_element()));
//         for (int i = 0; i < RS0.message_length(); i ++) {
//             Polynomial polMes(&FIE0, vector<Element>(RS0.message_length(), FIE0.zero_element()));
//             polMes.coefficients[i] = FIE0.unit_element();
//             Polynomial polCod = RS0.systematic_encode(polMes);
//             temp_0_elements[i] = polCod.coefficients;
//         }
//         Matrix temp_0_matrix = Matrix(&FIE0, temp_0_elements);

//         Polynomial polMes(&FIE0, vector<Element>(RS0.message_length(), FIE0.zero_element()));
//         polMes.coefficients[0] = ELE(0);
//         polMes.coefficients[1] = ELE(0);
//         polMes.coefficients[2] = ELE(1);
//         cout << left << setw(16) << "polMes" << polMes << endl;
//         vector<vector<Element>> temp_1_elements(1, vector<Element>(RS0.codeword_length(), FIE0.zero_element()));
//         temp_1_elements[0] = polMes.coefficients;
//         Matrix temp_1_matrix = Matrix(&FIE0, temp_1_elements);

//         cout << left << setw(16) << "polCodSys" << RS0.systematic_encode(polMes) << endl;
//         cout << left << setw(16) << "polCodMat" << Polynomial(&FIE0, (temp_1_matrix * temp_0_matrix).elements[0]) << endl;
        
//         // output_file.close();
//         cout << "Harw" << endl;
//     } catch (const char* error_message) {
//         cout << error_message << endl;
//     }
//     return 0;
// }

////////////////////////////////////////////////////////////////////////////////////////////////////
// This main function can be used to get all possible finite field over GF(2^n^m)

// #define BIN0 Element(false)
// #define BIN1 Element(true)
// int INT0 = 1;

// Polynomial POL0(&fields::default_field, {BIN1, BIN1, BIN0, BIN1});
// Field FIE0(&POL0);
// RS RS0(&FIE0, INT0);

// #define ELE(i) FIE0.general_elements[i]

// int main() {
//     for (int i = 0; i < FIE0.size(); i ++) {
//         for (int j = 0; j < FIE0.size(); j ++) {
//             try {
//                 cout << "Running:" << " i" << i << " j" << j << " ";
//                 Polynomial polExtGen(&FIE0, vector<Element>(RS0.parity_length() + 1, FIE0.zero_element()));
//                 polExtGen.coefficients[0] = ELE(i);
//                 polExtGen.coefficients[1] = ELE(j);
//                 polExtGen.coefficients[2] = ELE(0);
//                 Field fieExt(&polExtGen);
//                 cout << "g" << polExtGen << endl;
//                 ofstream output_file("Output/output_" + to_string(i) + "_" + to_string(j) + ".txt");
//                 for (int k = 0; k < fieExt.size(); k ++) {
//                     output_file << fieExt.general_elements[k] << endl;
//                 }
//                 output_file.close();
//             } catch (const char* error_message) {
//                 cout << error_message << endl;
//             }
//         }
//     }
//     return 0;
// }

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////