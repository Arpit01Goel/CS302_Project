// gtest_simplex_tests.cpp (reads input_#.txt files)
// Google Test file for Simplex implementation that reads LPs from input files and
// constructs canonical (A,B,C) then calls the Simplex class directly.
// -------------------------
// Usage notes:
// 1) Place this file next to your implementation.cpp (which must define Simplex and
//    provide: Simplex(const vector<vector<float>>&, const vector<float>&, const vector<float>&),
//    void CalculateSimplex(), float getMaximum() const ).
// 2) The test will open files named: tests/input_1.txt .. tests/input_10.txt (relative path).
//    If you put files elsewhere, update the BASE_PATH constant below.
// 3) Build with gtest (example):
//    g++ -std=c++11 gtest_simplex_tests.cpp implementation.cpp -lgtest -lpthread -O2 -o test_runner
//    ./test_runner
// -------------------------

#include <gtest/gtest.h>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>

#include "implementation.cpp"  // or change to "simplex.h" if you split files

using std::vector;
using std::string;

static const string BASE_PATH = "./tests/"; // where input files are located
static const double EPS = 1e-3;


// Parse the input file (format defined by assistant):
// Line1: max|min
// Line2: m n
// Line3: c1 c2 ... cn
// Next m lines: a1 a2 ... an REL b   where REL is <= or >=
// This function builds canonical A (with m slacks appended to n columns), B, C (internal sign conv)
static bool parse_input_file(const string &filename, vector<vector<float>> &A_out, vector<float> &B_out, vector<float> &C_out) {
    std::ifstream in(filename);
    if (!in.is_open()) return false;

    string token;
    if (!(in >> token)) return false;
    string opt = toLower(token);
    if (opt != "max" && opt != "min") return false;

    int m,n;
    if (!(in >> m >> n)) return false;

    vector<double> obj(n);
    for (int i=0;i<n;++i) if (!(in >> obj[i])) return false;

    string line;
    getline(in, line); // consume rest of line

    int totalCols = n + m;
    vector<vector<float>> A(m, vector<float>(totalCols, 0.0f));
    vector<float> B(m, 0.0f);

    for (int i=0;i<m;++i) {
        if (!std::getline(in, line)) return false;
        if (line.find_first_not_of(" ") == string::npos) { --i; continue; }
        std::stringstream ss(line);
        vector<double> coeffs(n);
        for (int j=0;j<n;++j) if (!(ss >> coeffs[j])) return false;
        string rel;
        if (!(ss >> rel)) return false;
        double rhs;
        if (!(ss >> rhs)) return false;
        if (rel == ">=") {
            for (int j=0;j<n;++j) coeffs[j] = -coeffs[j];
            rhs = -rhs;
        } else if (rel != "<=") {
            return false;
        }
        for (int j=0;j<n;++j) A[i][j] = static_cast<float>(coeffs[j]);
        A[i][n + i] = 1.0f; // slack
        B[i] = static_cast<float>(rhs);
    }

    vector<float> C(totalCols, 0.0f);
    if (opt == "max") {
        for (int i=0;i<n;++i) C[i] = static_cast<float>(-obj[i]);
    } else {
        for (int i=0;i<n;++i) C[i] = static_cast<float>(obj[i]);
    }

    A_out = std::move(A);
    B_out = std::move(B);
    C_out = std::move(C);
    return true;
}

// helper to run a single test case file and return the computed maximum
static float run_case_file_and_get_max(const string &filepath) {
    vector<vector<float>> A;
    vector<float> B,C;
    bool ok = parse_input_file(filepath, A, B, C);
    if (!ok) {
        std::cerr << "Failed to parse: " << filepath << " ";
        return NAN;
    }
    Simplex s(A,B,C);
    s.CalculateSimplex();
    return s.getMaximum();
}

// The 10 tests will read files tests/input_1.txt ... tests/input_10.txt
TEST(SimplexInputFiles, Input1) { EXPECT_NEAR(run_case_file_and_get_max(BASE_PATH + "input_1.txt"), 10.0, EPS); }
TEST(SimplexInputFiles, Input2) { EXPECT_NEAR(run_case_file_and_get_max(BASE_PATH + "input_2.txt"), 16.0, EPS); }
TEST(SimplexInputFiles, Input3) { EXPECT_NEAR(run_case_file_and_get_max(BASE_PATH + "input_3.txt"), 21.0, EPS); }
TEST(SimplexInputFiles, Input4) { EXPECT_NEAR(run_case_file_and_get_max(BASE_PATH + "input_4.txt"), 20.0, EPS); }
TEST(SimplexInputFiles, Input5) { EXPECT_NEAR(run_case_file_and_get_max(BASE_PATH + "input_5.txt"), 19.0, EPS); }
TEST(SimplexInputFiles, Input6) { EXPECT_NEAR(run_case_file_and_get_max(BASE_PATH + "input_6.txt"), 18.0, EPS); }
TEST(SimplexInputFiles, Input7) { EXPECT_NEAR(run_case_file_and_get_max(BASE_PATH + "input_7.txt"), 17.0, EPS); }
TEST(SimplexInputFiles, Input8) { EXPECT_NEAR(run_case_file_and_get_max(BASE_PATH + "input_8.txt"), 7.0, EPS); }
TEST(SimplexInputFiles, Input9) { EXPECT_NEAR(run_case_file_and_get_max(BASE_PATH + "input_9.txt"), 302.0, EPS); }
TEST(SimplexInputFiles, Input10) { EXPECT_NEAR(run_case_file_and_get_max(BASE_PATH + "input_10.txt"), 10.0, EPS); }

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
