// gtest_simplex_clrs.cpp
#include <gtest/gtest.h>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <iostream>
#include <cmath>

static const std::string EXECUTABLE = "./simplex_clrs"; // path to compiled solver
static const std::string BASE_PATH = "./tests/";        // path where input_1.txt ... input_10.txt are located
static const double EPS = 1e-3;

// Run `exec path` with inputFile as argument, redirect stdout to a temp file, return whole stdout as string.
// Returns empty string if execution failed.
static std::string runSolverAndCaptureOutput(const std::string &execPath, const std::string &inputFile) {
    // temporary output file (unique-ish)
    std::string tmpOut = "/tmp/solver_out.txt";
    std::string cmd = execPath + " " + inputFile + " > " + tmpOut + " 2>&1";
    int ret = std::system(cmd.c_str());
    if (ret != 0 && ret != 256 /* some systems return 256 for exit 1 */) {
        // still attempt to read file; test may print diagnostics even with non-zero exit
    }
    std::ifstream ifs(tmpOut);
    if (!ifs.is_open()) return std::string();
    std::stringstream ss;
    ss << ifs.rdbuf();
    return ss.str();
}

// Extract objective value from solver stdout. If the solver printed INFEASIBLE or UNBOUNDED,
// return pair(statusString, numericValue). For INF/UNB return status string and NAN as numeric.
static std::pair<std::string,double> parseSolverOutput(const std::string &out) {
    std::istringstream iss(out);
    std::string line;
    std::string status;
    double objective = NAN;
    while (std::getline(iss, line)) {
        // trim leading/trailing spaces
        auto l = line.find_first_not_of(" \t\r\n");
        if (l == std::string::npos) continue;
        line = line.substr(l);
        // detect status words
        if (line.rfind("INFEASIBLE", 0) == 0) {
            return {"INFEASIBLE", NAN};
        }
        if (line.rfind("UNBOUNDED", 0) == 0) {
            return {"UNBOUNDED", NAN};
        }
        if (line.rfind("OPTIMAL", 0) == 0) {
            status = "OPTIMAL";
            continue;
        }
        // objective line often like: "Objective = 10" or "Objective = 10.000"
        if (line.find("Objective") != std::string::npos) {
            // find '=' in line
            auto pos = line.find('=');
            if (pos != std::string::npos) {
                std::string numstr = line.substr(pos+1);
                // trim
                auto s = numstr.find_first_not_of(" \t");
                if (s != std::string::npos) numstr = numstr.substr(s);
                auto e = numstr.find_last_not_of(" \t\r\n");
                if (e != std::string::npos) numstr = numstr.substr(0, e+1);
                try {
                    objective = std::stod(numstr);
                } catch (...) {
                    objective = NAN;
                }
            }
        }
    }
    if (!status.empty()) return {status, objective};
    // fallback: if the program printed a line like "Objective =" without OPTIMAL, still return it
    if (!std::isnan(objective)) return {"OPTIMAL", objective};
    return {"", NAN};
}

// Expected results for the 10 test cases (from earlier test set).
// If a case is infeasible or unbounded set expected_value to NAN and expected_status accordingly.
struct Expected {
    std::string status; // "OPTIMAL", "INFEASIBLE", "UNBOUNDED"
    double value;       // objective when status == "OPTIMAL"
};

// Fill this with your expected outputs for input_1..input_10.
// These were the expected objective values earlier; adjust if your desired expected results differ.
static const Expected EXPECTEDS[10] = {
    { "OPTIMAL", 10.0 },   // input_1.txt
    { "OPTIMAL", 16.0 },   // input_2.txt
    { "OPTIMAL", 21.0 },   // input_3.txt
    { "OPTIMAL", 20.0 },   // input_4.txt
    { "OPTIMAL", 19.0 },   // input_5.txt
    { "OPTIMAL", 18.0 },   // input_6.txt
    { "OPTIMAL", 17.0 },   // input_7.txt
    { "OPTIMAL", 7.0  },   // input_8.txt
    { "OPTIMAL", 302.0 },  // input_9.txt
    { "OPTIMAL", 10.0 }    // input_10.txt
};

class SimplexCLRSTest : public ::testing::Test {
protected:
    // Helper to run a single numbered test (1-based)
    void runCase(int idx) {
        ASSERT_GE(idx, 1);
        ASSERT_LE(idx, 10);
        std::string infile = BASE_PATH + "input_" + std::to_string(idx) + ".txt";
        std::string out = runSolverAndCaptureOutput(EXECUTABLE, infile);
        ASSERT_FALSE(out.empty()) << "No output captured for test case " << idx << " (cmd failed or produced no output)";
        auto parsed = parseSolverOutput(out);
        const Expected &exp = EXPECTEDS[idx-1];

        if (exp.status == "OPTIMAL") {
            ASSERT_EQ(parsed.first, "OPTIMAL") << "Expected OPTIMAL for test " << idx << " but solver reported '" << parsed.first << "'\nSolver output:\n" << out;
            ASSERT_FALSE(std::isnan(parsed.second)) << "Solver did not print objective for test " << idx << "\nOutput:\n" << out;
            EXPECT_NEAR(parsed.second, exp.value, EPS) << "Objective mismatch for test " << idx << "\nSolver output:\n" << out;
        } else {
            // expect infeasible or unbounded
            EXPECT_EQ(parsed.first, exp.status) << "Expected status " << exp.status << " for test " << idx << " but got '" << parsed.first << "'.\nOutput:\n" << out;
        }
    }
};

TEST_F(SimplexCLRSTest, Input1) { runCase(1); }
TEST_F(SimplexCLRSTest, Input2) { runCase(2); }
TEST_F(SimplexCLRSTest, Input3) { runCase(3); }
TEST_F(SimplexCLRSTest, Input4) { runCase(4); }
TEST_F(SimplexCLRSTest, Input5) { runCase(5); }
TEST_F(SimplexCLRSTest, Input6) { runCase(6); }
TEST_F(SimplexCLRSTest, Input7) { runCase(7); }
TEST_F(SimplexCLRSTest, Input8) { runCase(8); }
TEST_F(SimplexCLRSTest, Input9) { runCase(9); }
TEST_F(SimplexCLRSTest, Input10){ runCase(10); }

// main for gtest
int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
