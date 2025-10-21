// testing_2.cpp
#include <gtest/gtest.h>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <cmath>

static const std::string EXECUTABLE = "./simplex_clrs";      // path to solver executable
static const std::string BASE_PATH  = "test_cases_2/";      // directory with input files
static const double EPS = 1e-3;                             // tolerance for objective comparisons

// Expected result struct
struct Expected {
    std::string status; // "OPTIMAL", "UNBOUNDED", "INFEASIBLE"
    double value;       // objective value if status == OPTIMAL (NaN otherwise)
};

// The 10 expected results for these edge tests
static const Expected EXPECTEDS[10] = {
    { "UNBOUNDED", NAN },  // input_1
    { "INFEASIBLE", NAN }, // input_2
    { "OPTIMAL", 0.0 },    // input_3
    { "OPTIMAL", 2.0 },    // input_4
    { "OPTIMAL", 10.0 },   // input_5
    { "OPTIMAL", 10.0 },   // input_6
    { "OPTIMAL", 0.0 },    // input_7
    { "UNBOUNDED", NAN },  // input_8
    { "INFEASIBLE", NAN }, // input_9
    { "OPTIMAL", 20.0 }    // input_10
};

// Run the solver on inputFile and return stdout as string (captures stderr too).
static std::string runSolver(const std::string &inputFile) {
    std::string tmpOut = "/tmp/solver_test_out.txt";
    std::string cmd = EXECUTABLE + " " + inputFile + " > " + tmpOut + " 2>&1";
    int rc = std::system(cmd.c_str());
    (void)rc; // we don't require rc==0; solver may return non-zero but still print
    std::ifstream in(tmpOut);
    if (!in.is_open()) return std::string();
    std::ostringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

// Parse solver output: find "INFEASIBLE" or "UNBOUNDED" or "OPTIMAL" and "Objective = <value>"
static std::pair<std::string,double> parseOutput(const std::string &out) {
    std::istringstream iss(out);
    std::string line;
    std::string status;
    double objective = NAN;
    while (std::getline(iss, line)) {
        // trim leading whitespace
        size_t p = line.find_first_not_of(" \t\r\n");
        if (p != std::string::npos) line = line.substr(p);
        if (line.rfind("INFEASIBLE", 0) == 0) return {"INFEASIBLE", NAN};
        if (line.rfind("UNBOUNDED", 0) == 0) return {"UNBOUNDED", NAN};
        if (line.rfind("OPTIMAL", 0) == 0) { status = "OPTIMAL"; continue; }
        auto pos = line.find("Objective");
        if (pos != std::string::npos) {
            auto eq = line.find('=', pos);
            if (eq != std::string::npos) {
                std::string num = line.substr(eq + 1);
                // trim
                size_t s = num.find_first_not_of(" \t");
                if (s != std::string::npos) num = num.substr(s);
                size_t e = num.find_last_not_of(" \t\r\n");
                if (e != std::string::npos) num = num.substr(0, e+1);
                try {
                    objective = std::stod(num);
                } catch (...) {
                    objective = NAN;
                }
            }
        }
    }
    if (!status.empty()) return {status, objective};
    if (!std::isnan(objective)) return {"OPTIMAL", objective};
    return {"", NAN};
}

class EdgeCaseTests : public ::testing::Test {
protected:
    void runCase(int idx) {
        ASSERT_GE(idx, 1);
        ASSERT_LE(idx, 10);
        std::string infile = BASE_PATH + "input_" + std::to_string(idx) + ".txt";
        std::string out = runSolver(infile);
        ASSERT_FALSE(out.empty()) << "No output captured (solver not found or produced no output).";
        auto parsed = parseOutput(out);
        const Expected &exp = EXPECTEDS[idx-1];

        if (exp.status == "OPTIMAL") {
            ASSERT_EQ(parsed.first, "OPTIMAL") << "Expected OPTIMAL for test " << idx << " but got '" << parsed.first << "'\nSolver output:\n" << out;
            ASSERT_FALSE(std::isnan(parsed.second)) << "Solver did not print an objective for test " << idx << "\nOutput:\n" << out;
            EXPECT_NEAR(parsed.second, exp.value, EPS) << "Objective mismatch for test " << idx << "\nOutput:\n" << out;
        } else {
            EXPECT_EQ(parsed.first, exp.status) << "Expected status " << exp.status << " for test " << idx << " but got '" << parsed.first << "'\nOutput:\n" << out;
        }
    }
};

TEST_F(EdgeCaseTests, Case1) { runCase(1); }
TEST_F(EdgeCaseTests, Case2) { runCase(2); }
TEST_F(EdgeCaseTests, Case3) { runCase(3); }
TEST_F(EdgeCaseTests, Case4) { runCase(4); }
TEST_F(EdgeCaseTests, Case5) { runCase(5); }
TEST_F(EdgeCaseTests, Case6) { runCase(6); }
TEST_F(EdgeCaseTests, Case7) { runCase(7); }
TEST_F(EdgeCaseTests, Case8) { runCase(8); }
TEST_F(EdgeCaseTests, Case9) { runCase(9); }
TEST_F(EdgeCaseTests, Case10){ runCase(10); }

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
