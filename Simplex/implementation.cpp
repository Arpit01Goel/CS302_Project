#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <limits>
using namespace std;


class Simplex
{
private:
    int rows, cols;
    // Stores coefficients of all the variables
    vector<vector<float>> A;
    // Stores constants of constraints
    vector<float> B;
    // Stores the coefficients of the objective function
    vector<float> C;

    float maximum;

    bool isUnbounded;
    vector<float> solution;

public:
    Simplex(vector<vector<float>> matrix, vector<float> b, vector<float> c)
    {
        maximum = 0;
        isUnbounded = false;
        rows = (int)matrix.size();
        cols = (rows > 0 ? (int)matrix[0].size() : 0);

        // initialize A, B, C
        A.assign(rows, vector<float>(cols, 0.0f));
        B.assign(b.begin(), b.end());
        C.assign(c.begin(), c.end());

        // copy matrix safely (in case provided matrix has consistent dimensions)
        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < cols; j++)
            {
                A[i][j] = matrix[i][j];
            }
        }

        // ensure B size matches rows
        if ((int)B.size() != rows)
        {
            B.assign(rows, 0.0f);
        }
        // ensure C size matches cols
        if ((int)C.size() != cols)
        {
            C.assign(cols, 0.0f);
        }
    }

    bool simplexAlgorithmCalculataion()
    {
        // Check whether the table is optimal, if optimal no need to process further
        if (checkOptimality() == true)
        {
            return true;
        }

        // Find the column which has the pivot. The least coefficient of the objective function(C array).
        int pivotColumn = findPivotColumn();

        // If no negative reduced cost found -> optimal
        if (pivotColumn == -1)
        {
            return true;
        }

        if (isUnbounded == true)
        {
            cout << "Error unbounded" << endl;
            return true;
        }

        // Find the row with the pivot value. The least value item's row in the B array
        int pivotRow = findPivotRow(pivotColumn);

        // If no valid pivot row (unbounded), stop
        if (pivotRow == -1)
        {
            cout << "Error unbounded or no valid pivot row" << endl;
            isUnbounded = true;
            return true;
        }

        // Form the next table according to the pivot value
        doPivotting(pivotRow, pivotColumn);

        return false;
    }

    bool checkOptimality()
    {
        // If there are no negative coefficients in C then it's optimal
        for (int i = 0; i < (int)C.size(); ++i)
        {
            if (C[i] < -1e-12f) // small tolerance
                return false;
        }
        // optimal: print & return true
        print();
        return true;
    }

    void doPivotting(int pivotRow, int pivotColumn)
    {
        // validate indices
        if (pivotRow < 0 || pivotRow >= rows || pivotColumn < 0 || pivotColumn >= cols)
        {
            cerr << "[doPivotting] invalid pivot indices: row=" << pivotRow << " col=" << pivotColumn << endl;
            return;
        }

        float pivetValue = A[pivotRow][pivotColumn]; // Gets the pivot value
        if (fabs(pivetValue) < 1e-12f)
        {
            cerr << "[doPivotting] pivot value is zero or nearly zero -> aborting pivot\n";
            return;
        }

        // Use vectors sized appropriately (avoid VLAs)
        vector<float> pivotRowVals(cols, 0.0f); // The pivot row
        vector<float> pivotColVals(rows, 0.0f); // The pivot column
        vector<float> rowNew(cols, 0.0f);       // The normalized pivot row

        // Update maximum -- the formula follows your earlier approach
        maximum = maximum - (C[pivotColumn] * (B[pivotRow] / pivetValue));

        // Get the row and column that has the pivot value
        for (int i = 0; i < cols; i++)
            pivotRowVals[i] = A[pivotRow][i];
        for (int j = 0; j < rows; j++)
            pivotColVals[j] = A[j][pivotColumn];

        // Normalize pivot row (divide by pivot)
        for (int k = 0; k < cols; k++)
            rowNew[k] = pivotRowVals[k] / pivetValue;

        // Normalize B[pivotRow]
        B[pivotRow] = B[pivotRow] / pivetValue;

        // Update other rows: A[m][p] = A[m][p] - pivotColVals[m] * rowNew[p]
        for (int m = 0; m < rows; m++)
        {
            if (m == pivotRow) continue;
            float multiplyValue = pivotColVals[m];
            for (int p = 0; p < cols; p++)
            {
                A[m][p] = A[m][p] - (multiplyValue * rowNew[p]);
            }
        }

        // Update B entries (except pivot row)
        for (int i = 0; i < (int)B.size(); i++)
        {
            if (i == pivotRow) continue;
            float multiplyValue = pivotColVals[i];
            B[i] = B[i] - (multiplyValue * B[pivotRow]);
        }

        // Update C: C = C - C[pivotColumn] * rowNew
        float multiplier = C[pivotColumn];
        for (int i = 0; i < cols; i++)
        {
            C[i] = C[i] - (multiplier * rowNew[i]);
        }

        // After row operations, set pivot column to 0 except pivot row
        for (int j = 0; j < rows; j++)
        {
            if (j == pivotRow) continue;
            A[j][pivotColumn] = 0.0f;
        }

        // Set the normalized pivot row back into A
        for (int i = 0; i < cols; i++)
            A[pivotRow][i] = rowNew[i];
    }

    // Print the current A array and B
    void print()
    {
        cout << "Current tableau (A | B):" << endl;
        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < cols; j++)
            {
                cout << A[i][j] << " ";
            }
            cout << "| " << B[i] << endl;
        }
        cout << "C: ";
        for (int j = 0; j < cols; ++j) cout << C[j] << " ";
        cout << endl;
    }

    // Find the least coefficients of constraints in the objective function's position
    // Return -1 if none negative (optimal)
    int findPivotColumn()
    {
        int location = -1;
        float minm = 0.0f;
        for (int i = 0; i < (int)C.size(); i++)
        {
            if (C[i] < minm - 1e-12f)
            {
                minm = C[i];
                location = i;
            }
        }
        return location;
    }

    // Find the row with the pivot value. Returns -1 for unbounded
    int findPivotRow(int pivotColumn)
    {
        if (pivotColumn < 0 || pivotColumn >= cols) return -1;

        vector<float> positiveValues(rows, 0.0f);
        vector<float> ratios(rows, std::numeric_limits<float>::infinity());
        int positiveCount = 0;

        for (int i = 0; i < rows; i++)
        {
            float val = A[i][pivotColumn];
            if (val > 1e-12f)
            {
                positiveValues[i] = val;
                // compute ratio
                if (fabs(val) > 1e-12f)
                    ratios[i] = B[i] / val;
                positiveCount++;
            }
            else
            {
                positiveValues[i] = 0.0f;
                ratios[i] = std::numeric_limits<float>::infinity();
            }
        }

        // If there are no positive entries in the pivot column -> unbounded
        if (positiveCount == 0)
        {
            isUnbounded = true;
            return -1;
        }

        // find minimal ratio
        float minimum = std::numeric_limits<float>::infinity();
        int location = -1;
        for (int i = 0; i < rows; ++i)
        {
            if (ratios[i] < minimum)
            {
                minimum = ratios[i];
                location = i;
            }
        }
        return location;
    }

    void CalculateSimplex()
    {
        bool end = false;

        cout << "initial array(Not optimal)" << endl;
        print();

        while (!end)
        {
            bool result = simplexAlgorithmCalculataion();
            if (result)
                end = true;
        }

        // compute solution vector (basic variables)
        solution.assign(cols, 0.0f);
        for (int col = 0; col < cols; ++col)
        {
            int oneIndex = -1;
            int zeroCount = 0;
            for (int r = 0; r < rows; ++r)
            {
                float v = A[r][col];
                if (fabs(v - 1.0f) < 1e-9f)
                {
                    if (oneIndex == -1) oneIndex = r;
                    else { oneIndex = -1; break; } // more than one '1', not a basic column
                }
                else if (fabs(v) < 1e-9f)
                {
                    zeroCount++;
                }
            }
            if (oneIndex != -1 && zeroCount == rows - 1)
            {
                solution[col] = B[oneIndex];
            }
            else
            {
                solution[col] = 0.0f;
            }
        }

        cout << "Answers for the Constraints of variables" << endl;
        for (int i = 0; i < cols; ++i)
        {
            cout << "variable " << (i + 1) << ": " << solution[i] << endl;
        }

        cout << endl;
        cout << "maximum value: " << maximum << endl; // print the maximum values
    }

    vector<float> getSolution() const { return solution; }
    float getMaximum() const { return maximum; }
};

// small helper: lowercase
static inline string toLower(const string &s)
{
    string r = s;
    transform(r.begin(), r.end(), r.begin(), ::tolower);
    return r;
}


int main(int argc, char** argv) {
    std::string filename = "input.txt";
    if (argc >= 2) filename = argv[1];

    std::ifstream in(filename);
    if (!in.is_open()) {
        std::cerr << "Failed to open input file: " << filename << std::endl;
        return 1;
    }

    auto toLower = [](std::string s){
        std::transform(s.begin(), s.end(), s.begin(), ::tolower);
        return s;
    };

    std::string token;
    if (!(in >> token)) {
        std::cerr << "Empty input file." << std::endl;
        return 1;
    }
    std::string opt = toLower(token);
    if (opt != "max" && opt != "min") {
        std::cerr << "First token must be 'max' or 'min' (found '" << token << "')." << std::endl;
        return 1;
    }

    int m, n;
    if (!(in >> m >> n)) {
        std::cerr << "Expecting two integers: m n (constraints, variables)." << std::endl;
        return 1;
    }
    if (m <= 0 || n <= 0) {
        std::cerr << "m and n must be positive." << std::endl;
        return 1;
    }

    std::vector<double> obj(n);
    for (int i = 0; i < n; ++i) {
        if (!(in >> obj[i])) {
            std::cerr << "Expecting " << n << " objective coefficients." << std::endl;
            return 1;
        }
    }

    std::string line;
    std::getline(in, line); // consume rest of line after objective coefficients

    int totalCols = n + m; // originals + one slack per constraint
    std::vector<std::vector<float>> A(m, std::vector<float>(totalCols, 0.0f));
    std::vector<float> B(m, 0.0f);

    for (int i = 0; i < m; ++i) {
        if (!std::getline(in, line)) {
            std::cerr << "Expecting " << m << " constraint lines but file ended early." << std::endl;
            return 1;
        }
        if (line.find_first_not_of(" \t\r\n") == std::string::npos) { --i; continue; } // skip blank lines

        std::stringstream ss(line);
        std::vector<double> coeffs(n);
        for (int j = 0; j < n; ++j) {
            if (!(ss >> coeffs[j])) {
                std::cerr << "Constraint " << (i+1) << ": expecting " << n << " coefficients." << std::endl;
                return 1;
            }
        }
        std::string rel;
        if (!(ss >> rel)) {
            std::cerr << "Constraint " << (i+1) << ": missing relation (<= or >=)." << std::endl;
            return 1;
        }
        double rhs;
        if (!(ss >> rhs)) {
            std::cerr << "Constraint " << (i+1) << ": missing RHS value." << std::endl;
            return 1;
        }

        if (rel == ">=") {
            for (int j = 0; j < n; ++j) coeffs[j] = -coeffs[j];
            rhs = -rhs;
        } else if (rel != "<=") {
            std::cerr << "Constraint " << (i+1) << ": relation must be '<=' or '>=' (found '" << rel << "')." << std::endl;
            return 1;
        }

        for (int j = 0; j < n; ++j) A[i][j] = static_cast<float>(coeffs[j]);
        A[i][n + i] = 1.0f; // slack variable
        B[i] = static_cast<float>(rhs);
    }

    // Build solver C vector using the internal sign convention:
    // for 'max' we use C[0..n-1] = -obj (slack coeffs are 0)
    std::vector<float> C(totalCols, 0.0f);
    if (opt == "max") {
        for (int i = 0; i < n; ++i) C[i] = static_cast<float>(-obj[i]);
    } else { // min -> convert to equivalent max
        for (int i = 0; i < n; ++i) C[i] = static_cast<float>(obj[i]);
    }

    std::cout << "Parsed problem: constraints=" << m << " originals=" << n
              << " totalColumns=" << totalCols << std::endl;

    // Construct and run Simplex (assumes Simplex available and API unchanged)
    Simplex solver(A, B, C);
    solver.CalculateSimplex();

    // If available, print getters (many implementations expose these)
    // (If your Simplex does not have getMaximum()/getSolution(), remove the lines below)
    try {
        std::vector<float> sol = solver.getSolution();
        float maxv = solver.getMaximum();
        std::cout << "\nResult (from getters):\n";
        for (size_t i = 0; i < sol.size(); ++i) {
            std::cout << "x" << (i+1) << " = " << sol[i] << "\n";
        }
        std::cout << "Objective = " << maxv << std::endl;
    } catch (...) {
        // if getters not present, ignore
    }

    return 0;
}
