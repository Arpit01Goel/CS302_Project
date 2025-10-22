#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <string>
#include <limits> // Required for numeric_limits

// To use this code, you must have the Eigen library installed and accessible.
// Compile with the appropriate flags for Eigen, e.g., g++ -I/path/to/eigen main.cpp
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

// --- Global Constants for Numerical Stability ---
// Tolerance to prevent division by zero in D_v^-2 calculation (Step 4)
const double EPSILON_FEASIBILITY = 1e-6;
// Tolerance for checking if the search direction is effectively zero (near optimum)
const double OPTIMALITY_TOL = 1e-6;


// Function to read the input file and populate the necessary parameters
bool readInput(const string& filename, MatrixXd& A, VectorXd& b, VectorXd& c, VectorXd& x0, double& gamma, int& m, int& n, int& max_iter) {
    ifstream infile(filename);
    if (!infile.is_open()) {
        cerr << "Error: Could not open input file " << filename << endl;
        return false;
    }

    string problem_type;
    infile >> problem_type; // Read and ignore the problem type (assuming maximization)

    if (!(infile >> m >> n)) {
        cerr << "Error reading m and n." << endl;
        return false;
    }

    if (!(infile >> gamma)) {
        cerr << "Error reading gamma." << endl;
        return false;
    }

    if (!(infile >> max_iter)) {
        cerr << "Error reading max_iter." << endl;
        return false;
    }

    // Resize matrices and vectors
    A.resize(m, n);
    b.resize(m);
    c.resize(n);
    x0.resize(n);

    // Read Matrix A
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            if (!(infile >> A(i, j))) {
                cerr << "Error reading matrix A at (" << i << ", " << j << ")." << endl;
                return false;
            }
        }
    }

    // Read Vector b
    for (int i = 0; i < m; ++i) {
        if (!(infile >> b(i))) {
            cerr << "Error reading vector b at index " << i << "." << endl;
            return false;
        }
    }

    // Read Vector c
    for (int i = 0; i < n; ++i) {
        if (!(infile >> c(i))) {
            cerr << "Error reading vector c at index " << i << "." << endl;
            return false;
        }
    }

    // Read Vector x0
    for (int i = 0; i < n; ++i) {
        if (!(infile >> x0(i))) {
            cerr << "Error reading vector x0 at index " << i << "." << endl;
            return false;
        }
    }

    cout << "Input read successfully. m=" << m << ", n=" << n << ", gamma=" << gamma << ", max_iter=" << max_iter << endl;
    return true;
}


// Simplified Affine-Scaling Karmarkar Algorithm (as per the image)
void simplifiedKarmarkar(const MatrixXd& A, const VectorXd& b, const VectorXd& c, const VectorXd& x0, double gamma, int max_iter) {
    int m = A.rows();
    int n = A.cols();

    // 1: k <- 0
    int k = 0;
    VectorXd x_k = x0; // x^(k)
    
    cout << "Starting Karmarkar's Algorithm..." << endl;

    // 2: while stopping criterion not met do
    while (k < max_iter) {
        cout << "\n--- Iteration " << k << " ---" << endl;

        // 3: v^(k) <- b - Ax^(k)
        VectorXd v_k = b - A * x_k;

        // --- NUMERICAL FIX 1: Prevent division by zero ---
        // Ensure v_k (slacks) are strictly positive (interior feasibility)
        if ((v_k.array() < EPSILON_FEASIBILITY).any()) {
            // Cap the slacks to ensure they are at least EPSILON_FEASIBILITY
            // This is crucial to prevent the 'nan' error.
            VectorXd v_k_capped = v_k.array().max(EPSILON_FEASIBILITY);
            v_k = v_k_capped;
            cout << "  * Adjusted slacks (v_k) to maintain interior feasibility." << endl;
        }

        // 4: D_v <- diag(v_1^(k), ..., v_m^(k)). Calculate D_v^-2 = diag(1/v_i^2)
        VectorXd v_k_sq_inv = v_k.array().pow(-2);
        MatrixXd D_v_sq_inv = v_k_sq_inv.asDiagonal();
        
        // 5: h_x <- (A^T D_v^-2 A)^-1 c
        MatrixXd M = A.transpose() * D_v_sq_inv * A;
        
        // Use LLT decomposition for solving the linear system M * h_x = c 
        // This is efficient and stable for the symmetric positive definite M.
        VectorXd h_x = M.llt().solve(c);

        // 6: h_v <- -A * h_x
        VectorXd h_v = -A * h_x;

        // --- NUMERICAL FIX 2: Check for Optimality/Unboundedness Robustly ---
        // Check objective change: c^T h_x should be close to zero at optimum
        double c_dot_hx = c.dot(h_x);

        // 7 & 8: Unbounded check (if h_v_i >= 0)
        // Check if all h_v components are non-negative (within tolerance)
        if ((h_v.array() >= -OPTIMALITY_TOL).all()) {
            if (c_dot_hx > OPTIMALITY_TOL) {
                 // c^T h_x > 0 means we can still move to infinity in a non-constraining direction
                 cout << "\nAlgorithm returned: Unbounded (h_v >= 0 and objective can increase)." << endl;
                 return;
            } else {
                 // h_v is all non-negative/zero AND objective increase is negligible (h_x is zero vector)
                 cout << "\nAlgorithm returned: Approximate Optimum reached (h_v and c^T h_x are near zero)." << endl;
                 break; // Exit loop for convergence
            }
        }
        
        // 9: end if
        
        // 10: alpha <- gamma * min_{i | (h_v)_i < 0} {-v_i^(k) / (h_v)_i}
        double alpha_max = numeric_limits<double>::max();
        bool found_negative_hv = false;

        for (int i = 0; i < m; ++i) {
            if (h_v(i) < -OPTIMALITY_TOL) { // Only consider significantly negative components
                double ratio = -v_k(i) / h_v(i);
                if (ratio < alpha_max) {
                    alpha_max = ratio;
                    found_negative_hv = true;
                }
            }
        }
        
        // The previous warning was here. With the EPSILON fix, this should rarely be needed,
        // but it provides a final safety net against numerical noise near optimum.
        if (!found_negative_hv) {
            // If the robust check in Step 7 didn't stop us, h_v must be numerically close to zero.
            // This is another indication of optimality. Stop.
            cout << "Warning: No sufficiently negative h_v found for ratio. Assuming convergence." << endl;
            break;
        }


        double alpha = gamma * alpha_max;
        cout << "  Alpha_max: " << alpha_max << ", Alpha: " << alpha << endl;


        // 11: x^(k+1) <- x^(k) + alpha * h_x
        x_k = x_k + alpha * h_x;
        
        // 12: k <- k + 1
        k++;

        // Calculate and print the current objective value
        double obj_value = c.dot(x_k);
        cout << "  Objective Value: " << obj_value << endl;
        
        // Optional: Check for negligible change in objective value for early exit
        // (Not strictly required by the image but good practice)
    }

    // 14: return x^(k) as approximate optimum
    cout << "\n\n--- Algorithm Finished ---" << endl;
    cout << "Stopping criterion (" << (k < max_iter ? "Convergence" : "Max Iterations") << ") reached." << endl;
    cout << "Approximate Optimum Solution x*: " << x_k.transpose() << endl;
    cout << "Maximum Objective Value c^T x*: " << c.dot(x_k) << endl;
}

// Main function
int main() {
    // Define the variables
    MatrixXd A;
    VectorXd b;
    VectorXd c;
    VectorXd x0;
    double gamma;
    int m, n;
    int max_iter;

    // IMPORTANT: The name of your input file
    string input_file = "input.txt";

    // Read the input
    if (!readInput(input_file, A, b, c, x0, gamma, m, n, max_iter)) {
        return 1; // Exit on error
    }

    // Run the Karmarkar Algorithm
    simplifiedKarmarkar(A, b, c, x0, gamma, max_iter);

    return 0;
}