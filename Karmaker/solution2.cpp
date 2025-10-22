#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <string>
#include <limits>
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

// --- Global Constants for Numerical Stability ---
const double EPSILON_FEASIBILITY = 1e-8;  // Increased from 1e-9
const double OPTIMALITY_TOL = 1e-6;
const double CONVERGENCE_TOL = 1e-8;

bool readInput(const string& filename, MatrixXd& A, VectorXd& b, VectorXd& c, VectorXd& x0, double& gamma, int& m, int& n, int& max_iter) {
    ifstream infile(filename);
    if (!infile.is_open()) {
        cerr << "Error: Could not open input file " << filename << endl;
        return false;
    }

    string problem_type;
    infile >> problem_type;

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

    A.resize(m, n);
    b.resize(m);
    c.resize(n);
    x0.resize(n);

    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            if (!(infile >> A(i, j))) {
                cerr << "Error reading matrix A at (" << i << ", " << j << ")." << endl;
                return false;
            }
        }
    }

    for (int i = 0; i < m; ++i) {
        if (!(infile >> b(i))) {
            cerr << "Error reading vector b at index " << i << "." << endl;
            return false;
        }
    }

    for (int i = 0; i < n; ++i) {
        if (!(infile >> c(i))) {
            cerr << "Error reading vector c at index " << i << "." << endl;
            return false;
        }
    }

    for (int i = 0; i < n; ++i) {
        if (!(infile >> x0(i))) {
            cerr << "Error reading vector x0 at index " << i << "." << endl;
            return false;
        }
    }

    cout << "Input read successfully. m=" << m << ", n=" << n << ", gamma=" << gamma << ", max_iter=" << max_iter << endl;
    return true;
}

void simplifiedKarmarkar(const MatrixXd& A, const VectorXd& b, const VectorXd& c, const VectorXd& x0, double gamma, int max_iter) {
    int m = A.rows();
    int n = A.cols();

    int k = 0;
    VectorXd x_k = x0;
    double prev_obj_value = c.dot(x_k);
    
    cout << "Starting Karmarkar's Algorithm..." << endl;
    cout << "Initial Objective Value: " << prev_obj_value << endl;

    while (k < max_iter) {
        cout << "\n--- Iteration " << k << " ---" << endl;

        // Step 3: v^(k) <- b - Ax^(k)
        VectorXd v_k = b - A * x_k;

        // Check if all slacks are positive (interior feasibility)
        double min_slack = v_k.minCoeff();
        
        // If we're too close to boundary, we need stronger adjustment
        if (min_slack <= EPSILON_FEASIBILITY) {
            cout << "Warning: Point is at/near boundary (min slack: " << min_slack << ")." << endl;
            
            // Project back to interior if needed
            if (min_slack <= 0) {
                cout << "  Projecting point back to interior..." << endl;
                v_k = v_k.array().max(1e-6); // Use larger epsilon for recovery
            } else {
                v_k = v_k.array().max(EPSILON_FEASIBILITY);
            }
        }
        
        // Additional check: if slacks are too small, matrix M becomes ill-conditioned
        double max_v_inv_sq = v_k.array().pow(-2).maxCoeff();
        if (max_v_inv_sq > 1e12) {
            cout << "Warning: Very small slack detected (conditioning issue). Adjusting..." << endl;
            // Increase minimum slack to prevent numerical issues
            v_k = v_k.array().max(1e-5);
        }

        // Step 4: D_v^-2 = diag(1/v_i^2)
        VectorXd v_k_sq_inv = v_k.array().pow(-2);
        MatrixXd D_v_sq_inv = v_k_sq_inv.asDiagonal();
        
        // Step 5: h_x <- (A^T D_v^-2 A)^-1 c
        MatrixXd M = A.transpose() * D_v_sq_inv * A;
        
        // Add regularization to ensure M is positive definite
        // This is crucial when constraints are nearly redundant or point is near boundary
        double regularization = 1e-8;
        MatrixXd M_reg = M + regularization * MatrixXd::Identity(n, n);
        
        // Try LLT decomposition first (for symmetric positive definite)
        LLT<MatrixXd> llt_of_M(M_reg);
        VectorXd h_x;
        
        if (llt_of_M.info() == Success) {
            h_x = llt_of_M.solve(c);
        } else {
            // If LLT fails, use more robust LDLT or full pivot LU
            cout << "  Warning: Using LDLT decomposition (M may be ill-conditioned)" << endl;
            LDLT<MatrixXd> ldlt_of_M(M_reg);
            if (ldlt_of_M.info() == Success) {
                h_x = ldlt_of_M.solve(c);
            } else {
                // Last resort: use FullPivLU (most robust but slower)
                cout << "  Warning: Using FullPivLU decomposition" << endl;
                h_x = M_reg.fullPivLu().solve(c);
            }
        }

        // Step 6: h_v <- -A * h_x
        VectorXd h_v = -A * h_x;

        // Check for optimality/convergence
        double c_dot_hx = c.dot(h_x);
        double h_x_norm = h_x.norm();
        
        cout << "  c^T h_x: " << c_dot_hx << ", ||h_x||: " << h_x_norm << endl;

        // If search direction is too small, we've converged
        if (h_x_norm < OPTIMALITY_TOL) {
            cout << "\nConverged: Search direction is negligible." << endl;
            break;
        }

        // Step 7-8: Check for unboundedness
        bool all_hv_nonnegative = true;
        for (int i = 0; i < m; ++i) {
            if (h_v(i) < -OPTIMALITY_TOL) {
                all_hv_nonnegative = false;
                break;
            }
        }
        
        if (all_hv_nonnegative) {
            if (c_dot_hx > OPTIMALITY_TOL) {
                cout << "\nAlgorithm returned: Unbounded (h_v >= 0 and objective can increase)." << endl;
                return;
            } else {
                cout << "\nConverged: Optimal point reached (h_v >= 0, negligible improvement)." << endl;
                break;
            }
        }

        // Step 10: alpha <- gamma * min_{i | (h_v)_i < 0} {-v_i^(k) / (h_v)_i}
        double alpha_max = numeric_limits<double>::max();
        bool found_negative_hv = false;

        for (int i = 0; i < m; ++i) {
            if (h_v(i) < -OPTIMALITY_TOL) {
                double ratio = -v_k(i) / h_v(i);
                if (ratio > 0 && ratio < alpha_max) {
                    alpha_max = ratio;
                    found_negative_hv = true;
                }
            }
        }
        
        if (!found_negative_hv) {
            cout << "Converged: No sufficiently negative h_v components." << endl;
            break;
        }

        double alpha = gamma * alpha_max;
        cout << "  Alpha_max: " << alpha_max << ", Alpha: " << alpha << endl;

        // Step 11: x^(k+1) <- x^(k) + alpha * h_x
        x_k = x_k + alpha * h_x;
        
        // Ensure x_k remains positive (shouldn't be needed but safety check)
        if ((x_k.array() < 0).any()) {
            cout << "Warning: Some variables became negative. Capping to zero." << endl;
            x_k = x_k.array().max(EPSILON_FEASIBILITY);
        }

        k++;

        // Calculate and print the current objective value
        double obj_value = c.dot(x_k);
        double obj_change = abs(obj_value - prev_obj_value);
        cout << "  Objective Value: " << obj_value << " (change: " << obj_change << ")" << endl;
        
        // Check for convergence based on objective value change
        if (obj_change < CONVERGENCE_TOL && k > 5) {
            cout << "\nConverged: Objective value change is negligible." << endl;
            break;
        }
        
        prev_obj_value = obj_value;
    }

    cout << "\n\n--- Algorithm Finished ---" << endl;
    cout << "Total Iterations: " << k << endl;
    cout << "Stopping Reason: " << (k < max_iter ? "Convergence" : "Max Iterations") << endl;
    
    // Final feasibility check
    VectorXd final_slack = b - A * x_k;
    double min_final_slack = final_slack.minCoeff();
    bool is_feasible = (x_k.array() >= -EPSILON_FEASIBILITY).all() && (final_slack.array() >= -EPSILON_FEASIBILITY).all();
    
    cout << "\nFinal Solution:" << endl;
    cout << "x* = " << x_k.transpose() << endl;
    cout << "Objective Value c^T x* = " << c.dot(x_k) << endl;
    cout << "Feasibility: " << (is_feasible ? "YES" : "NO") << endl;
    cout << "Minimum Slack: " << min_final_slack << endl;
}

int main() {
    MatrixXd A;
    VectorXd b;
    VectorXd c;
    VectorXd x0;
    double gamma;
    int m, n;
    int max_iter;

    string input_file = "input.txt";

    if (!readInput(input_file, A, b, c, x0, gamma, m, n, max_iter)) {
        return 1;
    }

    simplifiedKarmarkar(A, b, c, x0, gamma, max_iter);

    return 0;
}