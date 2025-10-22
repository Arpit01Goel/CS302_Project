// karmarkar_affine_fixed.cpp
// Simplified Karmarkar (affine-scaling) implementation (dense, educational).
// Fixes: correct Cholesky + robust regularization + better I/O errors.

#include <cmath>
#include <iostream>
#include <vector>
#include <fstream>
#include <random>
#include <limits>
#include <chrono>
#include <string>

using dvec = std::vector<double>;
using dmat = std::vector<double>; // row-major m x n stored as data[i*n + j]

static inline double EPS_DBL = 1e-15;

// --- Basic linear algebra helpers ---

// y = A * x  (A: m x n row-major)
void mat_vec_mul(const dmat &A, int m, int n, const dvec &x, dvec &y) {
    y.assign(m, 0.0);
    for (int i = 0; i < m; ++i) {
        double tmp = 0.0;
        const double *row = &A[i * n];
        for (int j = 0; j < n; ++j) tmp += row[j] * x[j];
        y[i] = tmp;
    }
}

// y = A^T * x  (A: m x n)
void matT_vec_mul(const dmat &A, int m, int n, const dvec &x, dvec &y) {
    y.assign(n, 0.0);
    for (int i = 0; i < m; ++i) {
        const double *row = &A[i * n];
        double xi = x[i];
        for (int j = 0; j < n; ++j) y[j] += row[j] * xi;
    }
}

// Build M = A^T * W * A  where W is diag(w[k]) (m x m). Result M is n x n symmetric row-major.
void build_AtWA(const dmat &A, int m, int n, const dvec &w, dmat &M) {
    M.assign(n * n, 0.0);
    // M_{ij} = sum_{k=0..m-1} A[k][i] * w[k] * A[k][j]
    for (int k = 0; k < m; ++k) {
        double wk = w[k];
        const double *row = &A[k * n];
        for (int i = 0; i < n; ++i) {
            double aik = row[i];
            double s = wk * aik;
            double *Mrow_i = &M[i * n];
            // accumulate into lower triangle and diagonal
            for (int j = 0; j <= i; ++j) {
                Mrow_i[j] += s * row[j];
            }
        }
    }
    // copy lower triangle to upper triangle
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < i; ++j) {
            M[j * n + i] = M[i * n + j];
        }
    }
}

// Standard Cholesky (in-place). On return, L stored in lower triangle of M.
// Returns true if successful (positive definite).
bool cholesky_inplace(dmat &M, int n, double tiny_reg = 0.0) {
    // optionally add tiny_reg to diagonal before factoring
    if (tiny_reg != 0.0) {
        for (int i = 0; i < n; ++i) M[i * n + i] += tiny_reg;
    }
    for (int k = 0; k < n; ++k) {
        double sum = M[k * n + k];
        for (int t = 0; t < k; ++t) {
            double Lkt = M[k * n + t];
            sum -= Lkt * Lkt;
        }
        if (sum <= 0.0) return false;
        double Lkk = std::sqrt(sum);
        M[k * n + k] = Lkk;
        for (int i = k + 1; i < n; ++i) {
            double s = M[i * n + k];
            for (int t = 0; t < k; ++t) s -= M[i * n + t] * M[k * n + t];
            M[i * n + k] = s / Lkk;
        }
    }
    // zero upper triangle (optional)
    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j)
            M[i * n + j] = 0.0;
    return true;
}

// Solve L * y = b where L is lower-triangular stored in "chol" (n x n)
void solve_lower(const dmat &L, int n, const dvec &b, dvec &y) {
    y.assign(n, 0.0);
    for (int i = 0; i < n; ++i) {
        double s = b[i];
        for (int j = 0; j < i; ++j) s -= L[i * n + j] * y[j];
        y[i] = s / L[i * n + i];
    }
}

// Solve L^T * x = y where L stored in chol (lower triangle)
void solve_upper_trans(const dmat &L, int n, const dvec &y, dvec &x) {
    x.assign(n, 0.0);
    for (int ii = n - 1; ii >= 0; --ii) {
        int i = ii;
        double s = y[i];
        for (int j = i + 1; j < n; ++j) s -= L[j * n + i] * x[j]; // L^T_ij = L_ji
        x[i] = s / L[i * n + i];
    }
}

// Solve symmetric positive-definite system M x = b using Cholesky in-place.
// Will try diagonal regularization if needed. Returns true if solved.
bool solve_spd_with_regularization(dmat M, int n, const dvec &b, dvec &x) {
    const double base_reg = 1e-12;
    const int MAX_TRIES = 6;
    double reg = 0.0;
    for (int attempt = 0; attempt < MAX_TRIES; ++attempt) {
        bool ok = cholesky_inplace(M, n, reg);
        if (ok) {
            dvec y;
            solve_lower(M, n, b, y);
            solve_upper_trans(M, n, y, x);
            return true;
        }
        // increase reg and retry
        reg = (attempt == 0) ? base_reg : reg * 10.0;
        // reload original diagonal + reg requires M be rebuilt by caller; here we already have copy M param so loop will add reg in cholesky
    }
    return false;
}

// Euclidean norm
double norm2(const dvec &v) {
    double s = 0.0;
    for (double x : v) s += x * x;
    return std::sqrt(s);
}

// read matrix A (m x n) from file; expects row-major with whitespace separation
bool read_matrix_file(const std::string &path, dmat &A, int m, int n) {
    std::ifstream in(path);
    if (!in) {
        std::cerr << "Error: cannot open matrix file '" << path << "'\n";
        return false;
    }
    A.assign(m * n, 0.0);
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            if (!(in >> A[i * n + j])) {
                std::cerr << "Error: failed reading A[" << i << "][" << j << "] from " << path << "\n";
                return false;
            }
        }
    }
    return true;
}

bool read_vector_file(const std::string &path, dvec &v, int len) {
    std::ifstream in(path);
    if (!in) {
        std::cerr << "Error: cannot open vector file '" << path << "'\n";
        return false;
    }
    v.assign(len, 0.0);
    for (int i = 0; i < len; ++i) {
        if (!(in >> v[i])) {
            std::cerr << "Error: failed reading element " << i << " from " << path << "\n";
            return false;
        }
    }
    return true;
}

// Simple demo problem generator (small)
void generate_demo(dmat &A, dvec &b, dvec &c, int m, int n) {
    std::mt19937_64 rng(123456);
    std::uniform_real_distribution<double> U(0.1, 2.0);
    A.assign(m * n, 0.0);
    for (int i = 0; i < m * n; ++i) A[i] = U(rng);
    c.assign(n, 0.0); for (int i = 0; i < n; ++i) c[i] = U(rng);
    // build b = A * x0 with x0 positive to ensure feasibility
    dvec x0(n, 1.0 / n);
    mat_vec_mul(A, m, n, x0, b);
    for (int i = 0; i < m; ++i) b[i] += 0.5;
}

// --- Karmarkar simplified affine-scaling algorithm ---
bool karmarkar_affine(const dmat &A, const dvec &b, const dvec &c,
                      dvec &x, int maxIter, double tol, double gamma) {
    int m = (int)b.size();
    int n = (int)c.size();
    if (n == 0 || m == 0) return false;

    dvec v(m), Dinv2(m), h_x(n), h_v(m);

    for (int iter = 0; iter < maxIter; ++iter) {
        // v = b - A x
        mat_vec_mul(A, m, n, x, v);
        for (int i = 0; i < m; ++i) v[i] = b[i] - v[i];

        // check strictly interior: v_i > 0
        for (int i = 0; i < m; ++i) {
            if (!(v[i] > 0.0)) {
                std::cerr << "Infeasible or on boundary at iter " << iter << ": v[" << i << "] = " << v[i] << "\n";
                return false;
            }
        }

        // Dinv2 = 1 / v_i^2
        for (int i = 0; i < m; ++i) {
            double vv = v[i];
            if (vv <= 0.0) { Dinv2[i] = 1e300; } else { Dinv2[i] = 1.0 / (vv * vv); }
        }

        // Build M = A^T * W * A  where W = diag(Dinv2)
        dmat M;
        build_AtWA(A, m, n, Dinv2, M);

        // Solve M * h_x = c with robust regularization
        bool ok = solve_spd_with_regularization(M, n, c, h_x);
        if (!ok) {
            std::cerr << "Failed to solve M h_x = c even after regularization (iter " << iter << ")\n";
            return false;
        }

        // h_v = -A * h_x
        mat_vec_mul(A, m, n, h_x, h_v);
        for (int i = 0; i < m; ++i) h_v[i] = -h_v[i];

        // If any h_v >= 0 then unbounded (per alg)
        for (int i = 0; i < m; ++i) {
            if (h_v[i] >= 0.0) {
                std::cout << "Unbounded detected (h_v[" << i << "] >= 0) at iter " << iter << "\n";
                return false;
            }
        }

        // alpha = gamma * min{-v_i / h_v_i | h_v_i < 0}
        double alpha = std::numeric_limits<double>::infinity();
        for (int i = 0; i < m; ++i) {
            if (h_v[i] < 0.0) {
                double cand = -v[i] / h_v[i];
                if (cand < alpha) alpha = cand;
            }
        }
        if (!std::isfinite(alpha) || alpha <= 0.0) {
            std::cerr << "Bad alpha (nonfinite or <=0) at iter " << iter << "\n";
            return false;
        }
        alpha *= gamma;

        // Update x <- x + alpha * h_x
        for (int j = 0; j < n; ++j) x[j] += alpha * h_x[j];

        // stopping criteria: small step (||alpha*h_x||) or small ||h_x||
        double step_norm = alpha * norm2(h_x);
        double hx_norm = norm2(h_x);
        if (step_norm < tol || hx_norm < tol) {
            std::cout << "Converged (iter=" << iter << ") step_norm=" << step_norm << " hx_norm=" << hx_norm << "\n";
            return true;
        }

        if ((iter & 15) == 0) {
            double obj = 0.0;
            for (int j = 0; j < n; ++j) obj += c[j] * x[j];
            std::cout << "iter " << iter << " alpha=" << alpha << " obj=" << obj << " ||h_x||=" << hx_norm << "\n";
        }
    }
    std::cout << "Reached maxIter\n";
    return true;
}

int main(int argc, char **argv) {
    int maxIter = 1000;
    double tol = 1e-8;
    double gamma = 0.9;
    dmat A;
    dvec b, c;

    if (argc == 1) {
        // demo with explicit b
        int m = 3, n = 3;
        A = {1,1,1,
             2,1,0,
             0,1,2};
        b = {5,8,6};
        c = {-3,-2,0};
        std::cout << "Running demo (hardcoded A,b,c) m=" << m << " n=" << n << "\n";
    } else if (argc >= 6) {
        // Usage: prog Afile bfile cfile m n [maxIter tol gamma]
        std::string Afile = argv[1], bfile = argv[2], cfile = argv[3];
        int m = std::stoi(argv[4]);
        int n = std::stoi(argv[5]);
        if (!read_matrix_file(Afile, A, m, n)) { std::cerr << "Failed read A\n"; return 1; }
        if (!read_vector_file(bfile, b, m)) { std::cerr << "Failed read b\n"; return 1; }
        if (!read_vector_file(cfile, c, n)) { std::cerr << "Failed read c\n"; return 1; }
        if (argc >= 7) maxIter = std::stoi(argv[6]);
        if (argc >= 8) tol = std::stod(argv[7]);
        if (argc >= 9) gamma = std::stod(argv[8]);
        std::cout << "Loaded A(" << m << "x" << n << "), maxIter=" << maxIter << " tol=" << tol << " gamma=" << gamma << "\n";
    } else {
        std::cerr << "Usage: " << argv[0] << " [Afile bfile cfile m n [maxIter tol gamma]]\n";
        return 1;
    }

    int m = (int)b.size();
    int n = (int)c.size();
    if (n == 0 || m == 0) { std::cerr << "Empty problem\n"; return 1; }

    // Initial interior point x0 strictly positive and (ideally) satisfying A x < b.
    dvec x(n, 1.0 / n);
    // ensure Ax < b: if not, try scaling down x
    dvec Ax;
    mat_vec_mul(A, m, n, x, Ax);
    double max_ratio = 0.0;
    for (int i = 0; i < m; ++i) {
        if (Ax[i] >= b[i]) {
            double r = (Ax[i] + 1e-12) / std::max(1e-12, b[i]);
            max_ratio = std::max(max_ratio, r);
        }
    }
    if (max_ratio > 0.0) {
        for (double &xi : x) xi *= (0.5 / (1.0 + max_ratio)); // shrink to get interior
    }

    auto t0 = std::chrono::high_resolution_clock::now();
    bool ok = karmarkar_affine(A, b, c, x, maxIter, tol, gamma);
    auto t1 = std::chrono::high_resolution_clock::now();

    std::cout << "Done. ok=" << ok << " time=" << std::chrono::duration<double>(t1 - t0).count() << "s\n";
    std::cout << "Final x (first 20):\n";
    for (int i = 0; i < (int)x.size() && i < 20; ++i) std::cout << x[i] << " ";
    std::cout << "\nObjective c^T x = ";
    double obj = 0.0;
    for (int j = 0; j < n; ++j) obj += c[j] * x[j];
    std::cout << obj << "\n";
    return 0;
}

