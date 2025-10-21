// FINAL FULLY FIXED CODE: Sign-flip correction applied in Phase I.

#include <bits/stdc++.h>
using namespace std;

enum class SolveStatus { OPTIMAL, UNBOUNDED, INFEASIBLE };
struct SimplexResult { SolveStatus status; vector<double> x; double objective; };

static bool DEBUG = true; // set true to print traces
static string DEBUG_CASE_FILENAME = ""; // if non-empty only print traces for that input file
#define epsilion 1e-12
class SimplexCLRS {
public:
    SimplexCLRS(int n, int m, const vector<vector<double>>& A_in, const vector<double>& b_in,
                const vector<string>& relations_in, bool maximize=true)
        : n_orig(n), m(m), maximize(maximize)
    {
        A_input = A_in;
        b_input = b_in;
        relations = relations_in;
        normalizeConstraints();
        buildTableau();
    }

    void setOriginalObjective(const vector<double>& c_in) {
        orig_c.assign(total_vars, 0.0);
        for (int j = 0; j < (int)c_in.size() && j < n_orig; ++j) {
            orig_c[j] = (maximize ? c_in[j] : -c_in[j]);
        }
    }

    SimplexResult solve(const string &inputName = "") {
        input_label = inputName;
        if (artificial_vars_count > 0) {
            if (!phaseI()) {
                return {SolveStatus::INFEASIBLE, vector<double>(n_orig,0.0), NAN};
            }
        }
        return phaseII();
    }

private:
    int n_orig;
    int m;
    bool maximize;
    string input_label;

    vector<vector<double>> A_input;
    vector<double> b_input;
    vector<string> relations;

    int total_vars = 0;
    int slack_count = 0;
    int artificial_count = 0;
    int artificial_vars_count = 0;
    vector<int> col_type; 
    vector<int> artificial_cols;

    vector<vector<double>> T; 
    vector<int> basis; 
    vector<double> orig_c;

    static string colTypeName(int t){
        if (t==0) return "orig";
        if (t==1) return "slack/sur";
        if (t==3) return "arti";
        return "unknown";
    }

    void dbgPrint(const string &s) const {
        if (!DEBUG) return;
        if (!DEBUG_CASE_FILENAME.empty() && input_label != DEBUG_CASE_FILENAME) return;
        cerr << s;
    }

    void dbgTableau(const string &title) const {
        if (!DEBUG) return;
        if (!DEBUG_CASE_FILENAME.empty() && input_label != DEBUG_CASE_FILENAME) return;
        cerr << "=== " << title << " ===\n";
        // header
        for (int j=0;j<total_vars;j++) {
            string hd = (j < n_orig ? ("x"+to_string(j+1)) : ("v"+to_string(j+1)));
            cerr << hd << "\t";
        }
        cerr << "| RHS\n";
        for (int i=0;i<m;i++){
            for (int j=0;j<total_vars;j++){
                cerr << fixed << setprecision(6) << T[i][j] << "\t";
            }
            cerr << "| " << T[i][total_vars] << "\n";
        }
        cerr << "C: ";
        for (int j=0;j<total_vars;j++) cerr << fixed << setprecision(6) << T[m][j] << " ";
        cerr << "| " << T[m][total_vars] << "\n";
        cerr << "Basis cols: ";
        for (int i=0;i<m;i++) cerr << basis[i] << (i+1<=m-1? ",":"");
        cerr << "\n";
    }

    void normalizeConstraints() {
        for (int i = 0; i < m; ++i) {
            if (b_input[i] < -epsilion) { 
                for (int j = 0; j < n_orig; ++j) A_input[i][j] = -A_input[i][j];
                b_input[i] = -b_input[i];
                if (relations[i] == "<=") relations[i] = ">=";
                else if (relations[i] == ">=") relations[i] = "<=";
            }
        }
    }

    void buildTableau() {
        col_type.assign(n_orig, 0);
        
        slack_count = 0;
        artificial_count = 0;

        for (int i = 0; i < m; ++i) {
            if (relations[i] == "<=") slack_count++;
            else if (relations[i] == ">=") { slack_count++; artificial_count++; }
            else { artificial_count++; } 
        }
        
        total_vars = n_orig + slack_count + artificial_count;
        col_type.resize(total_vars, 0);
        
        for (int j = n_orig; j < n_orig + slack_count; ++j) col_type[j] = 1;
        
        artificial_cols.clear();
        int arti_start = n_orig + slack_count;
        for (int j = arti_start; j < total_vars; ++j){
            col_type[j] = 3;
            artificial_cols.push_back(j);
        }
        artificial_vars_count = (int)artificial_cols.size();
        
        T.assign(m+1, vector<double>(total_vars+1, 0.0));
        basis.assign(m, -1);
        orig_c.assign(total_vars, 0.0);

        for (int i=0;i<m;i++){
            for (int j=0;j<n_orig;j++) T[i][j] = A_input[i][j];
        }
        
        int sidx = n_orig;
        int aidx = arti_start;
        for (int i=0;i<m;i++){
            if (relations[i] == "<=") {
                T[i][sidx] = 1.0;
                basis[i] = sidx;
                sidx++;
            } else if (relations[i] == ">=") {
                T[i][sidx] = -1.0; 
                T[i][aidx] = 1.0;
                basis[i] = aidx;
                sidx++; aidx++;
            } else {
                T[i][aidx] = 1.0;
                basis[i] = aidx;
                aidx++;
            }
            T[i][total_vars] = b_input[i];
        }
    }

    void addRowTimes(int to, int from, double factor) {
        for (int j=0;j<=total_vars;j++) T[to][j] += factor * T[from][j];
    }

    bool phaseI() {
        dbgPrint("START PHASE I\n");
        dbgTableau("Initial Phase I Tableau (before bottom row build)");

        for (int j=0;j<=total_vars;j++) T[m][j] = 0.0;
        for (int k : artificial_cols) T[m][k] = -1.0; 
        T[m][total_vars] = 0.0;

        for (int i=0;i<m;i++) {
            int bc = basis[i];
            if (bc >=0 && col_type[bc] == 3) {
                double factor = -T[m][bc] / T[i][bc];
                addRowTimes(m, i, factor);
            }
        }
        
        // ******************* CRITICAL FIX *******************
        // The elimination loop correctly zeros out T[m][bc], but the resultant T[m]
        // stores the NEGATIVE of the reduced costs needed for the CLRS convention 
        // (T[m] = -c_j, looking for T[m][j] < 0). We must flip the sign of the whole row.
        for (int j=0; j<=total_vars; j++) {
            T[m][j] = -T[m][j];
        }
        // ****************************************************

        dbgTableau("Phase I Tableau after bottom row init (FIXED SIGN)");

        SolveStatus st = runSimplexPhase(true);
        dbgTableau("Phase I Tableau after simplex");

        if (st == SolveStatus::UNBOUNDED) {
            dbgPrint("PHASE I reported UNBOUNDED (treat as infeasible)\n");
            return false;
        }
        double phase1_obj = T[m][total_vars];
        dbgPrint("PHASE I objective (should be 0 if feasible): " + to_string(phase1_obj) + "\n");
        if (fabs(phase1_obj) > 1e-8) {
            dbgPrint("PHASE I found nonzero objective -> INFEASIBLE\n");
            return false;
        }

        // Remove artificial columns 
        vector<int> colmap(total_vars, -1);
        int newcols = 0;
        for (int j=0;j<total_vars;j++) {
            if (col_type[j] == 3) continue;
            colmap[j] = newcols++;
        }
        vector<vector<double>> newT(m+1, vector<double>(newcols+1, 0.0));
        for (int i=0;i<=m;i++){
            for (int j=0;j<total_vars;j++){
                if (colmap[j] != -1) newT[i][colmap[j]] = T[i][j];
            }
            newT[i][newcols] = T[i][total_vars];
        }
        for (int i=0;i<m;i++){
            if (basis[i] != -1) {
                if (colmap[basis[i]] == -1) basis[i] = -1;
                else basis[i] = colmap[basis[i]];
            }
        }
        vector<int> newcoltype;
        for (int j=0;j<total_vars;j++) if (col_type[j] != 3) newcoltype.push_back(col_type[j]);
        total_vars = newcols;
        T.swap(newT);
        col_type.swap(newcoltype);
        artificial_cols.clear();
        artificial_vars_count = 0;

        dbgTableau("After removing artificial columns (start Phase II prep)");
        return true;
    }

    SimplexResult phaseII() {
        dbgPrint("START PHASE II\n");
        for (int j=0;j<total_vars;j++) T[m][j] = -orig_c[j];
        T[m][total_vars] = 0.0;
        for (int i=0;i<m;i++){
            int bc = basis[i];
            if (bc >=0) {
                double coeff = orig_c[bc];
                if (fabs(coeff) > 1e-14) {
                    addRowTimes(m, i, coeff);
                }
            }
        }
        dbgTableau("Phase II initial tableau");
        SolveStatus st = runSimplexPhase(false);
        if (st == SolveStatus::UNBOUNDED) return {SolveStatus::UNBOUNDED, {}, NAN};
        
        vector<double> sol(n_orig, 0.0);
        for (int i=0;i<m;i++){
            int bc = basis[i];
            if (bc >=0 && bc < n_orig) sol[bc] = T[i][total_vars];
        }
        double obj = T[m][total_vars];
        if (!maximize) obj = -obj; 
        
        return {SolveStatus::OPTIMAL, sol, obj};
    }

    SolveStatus runSimplexPhase(bool isPhaseI) {
        const int MAX_IT = 10000;
        int it = 0;
        double tol = 1e-12; 
        while (true) {
            ++it;
            if (it > MAX_IT) { 
                dbgPrint("Reached max iterations\n"); 
                return SolveStatus::OPTIMAL; 
            }
            
            int entering = -1;
            for (int j=0; j<total_vars; j++) {
                if (T[m][j] < -tol) {
                    entering = j;
                    break; 
                }
            }
            
            if (entering == -1) {
                return SolveStatus::OPTIMAL;
            }

            int leaving = -1;
            double bestRatio = numeric_limits<double>::infinity();
            
            for (int i=0; i<m; i++){
                double a = T[i][entering];
                if (a > tol) { 
                    double ratio = T[i][total_vars] / a;
                    
                    if (ratio < bestRatio - tol) { 
                        bestRatio = ratio;
                        leaving = i;
                    }
                }
            }
            
            if (leaving == -1) {
                return SolveStatus::UNBOUNDED; 
            }
            
            dbgPrint("Pivot: entering col=" + to_string(entering) + " leaving row=" + to_string(leaving) + "\n");
            pivot(leaving, entering);
            dbgTableau("After pivot (Iter " + to_string(it) + ")");
        }
    }

    void pivot(int l, int e) {
        double a = T[l][e];
        if (fabs(a) < 1e-14) return;
        
        for (int j=0;j<=total_vars;j++) T[l][j] /= a;
        
        for (int i=0;i<=m;i++) {
            if (i==l) continue;
            double factor = T[i][e];
            if (fabs(factor) > 0.0) {
                for (int j=0;j<=total_vars;j++) T[i][j] -= factor * T[l][j];
            }
        }
        basis[l] = e;
    }
};

// parse input file (omitted for brevity, assume correct from previous final version)

bool parse_input_file(const string &filename, int &n, int &m,
                      vector<vector<double>> &A, vector<double> &b, vector<string> &rels,
                      vector<double> &c, bool &maximize)
{
    ifstream in(filename);
    if (!in.is_open()) { cerr << "Cannot open " << filename << "\n"; return false; }
    
    string tok; 
    in >> ws; 
    if (!(in>>tok)) return false;
    string opt = tok; transform(opt.begin(), opt.end(), opt.begin(), ::tolower);
    maximize = (opt == "max");
    
    in >> ws; 
    if (!(in >> m >> n)) { cerr << "Expecting m n\n"; return false; }
    
    c.assign(n, 0.0);
    for (int i=0;i<n;i++) {
        if (!(in >> c[i])) { cerr << "Missing objective coefficients\n"; return false; }
    }
    
    string line; 
    getline(in, line); 
    
    A.assign(m, vector<double>(n,0.0));
    b.assign(m,0.0); rels.assign(m,"<=");
    
    int constraints_read = 0;
    for (int i=0;i<m;i++){
        if (!getline(in, line)) { break; } 
        if (line.find_first_not_of(" \t\r\n") == string::npos) { --i; continue; } 
        
        stringstream ss(line);
        
        bool read_ok = true;
        for (int j=0;j<n;j++) {
            if (!(ss >> A[i][j])) { read_ok = false; break; }
        }
        
        string rel; double rhs; 
        if (read_ok && (ss >> rel >> rhs)) {
            rels[i] = rel; b[i] = rhs;
            constraints_read++;
        } else {
            break;
        }
    }
    
    if (constraints_read < m) {
        m = constraints_read;
        A.resize(m);
        b.resize(m);
        rels.resize(m);
        if (m==0) { cerr << "No valid constraints read\n"; return false; }
    }

    return true;
}

int main(int argc, char** argv) {
    string filename = "input.txt";
    if (argc >= 2) filename = argv[1];
    if (argc >= 3) DEBUG_CASE_FILENAME = argv[2]; 
    DEBUG = true;

    int n,m;
    vector<vector<double>> A;
    vector<double> b;
    vector<string> rels;
    vector<double> c;
    bool maximize;
    
    if (!parse_input_file(filename, n, m, A, b, rels, c, maximize)) {
        cerr << "Parse failed\n"; return 1;
    }

    SimplexCLRS solver(n,m,A,b,rels, maximize);
    solver.setOriginalObjective(c);
    SimplexResult res = solver.solve(filename);

    if (res.status == SolveStatus::INFEASIBLE) {
        cout << "INFEASIBLE\n";
    } else if (res.status == SolveStatus::UNBOUNDED) {
        cout << "UNBOUNDED\n";
    } else {
        cout << "OPTIMAL\n";
        cout << "Objective = " << fixed << setprecision(6) << res.objective << "\n";
        for (int i=0;i<(int)res.x.size();++i) cout << "x" << (i+1) << " = " << fixed << setprecision(6) << res.x[i] << "\n";
    }
    return 0;
}