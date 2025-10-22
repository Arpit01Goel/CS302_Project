import os

# Define the directory for test cases
TEST_DIR = "test_cases"

# Define 10 rigorous test cases:
# Format: [M, N, Gamma, Max_Iter, A (flat list), B (list), C (list), X0 (list), Expected_Obj, Tolerance]
# Note: Max_Iter is set higher for safety in convergence. Gamma is 0.95.
TEST_CASES = [
    # TC 1: Simple 2x2 (Corner Solution) -> Expected Obj: 9.0
    # Max 3x1 + 2x2; x1+x2 <= 4, 2x1+x2 <= 5
    [2, 2, 0.95, 200, [1.0, 1.0, 2.0, 1.0], [4.0, 5.0], [3.0, 2.0], [0.5, 0.5], 9.0, 1e-4],

    # TC 2: Different Optimal Corner -> Expected Obj: 12.0
    # Max 4x1 + x2; x1 <= 2, x2 <= 5, x1+x2 <= 6 (Opt: x1=2, x2=4)
    [3, 2, 0.95, 200, [1.0, 0.0, 0.0, 1.0, 1.0, 1.0], [2.0, 5.0, 6.0], [4.0, 1.0], [1.0, 1.0], 12.0, 1e-4],

    # TC 3: Single Variable (1D problem) -> Expected Obj: 10.0
    # Max 5x1; 2x1 <= 4, x1 <= 2 (Opt: x1=2)
    [2, 1, 0.95, 200, [2.0, 1.0], [4.0, 2.0], [5.0], [1.0], 10.0, 1e-4],
    
    # TC 4: Redundant Constraints/Slack (Interior solution) -> Expected Obj: 7.0
    # Max 7x1 + 0x2; x1+x2 <= 10, x1 <= 1 (Opt: x1=1, x2=0)
    [2, 2, 0.95, 200, [1.0, 1.0, 1.0, 0.0], [10.0, 1.0], [7.0, 0.0], [0.1, 0.1], 7.0, 1e-4],
    
    # TC 5: More Variables than Constraints (2x3) -> Expected Obj: 5.0
    # Max x1 + x2 + x3; x1+x2 <= 2, x2+x3 <= 3 (Opt: x1=2, x2=0, x3=3)
    [2, 3, 0.95, 300, [1.0, 1.0, 0.0, 0.0, 1.0, 1.0], [2.0, 3.0], [1.0, 1.0, 1.0], [0.1, 0.1, 0.1], 5.0, 1e-4],
    
    # TC 6: Unbounded Case (h_v >= 0 check) -> Expected Obj: -1.0 (Sentinel Value for Unbounded)
    # Max x1 + x2; -x1 + x2 <= 1 (Only one upper bound, objective pushes out indefinitely)
    [1, 2, 0.95, 200, [-1.0, 1.0], [1.0], [1.0, 1.0], [0.1, 0.1], -1.0, 1e-9],
    
    # TC 7: High Coefficient Scaling -> Expected Obj: 1000.0
    # Max 1000x1 + 0x2; x1 <= 1, x2 <= 100 (Opt: x1=1, x2=0)
    [2, 2, 0.95, 300, [1.0, 0.0, 0.0, 1.0], [1.0, 100.0], [1000.0, 0.0], [0.1, 0.1], 1000.0, 1e-3],
    
    # TC 8: Highly Constrained/Tight (Tests the numerical EPSILON fix) -> Expected Obj: 7.0
    # Max 5x1 + 2x2; 3x1+x2 <= 7, x1+2x2 <= 3 (Opt: x1=11/5, x2=2/5) -> Fails LP intuition.
    # New TC8: Max 2x1 + x2; x1+x2 <= 2, 3x1+x2 <= 3 (Opt: x1=0.5, x2=1.5)
    [2, 2, 0.95, 300, [1.0, 1.0, 3.0, 1.0], [2.0, 3.0], [2.0, 1.0], [0.1, 0.1], 2.5, 1e-4],

    # TC 9: Near Optimal Initial Point (Tests fast convergence) -> Expected Obj: 7.0
    # Max 3x1 + x2; x1+x2 <= 3, 2x1+x2 <= 5 (Opt: x1=3, x2=0)
    [2, 2, 0.95, 50, [1.0, 1.0, 2.0, 1.0], [3.0, 5.0], [3.0, 1.0], [2.0, 0.5], 9.0, 1e-4], # Opt: x1=3, x2=0. Obj=9 (Constraint 1 active)

    # TC 10: Highly Asymmetric/Wide Feasible Region -> Expected Obj: 100.0
    # Max 100x1; x1 <= 1, 0.1x1 + 10x2 <= 10 (Opt: x1=1, x2=0.99)
    [2, 2, 0.95, 400, [1.0, 0.0, 0.1, 10.0], [1.0, 10.0], [100.0, 1.0], [0.5, 0.5], 100.99, 1e-3],
]

def write_test_case(test_id, data):
    """Writes a single test case to a file."""
    filepath = os.path.join(TEST_DIR, f"test_case_{test_id}.txt")
    m, n, gamma, max_iter, A, b, c, x0, expected_obj, tolerance = data

    with open(filepath, 'w') as f:
        # Line 1: Type
        f.write("max\n")
        # Line 2: m n
        f.write(f"{m} {n}\n")
        # Line 3: gamma
        f.write(f"{gamma}\n")
        # Line 4: Max_Iterations
        f.write(f"{max_iter}\n")
        
        # A matrix
        a_idx = 0
        for i in range(m):
            row = ' '.join(map(str, A[a_idx:a_idx + n]))
            f.write(f"{row}\n")
            a_idx += n
            
        # B vector
        f.write(' '.join(map(str, b)) + "\n")
        
        # C vector
        f.write(' '.join(map(str, c)) + "\n")
        
        # X0 vector
        f.write(' '.join(map(str, x0)) + "\n")

    # Write the expected output to a separate meta file for the C++ test runner
    meta_filepath = os.path.join(TEST_DIR, f"test_case_{test_id}.meta")
    with open(meta_filepath, 'w') as f:
        f.write(f"{expected_obj}\n")
        f.write(f"{tolerance}\n")
        # Note: We skip checking the full X vector as Karmarkar may converge to different points 
        # in the optimal set; checking the objective value is sufficient and more robust.


def generate_test_cases():
    """Main function to generate the test directory and files."""
    if not os.path.exists(TEST_DIR):
        os.makedirs(TEST_DIR)
        print(f"Created directory: {TEST_DIR}")
    
    print("Generating 10 test case files...")
    
    for i, test_data in enumerate(TEST_CASES):
        test_id = i + 1
        write_test_case(test_id, test_data)
        print(f"  - Generated test_case_{test_id}.txt and metadata.")

    print("\nGeneration complete. Run 'testing.cpp' now.")

if __name__ == "__main__":
    generate_test_cases()
