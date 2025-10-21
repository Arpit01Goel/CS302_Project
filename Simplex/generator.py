#!/usr/bin/env python3
import os

# Directory to store test case files
dir_name = "test_cases_2"
os.makedirs(dir_name, exist_ok=True)

# Each tuple: (filename, content)
test_cases = [
    # 1. Unbounded
    ("input_1.txt", """max
1 2
1 1
1 -1 >= 0
"""),

    # 2. Infeasible contradiction
    ("input_2.txt", """max
2 1
1
1 <= 10
1 >= 20
"""),

    # 3. Degenerate (zero RHS)
    ("input_3.txt", """max
2 2
1 1
1 1 <= 0
1 0 <= 0
"""),

    # 4. Numerical extremes
    ("input_4.txt", """max
2 2
1000000000 1
1 0 <= 0.000000001
0 1 <= 1
"""),

    # 5. Multiple optimal solutions
    ("input_5.txt", """max
3 2
1 1
1 1 <= 10
1 0 <= 10
0 1 <= 10
"""),

    # 6. Minimization
    ("input_6.txt", """min
1 2
2 3
1 1 >= 5
"""),

    # 7. Zero objective
    ("input_7.txt", """max
1 1
0
1 <= 100
"""),

    # 8. Zero coefficients, unbounded
    ("input_8.txt", """max
1 1
1
0 <= 0
"""),

    # 9. Infeasible (negative RHS)
    ("input_9.txt", """max
1 1
1
-1 >= 5
"""),

    # 10. Two-phase required
    ("input_10.txt", """max
3 2
1 2
1 1 >= 4
1 1 <= 10
1 0 <= 6
"""),
]

# Write each file
for filename, content in test_cases:
    path = os.path.join(dir_name, filename)
    with open(path, "w") as f:
        f.write(content.strip() + "\n")
    print(f"✅ Created {path}")

print(f"\nAll {len(test_cases)} edge-case test files generated in '{dir_name}/'")
