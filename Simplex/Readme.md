# Simplex Solver — `input.txt` format

This README explains the exact input format the Simplex program expects (so you can create `input.txt` files or the `tests/input_*.txt` files used by the gtest). It also includes examples and notes about limitations.

---

## Format summary (line-by-line)

1. **Line 1 — Objective type**
   A single token: `max` or `min` (case-insensitive).
   Example:

   ```
   max
   ```

2. **Line 2 — Dimensions**
   Two integers: `m n`

   * `m` = number of constraints
   * `n` = number of original decision variables
     Example:

   ```
   3 2
   ```

3. **Line 3 — Objective coefficients**
   `n` numeric coefficients (space-separated) for the objective function: `c1 c2 ... cn`
   These are the coefficients of the original decision variables (not including slack variables). Numbers may be integers or floating point.
   Example:

   ```
   6 5
   ```

4. **Next `m` lines — Constraints**
   Each constraint is one line with:

   ```
   a1 a2 ... an  REL  b
   ```

   * `a1..an` are the `n` coefficients of that constraint (space-separated).
   * `REL` is either `<=` or `>=` (exactly these two; `=` is **not** supported).
   * `b` is the right-hand side numeric value.

   Examples:

   ```
   2 1 <= 10
   1 3 >= 5
   ```

---

## How the parser interprets the input

* The solver **automatically adds one slack variable per constraint**. Internally the tableau has `n + m` columns (original variables then slack variables).
* A `>=` constraint is converted to `<=` by multiplying the entire constraint (coefficients and RHS) by `-1` before adding the slack variable.
* **Objective sign convention used by the solver code**:

  * For `max` problems the solver expects the internal C vector to contain the **negation** of the objective coefficients for the original `n` variables (i.e. `C[i] = -ci` for `i < n`) and zeros for slack columns. The parser constructs this automatically.
  * For `min` problems the parser converts to the equivalent max problem and sets the internal C accordingly.
* Blank lines are ignored when reading constraints (so you can add spacing).

---

## Allowed numeric formats & separators

* Numbers can be integers or decimals (e.g. `10`, `10.0`, `3.5`, `-2`).
* Fields must be separated by spaces or tabs.
* No commas inside numbers (use `3.5` not `3,5`).

---

## Example inputs

### Example A — small LP (used in examples & tests)

Maximize `x1 + x2` subject to `x1 + x2 <= 10`, `x1 <= 6`, `x2 <= 5`.

```
max
3 2
1 1
1 1 <= 10
1 0 <= 6
0 1 <= 5
```

### Example B — 3-variable LP

Maximize `6x1 + 5x2 + 4x3` subject to three `<=` constraints:

```
max
3 3
6 5 4
2 1 1 <= 180
1 3 2 <= 300
2 1 2 <= 240
```

---

## Program behavior & usage

* Default: the program looks for `input.txt` in the current working directory. Some builds accept a filename as `argv[1]`.
* Output: the program prints the tableau(s), the values of the decision variables (including slack variables), and the computed maximum (or minimum converted value).
* The parser builds canonical `A`, `B`, `C` and calls the Simplex routine — you do not need to add slack variables in the input; they are added automatically.

Run (example, if your binary is `simplex_solver`):

```bash
./simplex_solver                # uses input.txt
# or
./simplex_solver tests/case1.txt
```

---

## Limitations & important notes

* **Equality constraints (`=`) are not supported** by the parser. If you need `=`, convert them manually or extend the solver to use artificial variables / two-phase simplex.
* All variables are assumed to have `x >= 0` (non-negativity).
* The solver detects *unbounded* (no positive entries in pivot column) and will stop, but it does **not** implement a full infeasibility detection routine (two-phase) in the basic parser setup.
* Floating point arithmetic: due to numerical precision, small tolerances are used internally. Avoid tiny coefficients that may cause numerical instability.
* If you use the gtest harness included in the repo, input files are expected under `./tests/` and named `input_1.txt` … `input_10.txt` (or update the test file `BASE_PATH`).

---

## Troubleshooting tips

* If tests show `Failed to parse: ./tests/input_1.txt`, check:

  * Are you running the binary from the project root where `tests/` lives? (`pwd` and `ls ./tests`)
  * Is the first token `max` or `min`? (no BOM or stray characters)
  * Does the constraint line have exactly `n` coefficients followed by `<=`/`>=` and a numeric RHS?
* If you hit a segmentation fault, ensure the input produces a non-empty matrix (`m>0` and `n>0`) and that lines follow the exact format above. The current solver includes defensive checks for tiny pivot values and unbounded pivots; if you see crashes, paste the input and the error message.

---

## Next steps / Extras

If you want, I can:

* generate example `tests/input_*.txt` files, or
* modify the parser to accept `=` constraints by implementing a two-phase simplex, or
* provide a `Makefile` or `run_tests.sh` to compile and run the tests.

Tell me which and I will produce it.
