# Testing

```bash
make test
```

builds everything and runs all tests in a few seconds. Every line starts with `PASS` or `FAIL`, and
the command fails if any test fails.

| Test | File | Checks |
|---|---|---|
| Regression | `tests/test_regression.c` | results are bitwise identical to the stored "golden" results |
| Physics | `tests/test_physics.c` | the solver is correct: exact solutions, conservation, symmetry, order of accuracy |
| Snapshot I/O | `tests/test_io.c` | a frame written and read back gives the same data |
| Render | `tests/render_smoke.sh` | a small run produces the right frames, a correct image and a video |
| Docs | `tests/check_docs.sh` | every case and every option is documented |

## Regression test and the golden files

The files in `tests/golden/` are **stored answers**, not tests. Each holds the complete state (all
8 variables of every cell, at full precision) after a short reference run: the jet geometry at base
10 for 40 steps, once with first order and once with MUSCL-Hancock (van Leer).

The regression test reruns those simulations with the current code and compares every number bit
for bit. It also builds the same setup through the `jet` case, which checks that the case still
matches the original problem.

- **It passes** when your change did not alter the results. This is what you want after
  refactoring, speeding up code or reorganizing files.
- **It fails** when the results changed, even in the last digit. If you didn't mean to change the
  numerics, something is wrong.

The goldens are not analytic solutions: they only show that the results haven't changed. Correctness
is checked by the physics tests.

### Regenerating the goldens

Only after a **deliberate** change to the numerics (a new limiter formula, a different tolerance,
...), and only once the physics tests pass with the new code:

```bash
./build/bin/test_regression --write
```

Commit the new golden files together with the change that required them, and say so in the commit
message.

Results are reproducible across optimization levels because the build uses `-ffp-contract=off`
(no fused multiply-add). Changing compiler flags can still change the last digits; if the goldens
fail after a compiler or flag change only, regenerate them.

## Physics tests

| Check | Compared against |
|---|---|
| Ghost cells | every boundary type fills both ghost layers with the expected values |
| Conservation | total mass and energy in a closed domain stay constant to round-off, for every scheme |
| Symmetry | the same problem along x and along y gives identical results |
| Order of accuracy | a smooth wave advected once around a periodic domain: the error falls at rate ≈ 1 (first order) and ≈ 2 (MUSCL) as the grid is refined, with identical errors along x and y |
| Shock tubes | Sod's problem and Toro's five Riemann problems (Toro, Table 4.1) at CFL 0.4 and 0.9: every scheme must keep ρ > 0 and p > 0, and on Sod each MUSCL limiter must have a smaller error against the exact solution than first order. The errors for all problems are printed for inspection |
## Docs check

`tests/check_docs.sh` reads the case list and the options from `euler_2D --help` and `render --help`
and fails if a case has no `## <name>` section in `docs/cases.md` or an option is missing from
`docs/options.md`. When you add a case or an option, document it and this check passes again.

## Adding a test

Physics checks go in `tests/test_physics.c`: write a `static void test_...(void)` function that calls
`check(ok, name, format, value)` and call it from `main`. Helpers for small setups, running to a
time and comparing results are in `tests/test_util.h`. For a new test program, create
`tests/test_<name>.c` with its own `main`, add `$(BIN_DIR)/test_<name>` to `TESTS` in the Makefile
and a line running it to the `test` target.
