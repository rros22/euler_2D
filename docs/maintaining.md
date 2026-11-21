# Maintaining

Conventions that keep the code easy to read and the docs accurate as the project grows.

## Code conventions

- **Language:** C99, compiled with `-Wall -Wextra`; new code should build without warnings.
- **Layout:** one module per directory under `src/` (see [Architecture](architecture.md)). Lower
  layers must not depend on higher ones: `numerics` doesn't include `cases` or `io`.
- **Includes** are relative to `src/`: `#include "numerics/flux.h"`.
- **Names:** functions are prefixed by their module (`box_`, `case_`, `snapshot_`, `solver_`);
  variables say what they hold (`slope`, `w_minus`, `t_next_frame`), except standard physics symbols
  (`rho`, `u`, `p`, `dt`) and loop indices.
- **Style:** 4 spaces, braces on their own line, as in the existing files.
- **Public API:** the `Box` struct and the `box_*` functions are used by cases and tests; change them
  only deliberately.

## Comments

- **File header:** every file starts with a few lines saying what it does and how it fits in.
- **Contracts in headers:** every public function has a short comment in its header stating what it
  does, its inputs (and units or index ranges where relevant), who owns returned memory, and how it
  reports errors.
- **Inside functions:** explain *why*, or what a non-obvious block achieves; don't restate the code.
  Name the method and give a reference for numerics (e.g. "Toro, chapter 14").
- **No history in comments:** don't write "changed X" or "new version"; that belongs in commit
  messages.
- **One source of truth:** user-facing detail (options, case parameters, file formats) lives in
  `docs/`. Code comments point there instead of repeating it. Exception: the frame format is
  described in `src/io/snapshot.h` and summarized in the docs.

## Keeping the docs up to date

`make test` fails if a case or an option is missing from the docs, but it can't check that the
content is right. When you change the code, update the matching page:

| Change | Update |
|---|---|
| New case | `docs/cases.md`: a `## <name>` section with the diagram, box table, states and defaults |
| Changed case defaults or geometry | its section in `docs/cases.md` |
| New or changed option | `docs/options.md` (both the table and any affected explanation) |
| New boundary type or case helper | `docs/new-case.md` |
| New module or moved files | `docs/architecture.md` (code map and "where to change what") |
| Frame format change | `src/io/snapshot.h`, `tools/snapshot.py`, `docs/output-and-rendering.md`; bump `SNAP_VERSION` |
| New test | `docs/testing.md` |

The docs are plain Markdown, one topic per file, and render on GitHub. To publish them as a
searchable website later, add an `mkdocs.yml` listing these files; no rewriting is needed.

## Changing the numerics

1. Make the change and run `make test`. The regression test will fail; that's expected.
2. Confirm the physics tests still pass, and look at a case or two to see the results make sense.
3. Regenerate the goldens with `./build/bin/test_regression --write`.
4. Commit the code and the new goldens together, and say in the message why the results changed.

If the regression test fails and you did *not* mean to change the numerics, find the cause before
going further.

## Repository hygiene

- `build/` and `results/` are not tracked (`.gitignore`); never commit binaries, frames or videos.
- `make clean` removes build output, `make clean-results` removes simulation output.
- Commit messages: a short summary line, then what changed and why.
