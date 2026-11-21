#!/bin/sh
# Docs drift check: every case listed by "euler_2D --help" needs a "## <name>" section in
# docs/cases.md, and every --option printed by "euler_2D --help" or "render --help" must appear in
# docs/options.md. Usage: tests/check_docs.sh BIN_DIR (run from the repository root)
BIN=${1:-build/bin}
missing=""

# Case names: the indented lines under --case in the solver's usage text
cases=$("$BIN/euler_2D" --help 2>&1 | sed -n 's/^    \([a-z][a-z0-9_]*\) .*/\1/p')
for c in $cases; do
    grep -q "^## $c\$" docs/cases.md || missing="$missing case:$c"
done

# Option names of both programs (render exits non-zero after printing its usage)
flags=$({ "$BIN/euler_2D" --help; "$BIN/render" --help; } 2>&1 | grep -o -- '--[a-z][a-z-]*' | sort -u)
for f in $flags; do
    grep -Eq -- "$f([^a-z-]|\$)" docs/options.md || missing="$missing option:$f"
done

if [ -z "$cases" ] || [ -z "$flags" ]; then
    echo "FAIL docs: could not read the --help output from $BIN"
    exit 1
fi
if [ -n "$missing" ]; then
    echo "FAIL docs: not documented:$missing"
    echo "  (cases go in docs/cases.md as '## <name>', options in docs/options.md)"
    exit 1
fi
echo "PASS docs: $(echo $cases | wc -w | tr -d ' ') cases and $(echo $flags | wc -w | tr -d ' ') options documented"
