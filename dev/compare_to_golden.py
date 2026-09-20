"""Compare a result exported from Studio with the R operator's golden.

    tercenctl data export-csv -s <schemaId> --filePath /tmp/out.csv
    python dev/compare_to_golden.py /tmp/out.csv tests/table1.csv

The point of running this against a *platform* result rather than the in-process one: it checks
the whole path, including how the server read the TSON this operator wrote.
"""
import csv, sys, math

def load(path):
    with open(path) as f:
        rows = list(csv.DictReader(f))
    if not rows:
        sys.exit(f"{path} is empty")
    return rows

got_rows, want_rows = load(sys.argv[1]), load(sys.argv[2])
val_col = next(c for c in got_rows[0] if c.endswith(".asinh"))
want_col = next(c for c in want_rows[0] if c.endswith("asinh"))
key = lambda r: (int(float(r[".ri"])), int(float(r[".ci"])))
got = {key(r): float(r[val_col]) for r in got_rows}
want = {key(r): float(r[want_col]) for r in want_rows}

missing = set(want) - set(got)
extra = set(got) - set(want)
worst, worst_at = 0.0, None
for k, w in want.items():
    if k in got:
        rel = abs(got[k] - w) / max(abs(w), 1.0)
        if rel > worst:
            worst, worst_at = rel, k
print(f"rows: {len(got)} from Tercen, {len(want)} golden; missing {len(missing)}, extra {len(extra)}")
print(f"worst relative difference {worst:.3e}" + (f" at .ri={worst_at[0]} .ci={worst_at[1]}" if worst_at else ""))
ok = not missing and not extra and worst <= 1e-9
print("PARITY OK" if ok else "PARITY FAILED")
sys.exit(0 if ok else 1)
