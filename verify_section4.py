"""verify_section4.py -- check every numerical statement of Section 4 of the paper
against the tables written by the Maple experiments (tables/*.tex).

    python3 verify_section4.py
"""
import re, statistics, sys, os

TAB = os.path.join(os.path.dirname(os.path.abspath(__file__)), "tables")
npass = nfail = 0


def num(cell):
    cell = cell.strip()
    m = re.search(r"\$?(-?[\d.]+)\\times 10\^\{(-?\d+)\}\$?", cell)
    if m:
        return float(m.group(1)) * 10 ** int(m.group(2))
    m = re.search(r"\$?(-?[\d.]+)\$?", cell)
    return float(m.group(1)) if m else None


def rows(name):
    out = []
    for line in open(os.path.join(TAB, name)):
        line = line.strip()
        if line.endswith("\\\\") and "&" in line and not line.startswith("&") and "$\\ell$" not in line \
                and "example" not in line and "passes" not in line and "$\\delta$" not in line:
            out.append([c.strip() for c in line[:-2].split("&")])
    return out


def check(cid, ok, msg):
    global npass, nfail
    if ok:
        npass += 1
    else:
        nfail += 1
    print("[%s] %-5s %s" % ("PASS" if ok else "FAIL", cid, msg))


# ---------------------------------------------------------------- Table 1
T1 = rows("tab_examples.tex")
def t1(ex, eps, meth):
    for r in T1:
        if r[0] == ex and num(r[1]) == eps and r[2].startswith(meth):
            return [num(c) if "--" not in c else None for c in r[3:]]
    raise KeyError((ex, eps, meth))
check("4.1a", t1("Ex.~1", 0, "standard")[1] < 1e-14, "Ex.1 exact: GCD error at machine-epsilon level (%.1e)" % t1("Ex.~1", 0, "standard")[1])
r = [t1("Ex.~1", e, "standard")[1] / e for e in (1e-10, 1e-8)]
check("4.1b", all(5 <= x <= 7 for x in r), "Ex.1 noisy: GCD error / eps = %s (text: about 6)" % ["%.1f" % x for x in r])
check("4.1c", all(t1("Ex.~1", e, "B")[1] > t1("Ex.~1", e, "standard")[1] for e in (1e-10, 1e-8)),
      "Ex.1 noisy: Bezout lifting less accurate than proposed")
check("4.1d", all(t1("Ex.~2", e, m)[1] >= 0.99 for e in (0, 1e-10, 1e-8) for m in ("standard", "B")),
      "Ex.2: standard and Bezout lifting fail (error ~1) for every eps")
check("4.1e", abs(t1("Ex.~2", 0, "two-dim")[1] - 7.6e-12) < 0.05e-12, "Ex.2 exact: two-dim GCD error 7.6e-12")
r = [t1("Ex.~2", e, "two-dim")[1] / e for e in (1e-10, 1e-8)]
check("4.1f", all(1e3 <= x <= 1e4 for x in r), "Ex.2 noisy: two-dim error / eps = %s (text: 1e3--1e4)" % ["%.1e" % x for x in r])

# ---------------------------------------------------------------- Table 2
T2 = rows("tab_illcond.tex")
def t2(eps, dl):
    for r in T2:
        if num(r[0]) == eps and abs(num(r[1]) - dl) < 1e-12:
            return dict(ratio=num(r[2]), std=num(r[4]), ill=num(r[5]), ada=num(r[6]), which=r[6].split("(")[-1][:3], bez=num(r[7]))
deltas = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6]
check("4.2a", t2(0, 1e-1)["std"] < 1e-6 and t2(0, 1e-1)["bez"] < 1e-6 and all(t2(0, d)["std"] > 0.1 and t2(0, d)["bez"] > 0.1 for d in deltas[1:]),
      "eps=0: standard and Bezout accurate only for delta=1e-1")
check("4.2b", abs(t2(0, 1e-1)["ratio"] - 1e-2) < 1e-4, "delta=1e-1: sigma_{K-1}/sigma_1 = 1e-2")
check("4.2c", t2(1e-10, 1e-1)["std"] > 1e-3 and t2(1e-10, 1e-1)["bez"] > 1e-3, "eps=1e-10, delta=1e-1: standard / Bezout lose many digits (%.1e, %.1e)" % (t2(1e-10, 1e-1)["std"], t2(1e-10, 1e-1)["bez"]))
check("4.2d", all(t2(0, d)["ill"] < 1e-5 for d in deltas[1:]) and t2(0, 1e-1)["ill"] > 1e-2, "eps=0: two-dim accurate for delta<=1e-2, fails for 1e-1")
check("4.2e", all(t2(e, d)["ada"] == min(t2(e, d)["std"], t2(e, d)["ill"]) for e in (0, 1e-10) for d in deltas), "adaptive choice = better of the two in all rows")
check("4.2f", all(1e-7 <= t2(1e-10, d)["ill"] <= 1e-6 for d in deltas[2:]), "eps=1e-10: two-dim error in [1e-7, 1e-6] for delta<=1e-3")
T3 = [[num(c) for c in r] for r in rows("tab_passes.tex")]
check("4.2g", abs(T3[5][1] - 7.6e-12) < 0.05e-12 and abs(T3[5][3] - 1.1e-7) < 0.05e-7, "passes: stagnation at 7.6e-12 (eps=0) and 1.1e-7 (eps=1e-10)")
red = [T3[i + 1][1] / T3[i][1] for i in range(3)]
check("4.2h", all(5e-4 < x < 3e-3 for x in red), "passes: each pass gains about three digits %s" % ["%.1e" % x for x in red])
check("4.2i", T3[4][1] < 1e-11 and T3[3][3] < 1.2e-7, "four or five passes are sufficient")

# ---------------------------------------------------------------- Table 4
T4 = rows("tab_dense.tex")
exact = [r for r in T4 if num(r[2]) == 0]
noisy = [r for r in T4 if num(r[2]) == 1e-8]
check("4.3a", all(r[3] == "5/5" for r in exact), "dense exact: D correct in all cases")
check("4.3b", sum(int(r[3].split("/")[0]) for r in noisy) == 16, "dense eps=1e-8: D correct 16/20")
check("4.3c", all(num(r[4]) > 0.5 for r in T4), "fixed normalisation: error of order one")
check("4.3d", all(num(r[5]) < 1e-13 for r in exact), "dense exact: reconstruction to almost machine precision")
me = [num(r[5]) for r in noisy]
check("4.3e", abs(min(me) - 8.4e-8) < 0.5e-8 and abs(max(me) - 2.6e-5) < 0.5e-5, "dense eps=1e-8: median errors 8e-8 -- 3e-5 (%.1e .. %.1e)" % (min(me), max(me)))

# ---------------------------------------------------------------- Table 5
T5 = rows("tab_random.tex")
R5 = []
for r in T5:
    size = re.findall(r"\d+", r[3])
    R5.append(dict(l=int(r[0]), d=int(r[1]), eps=num(r[2]), rows=int(size[0]), cols=int(size[1]), td=float(r[4]),
                   ed=num(r[5]), t0=float(r[6]), e0m=num(r[7]), e0x=num(r[8]), esm=num(r[9]), esx=num(r[10])))
check("4.4a", max(r["rows"] for r in R5) == 1540 and max(r["cols"] for r in R5) == 120, "direct matrix up to 1540x120")
slower = [(r["l"], r["d"], r["eps"]) for r in R5 if r["t0"] >= r["td"]]
check("4.4b", all(l == 2 and d == 1 for l, d, e in slower), "proposed faster except smallest problems (slower rows: %s)" % slower)
fac = [r["td"] / r["t0"] for r in R5 if r["l"] == 4 and r["d"] == 2]
check("4.4c", abs(min(fac) - 5.6) < 0.05 and abs(max(fac) - 8.2) < 0.05, "l=4, d=2: speed-up %s (text 5.6--8.2)" % ["%.2f" % x for x in fac])
ex0 = [r for r in R5 if r["eps"] == 0]
check("4.4d", max(r["e0m"] for r in ex0) < 3e-14, "exact: median errors of order 1e-14 or smaller (max median %.1e)" % max(r["e0m"] for r in ex0))
check("4.4e", abs(max(r["e0x"] for r in ex0) - 1.3e-9) < 0.05e-9, "exact: occasional errors up to 1e-9 (%.1e)" % max(r["e0x"] for r in ex0))
nz = [r for r in R5 if r["eps"] > 0]
check("4.4f", all(r["ed"] < 3 * r["eps"] for r in nz), "noisy: direct SVD errors of order eps (max ratio %.1f)" % max(r["ed"] / r["eps"] for r in nz))
rat = [(r["e0m"] / r["eps"], r) for r in nz]
others = [x for x, r in rat if not (r["l"] == 3 and r["d"] == 2 and r["eps"] == 1e-6)]
exc = [r for x, r in rat if r["l"] == 3 and r["d"] == 2 and r["eps"] == 1e-6][0]
check("4.4g", 3 <= min(others) and max(others) <= 65 and abs(exc["e0m"] - 1.2e-3) < 0.05e-3,
      "noisy: median/eps in [%.1f, %.1f] except l=3,d=2,eps=1e-6 (%.1e)" % (min(others), max(others), exc["e0m"]))
big = sorted([(r["l"], r["d"], r["eps"]) for r in R5 if r["e0x"] >= 1e-1])
check("4.4h", big == sorted([(4, 2, 1e-8), (2, 3, 1e-8), (3, 2, 1e-6), (4, 1, 1e-6)]), "rows with max error >= 1e-1 without shift: %s" % big)
pairs = {(r["l"], r["d"], r["eps"]): (r["e0x"], r["esx"]) for r in R5}
want = {(4, 2, 1e-8): (9.9e-1, 6.8e-6), (2, 3, 1e-8): (1.6e-1, 4.4e-6), (3, 2, 1e-6): (1.0, 5.1e-3), (4, 1, 1e-6): (4.2e-1, 4.0e-4)}
ok = all(abs(pairs[k][0] - v[0]) < 0.05 * v[0] and abs(pairs[k][1] - v[1]) < 0.05 * v[1] for k, v in want.items())
check("4.4i", ok and abs(exc["esm"] - 2.8e-5) < 0.05e-5, "shift rule: worst errors as quoted, median 1.2e-3 -> %.1e" % exc["esm"])

print("\n%d PASS, %d FAIL" % (npass, nfail))
sys.exit(1 if nfail else 0)
