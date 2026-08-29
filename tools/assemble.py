#!/usr/bin/env python3
"""The Cooper-Wallis assembler shared by the certificates in certs/.

Standard library only. Exact integer arithmetic only: matrices are carried as
lists of '+'/'-' strings and streamed to disk, and the only sums taken are
integer autocorrelation sums. This file is a builder, not a verifier. The
trust chain is verify/verify.py alone, which shares no code with this module.

THE PRODUCT. Inputs are

  * four T-matrices of order t: circulant (0, +-1) matrices T_1..T_4 whose
    first rows have pairwise disjoint supports partitioning Z_t (exactly one
    T_i is nonzero at each of the t positions) and satisfy the PERIODIC
    identity  sum_i PAF(T_i)(s) = t*delta_0(s), where
    PAF(x)(s) = sum_q x[q]*x[(q+s) mod n];

  * a Williamson quadruple of order w: symmetric circulant +-1 first rows
    a, b, c, d generating A, B, C, D with sum_q PAF_q(s) = 4w*delta_0(s),
    i.e. A^2 + B^2 + C^2 + D^2 = 4w*I_w.

Step 1 builds four +-1 matrices of order tw, group-developed over Z_t x Z_w
by the Williamson-array sign table:

    X1 =  T1(x)A + T2(x)B + T3(x)C + T4(x)D
    X2 = -T1(x)B + T2(x)A + T3(x)D - T4(x)C
    X3 = -T1(x)C - T2(x)D + T3(x)A + T4(x)B
    X4 = -T1(x)D + T2(x)C - T3(x)B + T4(x)A

with (x) the Kronecker product of circulants and the index on Z_t x Z_w
lexicographic, u = p*w + r. Every entry is +-1 because exactly one T_i is
nonzero at each position. Each column of the sign table is a permutation of
{A, B, C, D}, so the diagonal coefficient of sum_j X_j X_j^T is
A^2 + B^2 + C^2 + D^2 = 4w*I_w; the off-diagonal coefficients telescope to
zero because symmetric circulants commute. Hence

    sum_j X_j X_j^T = (sum_i T_i T_i^T) (x) 4w*I_w = 4tw*I_{tw}.

Step 2 is the Goethals-Seidel array of order 4tw, R the back-identity of
order tw (X R is column reversal of X):

    [   X1      X2 R     X3 R     X4 R  ]
    [ -X2 R      X1    -X4^T R   X3^T R ]
    [ -X3 R    X4^T R     X1    -X2^T R ]
    [ -X4 R   -X3^T R   X2^T R     X1   ]

R = R_t (x) R_w is the back-identity of order tw under the lexicographic
index, and each X_j is group-developed over the abelian group Z_t x Z_w, so
X_j^T R = R X_j -- the only property the array needs beyond
sum_j X_j X_j^T = 4tw*I. build() does not assume that identity: it builds
X_j^T R from its own closed formula and checks it against R X_j on all four
blocks and all tw rows.

This is Cooper and Wallis 1972: T-matrices of order n, plus four +-1 matrices
of order w pairwise satisfying M N^T = N M^T with sum M M^T = 4w*I, give a
Hadamard matrix of order 4nw. Symmetric circulants satisfy the pairwise
condition.

TERMINOLOGY. T-sequences are defined by the NON-PERIODIC autocorrelation
vanishing; T-matrices by sum_i X_i X_i^T = nI with disjoint supports, i.e. by
the PERIODIC identity. The periodic form is what Cooper-Wallis needs and what
this module assumes.
"""

import hashlib
import json
import os
import subprocess
import sys
import tempfile
import time
from operator import mul

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
VERIFY = os.path.join(ROOT, "verify", "verify.py")
W_JSON = os.path.join(DATA, "williamson-quadruples.json")
CP_TABLE4_JSON = os.path.join(DATA, "cp-table4-v2.json")

# The transcription's shape: 195 entries, all at odd n <= 2999.
CP_ENTRY_COUNT = 195
CP_N_MAX = 2999

# Orders above this are skipped under --quick.
QUICK_MAX_ORDER = 1200

# The Williamson-array sign table. K[j][i] = (sign, which Williamson matrix),
# matrices indexed 0,1,2,3 = A,B,C,D. Row j spells out X_{j+1} as above.
K = [
    [(+1, 0), (+1, 1), (+1, 2), (+1, 3)],
    [(-1, 1), (+1, 0), (+1, 3), (-1, 2)],
    [(-1, 2), (-1, 3), (+1, 0), (+1, 1)],
    [(-1, 3), (+1, 2), (-1, 1), (+1, 0)],
]

NEG = str.maketrans("+-", "-+")
FAIL = []


def check(label, ok):
    """Record one check. Failures accumulate in FAIL; the cert exits 1."""
    print("      [%s] %s" % ("ok" if ok else "FAIL", label))
    if not ok:
        FAIL.append(label)
    return ok


# ------------------------------------------------------------ own exact loops


def paf(seq, s):
    """Periodic autocorrelation of one integer sequence at lag s."""
    n = len(seq)
    return sum(seq[q] * seq[(q + s) % n] for q in range(n))


def npaf(seq, s):
    """Non-periodic (aperiodic) autocorrelation at lag s."""
    return sum(seq[q] * seq[q + s] for q in range(len(seq) - s))


def aggregate_paf(quad):
    """[sum_i PAF(quad_i)(s) for s] -- own O(4n^2) exact-integer loop."""
    n = len(quad[0])
    return [sum(paf(r, s) for r in quad) for s in range(n)]


def aggregate_npaf(quad):
    """[sum_i NPAF(quad_i)(s) for s] -- non-periodic, for the record."""
    n = len(quad[0])
    return [sum(npaf(r, s) for r in quad) for s in range(n)]


def group_aggregate_paf(x, t, w):
    """sum_j sum_u x_j[u] x_j[u+d] over the group Z_t x Z_w, exact ints.

    Returns (peak, offpeak_all_zero).
    """
    total = [[0] * w for _ in range(t)]
    for j in range(4):
        xj = x[j]
        for dr in range(w):
            rot = [r[dr:] + r[:dr] for r in xj]
            col = total
            for dp in range(t):
                acc = 0
                for p in range(t):
                    acc += sum(map(mul, xj[p], rot[(p + dp) % t]))
                col[dp][dr] += acc
    peak = total[0][0]
    total[0][0] = 0
    return peak, not any(v for row in total for v in row)


# ------------------------------------------------------------------- loading


def load_T(path):
    """(classes, signs, T as four rows of {0,+-1}, t, the raw record)."""
    d = json.load(open(path, encoding="ascii"))
    t, cls, sgn = d["n"], d["classes"], d["signs"]
    if not (len(cls) == len(sgn) == t):
        raise SystemExit("%s: length mismatch" % path)
    Tq = [[0] * t for _ in range(4)]
    for q in range(t):
        Tq[cls[q]][q] = sgn[q]
    return cls, sgn, Tq, t, d


def load_W(bank, w):
    """Four +-1 first rows of length w."""
    rows = bank["quadruples"][str(w)]["rows"]
    return [[1 if c == "+" else -1 for c in r] for r in rows]


def load_bank(*extra):
    """data/williamson-quadruples.json, merged with the named extra files."""
    bank = json.load(open(W_JSON, encoding="ascii"))
    for name in extra:
        more = json.load(open(os.path.join(DATA, name), encoding="ascii"))
        for k, v in more["quadruples"].items():
            bank["quadruples"][k] = v
    return bank


def load_cp_table4():
    """data/cp-table4-v2.json -> ({n: m}, provenance).

    The transcribed Cati-Pasechnik Table 4. An entry n(m) means m is the least
    exponent for which a Hadamard matrix of order 2^m*n is recorded known, so
    H(4n) is NOT recorded known while the entry stands. An odd n <= 2999 that
    is ABSENT from the map has m(n) = 2, i.e. H(4n) was already known. The
    file's own entry_count and range are re-checked against the map it carries.
    """
    d = json.load(open(CP_TABLE4_JSON, encoding="ascii"))
    prov = d["provenance"]
    table = {int(k): int(v) for k, v in d["entries"].items()}
    if not (len(table) == prov["entry_count"] == CP_ENTRY_COUNT):
        raise SystemExit("data/cp-table4-v2.json: %d entries, expected %d"
                         % (len(table), CP_ENTRY_COUNT))
    bad = sorted(n for n, m in table.items()
                 if n % 2 == 0 or n > prov["range"]["n_max"] or m < 3)
    if bad:
        raise SystemExit("data/cp-table4-v2.json: malformed entries %s" % bad)
    return table, prov


def cp_open_w(table, t, w_max=63):
    """The odd w <= w_max for which n = t*w carries a Table 4 entry."""
    return sorted(w for w in range(1, w_max + 1, 2)
                  if t * w <= CP_N_MAX and t * w in table)


def cp_agrees(table, static):
    """True iff every (n, m) in the cert's static reading matches the file."""
    return all(table.get(n) == m for n, m in static.items())


# ------------------------------------------------------------- the product


def make_x(Tq, W, t, w):
    """x[j][p][r] -- the four group-developed first rows over Z_t x Z_w.

    x_j[(p,r)] = sum_i K[j][i].sign * T_i[p] * W[K[j][i].idx][r]
    """
    x = []
    for j in range(4):
        xj = []
        for p in range(t):
            row = [0] * w
            for i in range(4):
                sg, mi = K[j][i]
                c = sg * Tq[i][p]
                if c:
                    wm = W[mi]
                    for r in range(w):
                        row[r] += c * wm[r]
            xj.append(row)
        x.append(xj)
    return x


def row_tables(xj, t, w):
    """Per-offset column tables for the developed matrix and its R-transpose.

    COL[r][a]  = the length-w string of x_j[a][(s-r) mod w], s = 0..w-1
    COLB[r][a] = the length-w string of x_j[a][(r+s-w+1) mod w], s = 0..w-1

    Then, with u = (p, r) and v = (q, s) and m = t*w:
      X_j[u][v]        = x_j[(q-p) mod t][(s-r) mod w]
                       = ''.join(COL[r] rotated to start at (-p) mod t)
      (X_j^T R)[u][v]  = X_j[m-1-v][u] = x_j[(p+q+1) mod t][(r+s+1) mod w]
                       = ''.join(COLB[r] rotated to start at (p+1) mod t)
    """
    sym = "+-"
    COL, COLB = [], []
    for r in range(w):
        COL.append([("".join(sym[0] if xj[a][(s - r) % w] > 0 else sym[1]
                             for s in range(w))) for a in range(t)])
        COLB.append([("".join(sym[0] if xj[a][(r + s - w + 1) % w] > 0
                              else sym[1]
                              for s in range(w))) for a in range(t)])
    return COL, COLB


def developed_rows(xj, t, w):
    """(Xrows, XtRrows): each a list of m strings of length m, u = p*w + r."""
    COL, COLB = row_tables(xj, t, w)
    Xr, Xt = [], []
    for p in range(t):
        k = (-p) % t
        k2 = (p + 1) % t
        for r in range(w):
            c = COL[r]
            Xr.append("".join(c[k:] + c[:k]))
            cb = COLB[r]
            Xt.append("".join(cb[k2:] + cb[:k2]))
    return Xr, Xt


def write_gs_array(rows, path):
    """Stream the Goethals-Seidel array to path; return the SHA-256 written.

    The bytes are exactly verify.py's canonical serialisation, so the returned
    digest is exactly the canonical_sha256 verify.py reports for the file.
    """
    X1r, X2r, X3r, X4r = (rows[j][0] for j in range(4))
    X2t, X3t, X4t = (rows[j][1] for j in (1, 2, 3))
    m = len(X1r)
    h = hashlib.sha256()
    with open(path, "wb") as fh:
        for band in range(4):
            for u in range(m):
                A = X1r[u]
                BR = X2r[u][::-1]
                CR = X3r[u][::-1]
                DR = X4r[u][::-1]
                if band == 0:
                    line = A + BR + CR + DR
                elif band == 1:
                    line = (BR.translate(NEG) + A
                            + X4t[u].translate(NEG) + X3t[u])
                elif band == 2:
                    line = (CR.translate(NEG) + X4t[u] + A
                            + X2t[u].translate(NEG))
                else:
                    line = (DR.translate(NEG) + X3t[u].translate(NEG)
                            + X2t[u] + A)
                b = (line + "\n").encode("ascii")
                h.update(b)
                fh.write(b)
    return h.hexdigest()


def scratch_path(tag):
    """A unique scratch filename. Callers must delete it.

    Written under out/ when out/ exists at the repo root, otherwise under the
    system temp directory. The name is unique per process and per call, so
    concurrent runs cannot delete each other's matrix.
    """
    outdir = os.path.join(ROOT, "out")
    d = outdir if os.path.isdir(outdir) else None
    fd, path = tempfile.mkstemp(prefix="hadamard-t-%s-" % tag, suffix=".txt",
                                dir=d)
    os.close(fd)
    return path


def run_verify(path):
    """verify/verify.py in a subprocess. Returns (returncode, stdout, secs)."""
    t0 = time.time()
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1")
    p = subprocess.run([sys.executable, VERIFY, path],
                       capture_output=True, text=True, env=env)
    return p.returncode, p.stdout.strip(), time.time() - t0


# ---------------------------------------------------------------- one build


def build(Tq, W, t, w, tag, verbose=True):
    """Build H(4tw), verify it, delete it. Returns (rc, verdict, sha, secs)."""
    n = 4 * t * w
    m = t * w
    t0 = time.time()
    x = make_x(Tq, W, t, w)
    pm1 = all(v in (1, -1) for xj in x for row in xj for v in row)
    peak, offpeak_zero = group_aggregate_paf(x, t, w)
    if verbose:
        check("%s: X_j entries all +-1" % tag, pm1)
        check("%s: aggregate group-PAF = %d*delta_0 (target 4tw = %d)"
              % (tag, peak, n), peak == n and offpeak_zero)
    if not (pm1 and peak == n and offpeak_zero):
        return 1, "PRE-CHECK FAILED", None, time.time() - t0
    rows = [developed_rows(x[j], t, w) for j in range(4)]
    ident = all(rows[j][1][u] == rows[j][0][m - 1 - u]
                for j in range(4) for u in range(m))
    if verbose:
        check("%s: X_j^T R = R X_j for all four blocks, all %d rows "
              "(the only property the GS array needs beyond 4-complementarity)"
              % (tag, m), ident)
    if not ident:
        return 1, "GS PRECONDITION FAILED", None, time.time() - t0
    path = scratch_path("H%d" % n)
    try:
        sha = write_gs_array(rows, path)
        del rows, x
        size = os.path.getsize(path)
        rc, verdict, vsecs = run_verify(path)
    finally:
        os.unlink(path)
    if verbose:
        print("      built %d x %d (%.1f MB) | verify.py exit=%d in %.1fs"
              % (n, n, size / 1e6, rc, vsecs))
    return rc, verdict, sha, time.time() - t0


def build_and_pin(Tq, W, t, w, order, pin, tag, verbose=True):
    """build() plus the two standard acceptance checks. Returns (ok, sha)."""
    rc, verdict, sha, secs = build(Tq, W, t, w, tag, verbose=verbose)
    print("      %s" % verdict[:170])
    ok1 = check("H(%d): verify.py exit 0, order=%d" % (order, order),
                rc == 0 and "HADAMARD" in verdict
                and "order=%d " % order in verdict)
    ok2 = check("H(%d): canonical sha256 pinned three ways (own stream, "
                "verify.py's verdict field, the constant in run.py) = %s"
                % (order, pin),
                sha == pin and ("canonical_sha256=" + pin) in verdict)
    return ok1 and ok2, sha


def tamper_reject(Tq, W, t, w, order):
    """Flip one character of a green H(4tw); verify.py must exit 1."""
    x = make_x(Tq, W, t, w)
    rows = [developed_rows(x[j], t, w) for j in range(4)]
    path = scratch_path("tamper-H%d" % order)
    try:
        write_gs_array(rows, path)
        with open(path, "r+b") as fh:
            ch = fh.read(1)
            fh.seek(0)
            fh.write(b"-" if ch == b"+" else b"+")
        rc, verdict, _ = run_verify(path)
    finally:
        os.unlink(path)
    return check("flipping ONE character of the green H(%d) makes verify.py "
                 "exit 1" % order, rc == 1 and "FAIL" in verdict)
