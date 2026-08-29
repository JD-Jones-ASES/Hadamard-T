#!/usr/bin/env python3
"""Cert 01 -- the public H(2060), pinned and decoded. See NOTES.md.

Exit 0 iff every check passes.

  [A] the banked artifact data/H2060.txt.gz: gunzip, payload SHA-256 against
      the pin, shape pin (2060 lines x 2060 chars over {+,-}), then
      verify/verify.py on the expansion with the canonical sha in its verdict;
      plus an optional cross-check outside the trust chain (full int64 Gram
      H.H^T = 2060*I in numpy, skipped gracefully when numpy is absent);
  [B] the banked length-515 seed: row sums (-17, 5, 39, 15), sum of squares
      2060 = 4v, aggregate PAF = 2060*delta_0 by an own O(v^2) double loop;
  [C] strict invariance of the normalised seed under the order-6 multiplier
      group <149> = {1, 46, 56, 104, 149, 159} <= U(Z_515), and recovery of
      the raw seed by the banked translations (0, 206, 206, 206);
  [D] the T-matrix readout: under the CRT split Z_515 = Z_5 x Z_103 one 4x4
      Hadamard combination of the raw seed's fibres yields four (0,+-1)
      circulant-generating rows of length 103 with disjoint supports of sizes
      (33, 30, 19, 21) partitioning Z_103 and aggregate PAF = 103*delta_0 by
      an own loop, and every fibre is a sign combination of that ONE quadruple;
  [E] the PLAIN Goethals-Seidel array over the same seed is a DIFFERENT exact
      H(2060), built here, handed to verify/verify.py, sha pinned;
  [F] the x104-TWISTED Goethals-Seidel array over the same seed, relabelled
      back through ORD[g] = 5*(g mod 103) + (g mod 5), reproduces the banked
      artifact byte for byte;
  [G] positive control: the T(t) -> H(20t) template, re-run on an
      independently brute-forced quadruple of length 3, yields a
      verify.py-green H(60).

Standard library only (numpy optional in [A], never required). Exact integer
arithmetic only: the 2060 x 2060 matrices are carried as lists of '+'/'-'
strings, and the only sums taken are integer PAF sums and the small template.

    python certs/01-h2060-artifact/run.py            # everything
    python certs/01-h2060-artifact/run.py --quick    # skips orders > 1200
"""

import argparse
import gzip
import hashlib
import itertools
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.join(ROOT, "tools"))
import assemble as asm  # noqa: E402

DATA = os.path.join(ROOT, "data")
GZ = os.path.join(DATA, "H2060.txt.gz")
SEED_JSON = os.path.join(DATA, "seed515.json")
T103_JSON = os.path.join(DATA, "T103.json")
T3_JSON = os.path.join(DATA, "T3-control.json")

V, N, M = 515, 2060, 103
EXPECT_ROW_SUMS = [-17, 5, 39, 15]
EXPECT_SUPPORTS = [33, 30, 19, 21]
ARTIFACT_SHA = "c7a145d86210740dd3f8ea21ca896a54d6916007a042638f17c8c47f097200f7"
PLAIN_SHA = "510f89b7b423c85da0c7cada52cb0f62e0415d736d2040d6701f66c4e524cf6a"
H60_SHA = "70f979f207bbd1d05fe2e28f7730206d797f44b972ac649466b21e39f53a1b13"

FLIP = str.maketrans("+-", "-+")
check = asm.check


# ------------------------------------------------------------------ arithmetic

def own_paf(a, n):
    """Exact integer PAF vector by an O(n^2) double loop written here."""
    out = []
    for s in range(n):
        t = 0
        for j in range(n):
            t += a[j] * a[(j + s) % n]
        out.append(t)
    return out


def own_aggregate_paf(seqs, n):
    total = [0] * n
    for a in seqs:
        for i, val in enumerate(own_paf(a, n)):
            total[i] += val
    return total


def cyclic_group(g, n):
    """<g> <= U(Z_n), closed by repeated multiplication."""
    S, x = set(), 1
    while True:
        S.add(x)
        x = (x * g) % n
        if x == 1:
            return sorted(S)


# ------------------------------------------- 515x515 blocks as lists of strings

def circ(seq):
    """Plain circulant M[i][j] = seq[(i-j) mod v], rows as strings."""
    v = len(seq)
    r0 = "".join(seq[(-j) % v] for j in range(v))
    return [r0[v - i:] + r0[:v - i] for i in range(v)]


def rev_seq(seq):
    v = len(seq)
    return "".join(seq[(-k) % v] for k in range(v))


def twist(seq, k):
    """Twisted circulant M[i][j] = seq[(i - k*j) mod v]."""
    v = len(seq)
    sig = [(k * j) % v for j in range(v)]
    return ["".join(row[t] for t in sig) for row in circ(seq)]


def neg(rows):           # -M
    return [r.translate(FLIP) for r in rows]


def mulJ(rows):          # M @ J (J = back identity): reverse each row
    return [r[::-1] for r in rows]


def conjJ(rows):         # J @ M @ J: reverse the row order and each row
    return [r[::-1] for r in rows[::-1]]


def relabel(rows, inv):
    """M'[i][j] = M[inv[i]][inv[j]]."""
    out = []
    for i in inv:
        src = rows[i]
        out.append("".join(src[t] for t in inv))
    return out


def write_gs(path, grid, v):
    """Write the 4x4 block grid in canonical serialisation; return its sha."""
    h = hashlib.sha256()
    with open(path, "w", encoding="ascii", newline="\n") as fh:
        for r in range(4):
            b0, b1, b2, b3 = grid[r]
            for p in range(v):
                line = b0[p] + b1[p] + b2[p] + b3[p] + "\n"
                fh.write(line)
                h.update(line.encode("ascii"))
    return h.hexdigest()


def verdict_line(path):
    rc, out, _ = asm.run_verify(path)
    return rc, (out.splitlines() or [""])[-1]


# ---------------------------- the T(t) -> H(20t) plug-in template, read off
# the artifact as pure index arithmetic: a core of order 5t is indexed by
# (c, q) in Z_5 x Z_t and built from a level triple of length-t sequences
# chosen by the Z_5 displacement |delta|; the A-core uses the DIFFERENCE
# c'-c, the B/C/D cores the reflected SUM 4-c-c' (the kron(R5, I_t) factor).

FAM = {"A": ("W", "X", "Y"), "B": ("V", "Z", "U"),
       "C": ("P", "U", "Z"), "D": ("Q", "Y", "X")}
SIG = {"A": (1, 1, 1), "B": (1, 1, 1), "C": (1, -1, 1), "D": (1, -1, 1)}
LEV = (0, 1, 2, 2, 1)                      # |delta| class of delta in Z_5
# the Goethals-Seidel grid: name, star = conjugate by R = kron(R5, R_t), sign
GRID = (("A", 0, 1), ("B", 0, 1), ("C", 0, 1), ("D", 0, 1),
        ("B", 0, -1), ("A", 0, 1), ("D", 1, -1), ("C", 1, 1),
        ("C", 0, -1), ("D", 1, 1), ("A", 0, 1), ("B", 1, -1),
        ("D", 0, -1), ("C", 1, -1), ("B", 1, 1), ("A", 0, 1))


def template(seq, t):
    """Seeds of length t derived from a T-quadruple -> a 20t x 20t matrix."""
    def core(k, c, q, cp, qp):
        dc = (cp - c) % 5 if k == "A" else (4 - cp - c) % 5
        i = LEV[dc]
        return SIG[k][i] * seq[FAM[k][i]][(qp - q) % t]

    def blk(k, star, c, q, cp, qp):
        if star:                                    # R . blk . R
            c, q, cp, qp = (4 - c) % 5, t - 1 - q, (4 - cp) % 5, t - 1 - qp
        if k != "A":                                # blk = core @ R
            cp, qp = (4 - cp) % 5, t - 1 - qp
        return core(k, c, q, cp, qp)

    n5 = 5 * t
    H = [[0] * (20 * t) for _ in range(20 * t)]
    for I in range(4):
        for J in range(4):
            k, star, sg = GRID[4 * I + J]
            for c in range(5):
                for q in range(t):
                    row = H[I * n5 + c * t + q]
                    for cp in range(5):
                        for qp in range(t):
                            row[J * n5 + cp * t + qp] = sg * blk(
                                k, star, c, q, cp, qp)
    return H


# ------------------------------------------------------------------------ main

def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--quick", action="store_true",
                    help="skip every matrix of order > %d"
                         % asm.QUICK_MAX_ORDER)
    args = ap.parse_args()
    big = args.quick and N > asm.QUICK_MAX_ORDER
    t0 = time.time()
    print("=" * 74)
    print("CERT 01 -- the public H(2060): pinned, and decoded to a T-matrix "
          "quadruple of\n           order 103")
    print("           trust chain: verify/verify.py")
    print("=" * 74)
    print("python %s" % sys.version.split()[0])
    seed = json.load(open(SEED_JSON, "r", encoding="ascii"))
    tq = json.load(open(T103_JSON, "r", encoding="ascii"))
    ctl = json.load(open(T3_JSON, "r", encoding="ascii"))

    # ------------------------------------------------ [A] the banked artifact
    print("\n[A] the banked artifact data/H2060.txt.gz")
    data = gzip.open(GZ, "rb").read()
    sha = hashlib.sha256(data).hexdigest()
    print("      payload sha256: %s" % sha)
    check("payload sha256 equals the pin %s..." % ARTIFACT_SHA[:16],
          sha == ARTIFACT_SHA)
    text = data.decode("ascii")
    lines = [ln for ln in text.split("\n") if ln]
    check("shape pin: %d lines x %d chars over the alphabet {+,-}" % (N, N),
          len(lines) == N
          and all(len(ln) == N and set(ln) <= {"+", "-"} for ln in lines))
    if big:
        print("      [--] skipped under --quick: verify.py on order %d" % N)
    else:
        path = asm.scratch_path("banked-H2060")
        try:
            with open(path, "w", encoding="ascii", newline="\n") as fh:
                fh.write(text)
            rc, line = verdict_line(path)
        finally:
            os.unlink(path)
        print("      verify.py: %s" % line[:118])
        check("verify.py exit 0, HADAMARD order=2060, canonical sha in the "
              "verdict", rc == 0 and "HADAMARD" in line
              and "order=2060 " in line and ARTIFACT_SHA in line)
        try:
            import numpy as np
        except ImportError:
            print("      [--] numpy cross-check skipped (numpy absent; it is "
                  "not required and is outside the trust chain)")
        else:
            H = np.array([[1 if c == "+" else -1 for c in ln] for ln in lines],
                         dtype=np.int64)
            G = H @ H.T
            ok = ((np.diag(G) == N).all()
                  and (G - np.diag(np.diag(G)) == 0).all())
            check("cross-check outside the trust chain: full int64 Gram "
                  "H.H^T = %d*I in numpy" % N, bool(ok))

    # ---------------------------------------------------- [B] the banked seed
    print("\n[B] the banked length-515 seed")
    ns = seed["normalised_seed"]
    check("four banked strings of length 515 over {+,-}",
          len(ns) == 4 and all(len(s) == V and set(s) <= {"+", "-"}
                               for s in ns))
    nrm = [[1 if ch == "+" else -1 for ch in s] for s in ns]
    sums = [sum(a) for a in nrm]
    print("      row sums recomputed: %s   sum of squares %d   4v = %d"
          % (sums, sum(x * x for x in sums), 4 * V))
    check("row sums = %s (equal to the banked field)" % EXPECT_ROW_SUMS,
          sums == EXPECT_ROW_SUMS == list(seed["row_sums"]))
    check("sum of squares of row sums = 4v = %d (four-square condition)"
          % (4 * V), sum(x * x for x in sums) == 4 * V)
    agg = own_aggregate_paf(nrm, V)
    print("      aggregate PAF (own O(v^2) loop): agg[0] = %d, "
          "max|agg[s>0]| = %d" % (agg[0], max(abs(x) for x in agg[1:])))
    check("sum_q PAF_q = 2060*delta_0 -- the Goethals-Seidel condition",
          agg[0] == 4 * V and all(x == 0 for x in agg[1:]))

    # ---------------------------------------- [C] the multiplier tier <149>
    print("\n[C] the multiplier tier <149> and the raw seed")
    G149 = cyclic_group(149, V)
    print("      <149> mod 515 = %s (order %d)" % (G149, len(G149)))
    check("<149> = {1,46,56,104,149,159}, order 6",
          G149 == [1, 46, 56, 104, 149, 159] and len(G149) == 6)
    strict = all(nrm[i][u] == nrm[i][(g * u) % V]
                 for i in range(4) for g in G149 for u in range(V))
    check("all four banked sequences STRICTLY <149>-invariant (a[u] = a[g u])",
          strict)
    e = list(seed["normalising_shifts"])
    check("banked normalising shifts = [0, 206, 206, 206]",
          e == [0, 206, 206, 206])
    rawseq = ["".join(ns[i][(u + e[i]) % V] for u in range(V))
              for i in range(4)]
    rawv = [[1 if ch == "+" else -1 for ch in s] for s in rawseq]
    ragg = own_aggregate_paf(rawv, V)
    check("un-normalised (raw) seed keeps the row sums and the PAF identity",
          [sum(a) for a in rawv] == EXPECT_ROW_SUMS
          and ragg[0] == 4 * V and all(x == 0 for x in ragg[1:]))
    G3 = cyclic_group(46, V)
    raw3 = all(rawv[i][u] == rawv[i][(g * u) % V]
               for i in range(4) for g in G3 for u in range(V))
    raw6 = all(rawv[i][u] == rawv[i][(g * u) % V]
               for i in range(4) for g in G149 for u in range(V))
    print("      <46> mod 515 = %s (order %d)" % (G3, len(G3)))
    check("the raw seed is strictly <46>-invariant (order 3) but NOT "
          "<149>-invariant -- the tier is reached only after translation",
          raw3 and not raw6)

    # ------------------------------------------------- [D] the T(103) readout
    print("\n[D] the T-matrix quadruple of order 103")
    check("515 = 5 * 103 with gcd(5,103) = 1 (CRT split of Z_515)",
          5 * 103 == V)
    fib = {}
    for i, nm in enumerate("abcd"):
        F = [[0] * M for _ in range(5)]
        for u in range(V):
            F[u % 5][u % 103] = rawv[i][u]
        fib[nm] = F
    pair_a = fib["a"][1] == fib["a"][4] and fib["a"][2] == fib["a"][3]
    pair_r = all(fib[nm][0] == fib[nm][3] and fib[nm][1] == fib[nm][2]
                 for nm in "bcd")
    check("fibre pairing law: F_a[1]=F_a[4], F_a[2]=F_a[3] (untwisted core) "
          "and F_s[0]=F_s[3], F_s[1]=F_s[2] for s in {b,c,d} (twisted cores)",
          pair_a and pair_r)
    H4 = [[1, 1, 1, 1], [1, 1, -1, -1], [1, -1, 1, -1], [1, -1, -1, 1]]
    Q = [fib["a"][1], fib["a"][2], fib["b"][0], fib["b"][1]]
    raw4 = [[sum(H4[i][k] * Q[k][q] for k in range(4)) for q in range(M)]
            for i in range(4)]
    div4 = all(x % 4 == 0 for row in raw4 for x in row)
    check("the 4x4 Hadamard combination of (F_a[1],F_a[2],F_b[0],F_b[1]) is "
          "divisible by 4 throughout", div4)
    T = [[x // 4 for x in row] for row in raw4] if div4 else [[0] * M] * 4
    check("four circulant-generating rows over the alphabet {0,+1,-1}",
          set(x for row in T for x in row) <= {0, 1, -1})
    supp = [set(q for q in range(M) if T[i][q]) for i in range(4)]
    sizes = [len(s) for s in supp]
    print("      support sizes %s (sum %d)" % (sizes, sum(sizes)))
    check("supports pairwise disjoint and covering Z_103, sizes %s"
          % EXPECT_SUPPORTS,
          sizes == EXPECT_SUPPORTS and sum(sizes) == M
          and len(set().union(*supp)) == M
          and all(sum(1 for i in range(4) if T[i][q]) == 1 for q in range(M)))
    tp = own_aggregate_paf(T, M)
    print("      sum_i PAF(T_i) (own loop): [0] = %d, max|[s>0]| = %d"
          % (tp[0], max(abs(x) for x in tp[1:])))
    check("sum_i PAF(T_i) = 103*delta_0 -- the T-matrix (PERIODIC) axiom",
          tp[0] == M and all(x == 0 for x in tp[1:]))
    rs = [sum(r) for r in T]
    check("derived identity sum_i (row sum)^2 = %d = t"
          % sum(v * v for v in rs), sum(v * v for v in rs) == M)
    nagg = asm.aggregate_npaf(T)
    nz = sum(1 for s in range(1, M) if nagg[s])
    print("      [--] recorded, NOT a failure: the aggregate NON-PERIODIC "
          "autocorrelation is\n           nonzero at %d of the %d positive "
          "lags. This object is a T-MATRIX\n           quadruple, which is "
          "exactly the Cooper-Wallis hypothesis. It is NOT an\n           "
          "aperiodic T-sequence of length 103, and none is claimed."
          % (nz, M - 1))
    check("the NPAF peak still reads t (%d)" % nagg[0], nagg[0] == M)
    bank_T = [[tq["signs"][q] if tq["classes"][q] == i else 0
               for q in range(M)] for i in range(4)]
    check("derived quadruple == the banked data/T103.json quadruple",
          T == bank_T)
    combos = {}
    for bits in range(16):
        ev = tuple(1 if (bits >> k) & 1 == 0 else -1 for k in range(4))
        combos[tuple(sum(ev[k] * T[k][q] for k in range(4))
                     for q in range(M))] = ev
    fl = [(nm, p) for nm in "abcd" for p in range(5)]
    got = [combos.get(tuple(fib[nm][p])) for nm, p in fl]
    dis = set(g for g in got if g is not None)
    upto = set(frozenset([g, tuple(-x for x in g)]) for g in dis)
    print("      the 20 fibres realise %d distinct sign vectors, %d classes "
          "under global negation" % (len(dis), len(upto)))
    check("every one of the 20 CRT fibres of the raw seed is a sign "
          "combination sum_i eps_i T_i of the SINGLE quadruple",
          all(g is not None for g in got))
    check("exactly 10 distinct sub-block circulants = 8 classes up to global "
          "sign -- 824 +-1 cells generated by one length-103 object",
          len(dis) == 10 and len(upto) == 8)
    SLOT = {"W": ("a", 0), "X": ("a", 1), "Y": ("a", 2), "V": ("b", 4),
            "Z": ("b", 0), "U": ("b", 1), "P": ("c", 4), "Q": ("d", 4)}
    eps = {k: combos.get(tuple(fib[nm][p])) for k, (nm, p) in SLOT.items()}
    check("the eight template slots W,V,P,Q,X,Y,Z,U all resolve to sign "
          "vectors", all(v is not None for v in eps.values()))
    print("      eps: %s" % {k: "".join("+" if x > 0 else "-" for x in eps[k])
                             for k in ["W", "V", "P", "Q", "X", "Y", "Z", "U"]})

    # -------------------------------- [E] the plain Goethals-Seidel array
    print("\n[E] the PLAIN Goethals-Seidel array over the same seed")
    a, b, c, d = rawseq
    if big:
        print("      [--] skipped under --quick: order %d" % N)
    else:
        A = circ(a)
        Bj, Cj, Dj = mulJ(circ(b)), mulJ(circ(c)), mulJ(circ(d))
        Bt, Ct, Dt = (mulJ(circ(rev_seq(b))), mulJ(circ(rev_seq(c))),
                      mulJ(circ(rev_seq(d))))
        grid = [[A, Bj, Cj, Dj],
                [neg(Bj), A, Dt, neg(Ct)],
                [neg(Cj), neg(Dt), A, Bt],
                [neg(Dj), Ct, neg(Bt), A]]
        p_plain = asm.scratch_path("plain-H2060")
        try:
            sha_plain = write_gs(p_plain, grid, V)
            rc, line = verdict_line(p_plain)
        finally:
            os.unlink(p_plain)
        print("      verify.py: %s" % line[:118])
        check("plain GS array: verify/verify.py exit 0, HADAMARD order=2060",
              rc == 0 and "HADAMARD" in line and "order=2060 " in line)
        check("plain GS canonical sha256 == pinned %s..." % PLAIN_SHA[:16],
              sha_plain == PLAIN_SHA and PLAIN_SHA in line
              and seed["plain_gs_canonical_sha256"] == PLAIN_SHA)

    # ------------------------ [F] the x104-twisted array is the artifact
    print("\n[F] the x104-TWISTED array, relabelled, IS the banked artifact")
    psi = [(104 * g) % V for g in range(V)]
    check("psi(g) = 104g is an involution of Z_515 with 104 = 1 (mod 103) and "
          "104 = -1 (mod 5)",
          all(psi[psi[g]] == g for g in range(V))
          and 104 % 103 == 1 and 104 % 5 == 4)
    ORD = [5 * (g % 103) + (g % 5) for g in range(V)]
    check("ORD[g] = 5*(g mod 103) + (g mod 5) is a bijection Z_515 -> Z_515",
          sorted(ORD) == list(range(V)))
    if big:
        print("      [--] skipped under --quick: order %d" % N)
    else:
        inv = [0] * V
        for g in range(V):
            inv[ORD[g]] = g
        A_ = circ(a)
        B_, C_, D_ = mulJ(twist(b, 104)), mulJ(twist(c, 104)), mulJ(twist(d, 104))
        X12, X13, X23 = neg(conjJ(D_)), conjJ(C_), neg(conjJ(B_))
        grid = [[A_, B_, C_, D_],
                [neg(B_), A_, X12, X13],
                [neg(C_), neg(X12), A_, X23],
                [neg(D_), neg(X13), neg(X23), A_]]
        grid = [[relabel(blk, inv) for blk in row] for row in grid]
        p_art = asm.scratch_path("twisted-H2060")
        try:
            sha_art = write_gs(p_art, grid, V)
            print("      built canonical sha256: %s" % sha_art)
            check("twisted array canonical sha256 == the banked artifact sha "
                  "%s..." % ARTIFACT_SHA[:16],
                  sha_art == ARTIFACT_SHA
                  and seed["gist_canonical_sha256"] == ARTIFACT_SHA)
            ours = open(p_art, "rb").read()
        finally:
            os.unlink(p_art)
        check("byte-for-byte equal to the decompressed banked artifact "
              "(%d bytes)" % len(data), ours == data)

    # --------------------------------------- [G] template control at t = 3
    print("\n[G] positive control: the T(t) -> H(20t) template at t = 3")
    sols = []
    for cl in itertools.product(range(4), repeat=3):
        for sg in itertools.product((1, -1), repeat=3):
            Tq3 = [[0] * 3 for _ in range(4)]
            for i in range(3):
                Tq3[cl[i]][i] = sg[i]
            p = own_aggregate_paf(Tq3, 3)
            if p[0] == 3 and p[1] == 0 and p[2] == 0:
                sols.append([list(cl), list(sg)])
    print("      brute force over all 4^3 * 2^3 = 512 candidates: "
          "%d quadruples" % len(sols))
    check("the banked length-3 control is among the brute-forced solutions "
          "(%d found, banked field says %d)"
          % (len(sols), ctl["total_quadruples_of_length_3"]),
          [ctl["classes"], ctl["signs"]] in sols
          and len(sols) == ctl["total_quadruples_of_length_3"])
    T3 = [[ctl["signs"][q] if ctl["classes"][q] == i else 0 for q in range(3)]
          for i in range(4)]
    p3 = own_aggregate_paf(T3, 3)
    check("banked control satisfies the axiom sum_i PAF(T_i) = 3*delta_0",
          p3[0] == 3 and all(x == 0 for x in p3[1:]))
    seq3 = {k: [sum(eps[k][i] * T3[i][q] for i in range(4)) for q in range(3)]
            for k in SLOT}
    check("the eight seeds sum_i eps_i T_i are all +-1 at t = 3",
          all(set(v) <= {1, -1} for v in seq3.values()))
    H60 = template(seq3, 3)
    check("template output is 60 x 60 with entries in {+1,-1}",
          len(H60) == 60 and all(len(r) == 60 for r in H60)
          and all(x in (1, -1) for r in H60 for x in r))
    txt = "\n".join("".join("+" if x == 1 else "-" for x in r)
                    for r in H60) + "\n"
    sha60 = hashlib.sha256(txt.encode("ascii")).hexdigest()
    p60 = asm.scratch_path("H60")
    try:
        with open(p60, "w", encoding="ascii", newline="\n") as fh:
            fh.write(txt)
        rc, line = verdict_line(p60)
    finally:
        os.unlink(p60)
    print("      verify.py: %s" % line[:118])
    check("H(60) from the independent length-3 quadruple: verify.py exit 0, "
          "order=60", rc == 0 and "HADAMARD" in line and "order=60 " in line)
    check("H(60) canonical sha256 == pinned %s..." % H60_SHA[:16],
          sha60 == H60_SHA and ctl["h60_canonical_sha256"] == H60_SHA)

    # ---------------------------------------------------------- verdicts
    print("\n[H] verdicts")
    if big:
        print("        --quick: the order-2060 work in [A], [E] and [F] did "
              "not run, so the two\n        verdicts on the matrix itself are "
              "not established by this run. Run without\n        --quick for "
              "them.")
    else:
        print("VERDICT the banked artifact is a Hadamard matrix of order 2060 "
              "at canonical sha\n        %s  |  PROVEN-BY-CERTIFICATE"
              % ARTIFACT_SHA)
        print("VERDICT it is the x104-twisted Goethals-Seidel array over a "
              "<149>-invariant seed on\n        Z_515 generated by ONE "
              "T-matrix quadruple of order 103  |  "
              "PROVEN-BY-CERTIFICATE")
    print("VERDICT the order-103 object satisfies the PERIODIC identity. It is "
          "a T-matrix\n        quadruple, not an aperiodic T-sequence; TS(103) "
          "is not claimed.")

    print("\n" + "=" * 74)
    print("-- elapsed %.1fs --" % (time.time() - t0))
    if asm.FAIL:
        print("CERT 01: FAIL (%d): %s" % (len(asm.FAIL), asm.FAIL))
        return 1
    print("CERT 01 GREEN%s" % (" (--quick)" if args.quick else ""))
    print("=" * 74)
    return 0


if __name__ == "__main__":
    sys.exit(main())
