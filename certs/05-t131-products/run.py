#!/usr/bin/env python3
"""Cert 05 -- T-sequences of order 131, and five Hadamard matrices.

See NOTES.md. Exit 0 iff every check passes.

    TT(44)  ->  BS(87,44)  ->  TS(131)  ->  H(4*131*w)

Every arrow is re-derived and re-checked here, from banked data, with exact
integer arithmetic and the standard library only.

  [A] the Turyn-type payloads (data/turyn-type-payloads.json): decoded under
      the published nibble map, then checked against the exact Turyn-type
      identity  N_A + N_B + 2N_C + 2N_D = 0 (i >= 1)  by an own O(n^2) integer
      loop, at all three printed lengths 40, 42, 44. The naive left-aligned
      reading of the last digit is also run, and shown to fail;
  [A'] a from-scratch control that the composition code is not tuned to 44:
      TT(6) and TT(8) are found here by exhaustive search over 2^(4n-1), with
      no external input, and pushed through the same two arrows;
  [B] arrow 1, TT(n) -> BS(2n-1, n) = (C;D, C;-D, A, B): the base-sequence
      identity checked by own loops;
  [C] arrow 2, BS(m,n) -> TS(m+n): the four (0,+-1) sequences built, their
      supports shown to partition Z_131, and BOTH the nonperiodic and the
      periodic aggregate autocorrelations checked. The periodic one is the
      T-matrix axiom, which is what Cooper-Wallis needs. Also recorded: the
      full multiplier group is trivial;
  [C'] the second, independent witness: a K-invariant T-matrix quadruple of
      order 131 from this laboratory's own multiplier-tier search, the tier
      lattice that framed it, and the invariant separation of the two
      witnesses;
  [D] the Williamson quadruples the cash-in assumes;
  [E] the openness of the five target orders, reported from an audited table;
  [F] calibration gates t=131 x w in {1,3} -> H(524), H(1572), known orders;
  [G] the five targets from witness 1, built one at a time, verified by
      verify/verify.py in a subprocess, digest pinned three ways, deleted;
  [G'] the same five orders from witness 2 -- ten matrices in all, in five
      pairs with distinct digests;
  [H] a negative control: one flipped character must make verify.py exit 1.

The trust chain is verify/verify.py alone; the assembler is tools/assemble.py.

    python certs/05-t131-products/run.py           # everything (~4 min)
    python certs/05-t131-products/run.py --quick   # no order > 1200
"""

import argparse
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
PAYLOADS = os.path.join(DATA, "turyn-type-payloads.json")
T131_LIT = os.path.join(DATA, "T131-literature-route.json")
T131_MUL = os.path.join(DATA, "T131-multiplier-tier.json")

T = 131
WS = (1, 3, 7, 11, 13, 17, 19)
check = asm.check

# Canonical SHA-256 of verify.py's '+'/'-' serialisation.
GATE_PINS = {
    (131, 1): (524,
               "97aca328b00f3715fcce2a7e15c8637da09d279e383a5c3d41b6587be60d1205"),
    (131, 3): (1572,
               "4ce087b59bf79e0bfaf1d99c13c9f1a39c397ba24eb518b6e85a6fa9cb908196"),
}
TARGET_PINS = {
    3668: (7,
           "0dee08275687490a49df9c43bfc276ec93f95107365d05e3f3b48cfa31991ddf"),
    5764: (11,
           "b391b9bad54b235f6f1d9c778cd72256760a46f6e04cb26f2193e016d326ba02"),
    6812: (13,
           "eec46e29ad9b7933053c5e62912fb01c4113289e167a8b0ecfc9d9bb287a9d42"),
    8908: (17,
           "0107ba27c6c10f4ff65011f79e4e564a56a7ab4359535a634dc8c4b6bb2b9afe"),
    9956: (19,
           "a14f4b6991fbad498cebc568383119a0f4afc38eafc4be3669bf503619280496"),
}
# The same orders built from the second witness. Every digest must DIFFER from
# the first witness's: two routes, two sets of matrices.
WITNESS2_PINS = {
    524: (1,
          "003385ed9d3e528cea9fdfb85835d84b6cd92169f2062742547f43f48892851f"),
    3668: (7,
           "defab1069bc0af50b5ec85dc913818f1d0763c69cf9c3eaefbedefb4a33cd90a"),
    5764: (11,
           "3cbaf9fbcb28f601ef490e4b631961b4e5d60f463069188a0f830f265881d7cd"),
    6812: (13,
           "4fa5d4d6961eaa699590c5f45264cda2eaf6f42f1f8fd8c60e1ee899138f368b"),
    8908: (17,
           "88c6d5e9cd97eca3e3405779b68fb55f7c1d902d74719a4e8dfffc3979fdc0af"),
    9956: (19,
           "2799b6154c4768b1b2df72b8ab472ddf6565f54c7d4d2d99fde83afbead41387"),
}
# Cati-Pasechnik arXiv:2411.18897v2 Table 4, entries n(m), for the five
# targets: n = 131*w.
CP_TABLE4 = {3668: (917, 3), 5764: (1441, 3), 6812: (1703, 3),
             8908: (2227, 3), 9956: (2489, 3)}
EXPECT_SUPPORTS = [44, 43, 23, 21]
EXPECT_ROWSUMS = [8, 3, -3, -7]


# --------------------------------------------------------- own exact loops
def npaf(x, i):
    """Nonperiodic (aperiodic) autocorrelation at lag i."""
    return sum(x[j] * x[j + i] for j in range(len(x) - i))


def paf(x, i):
    """Periodic autocorrelation at lag i."""
    n = len(x)
    return sum(x[q] * x[(q + i) % n] for q in range(n))


def turyn_residual(A, B, C, D, n):
    """[N_A + N_B + 2N_C + 2N_D](i) for i = 0..n-1. Zero at every i >= 1 and
    6n-2 at i = 0 is the Turyn-type axiom."""
    return [npaf(A, i) + npaf(B, i) + 2 * npaf(C, i) + 2 * npaf(D, i)
            for i in range(n)]


def decode(h, n, right_align_last=True):
    """Column-major hex -> (A, B, C, D) of lengths n, n, n, n-1.

    Hex digit j carries column j of the 4-row array, bits MSB..LSB = the four
    rows, '+' -> 0. The last row has length n-1, so the final digit carries
    only the first three rows and is right-aligned: its high bit is the pad.
    """
    A, B, C, D = [], [], [], []
    for j, ch in enumerate(h):
        v = int(ch, 16)
        if j < n - 1:
            bits = [(v >> k) & 1 for k in (3, 2, 1, 0)]
            A.append(bits[0])
            B.append(bits[1])
            C.append(bits[2])
            D.append(bits[3])
        else:
            order = (2, 1, 0) if right_align_last else (3, 2, 1)
            bits = [(v >> k) & 1 for k in order]
            A.append(bits[0])
            B.append(bits[1])
            C.append(bits[2])
    f = lambda s: [1 - 2 * x for x in s]        # noqa: E731
    return f(A), f(B), f(C), f(D)


def brute_TT(n):
    """A TT(n) found by a complete search of the full 2^(4n-1) space, with no
    external input: every (C,D) half is tabulated by its required partner
    vector -2(N_C + N_D), every (A,B) half is scanned, and the two are joined
    exactly. Complete (a solution exists iff one is returned) at cost
    2^(2n-1) + 2^(2n) instead of 2^(4n-1)."""
    def seqs(L):
        for m in range(1 << L):
            yield [1 - 2 * ((m >> k) & 1) for k in range(L)]
    tab = {}
    for C in seqs(n):
        nc = [npaf(C, i) for i in range(1, n)]
        for D in seqs(n - 1):
            key = tuple(-2 * (nc[i] + npaf(D, i + 1)) for i in range(n - 1))
            if key not in tab:
                tab[key] = (tuple(C), tuple(D))
    for A in seqs(n):
        na = [npaf(A, i) for i in range(1, n)]
        for B in seqs(n):
            key = tuple(na[i] + npaf(B, i + 1) for i in range(n - 1))
            hit = tab.get(key)
            if hit:
                return A, B, list(hit[0]), list(hit[1])
    return None


# ------------------------------------------------------ the two arrows
def TT_to_BS(A, B, C, D):
    """Arrow 1. (C;D, C;-D, A, B) has lengths 2n-1, 2n-1, n, n."""
    return C + D, C + [-v for v in D], A, B


def BS_to_TS(P, Q, X, Y):
    """Arrow 2. ((P+Q)/2;0_n, (P-Q)/2;0_n, 0_m;(X+Y)/2, 0_m;(X-Y)/2)."""
    m, n = len(P), len(X)
    return ([(P[i] + Q[i]) // 2 for i in range(m)] + [0] * n,
            [(P[i] - Q[i]) // 2 for i in range(m)] + [0] * n,
            [0] * m + [(X[i] + Y[i]) // 2 for i in range(n)],
            [0] * m + [(X[i] - Y[i]) // 2 for i in range(n)])


def check_BS(P, Q, X, Y):
    m, n = len(P), len(X)
    if len(Q) != m or len(Y) != n:
        return False
    if any(v not in (1, -1) for s in (P, Q, X, Y) for v in s):
        return False
    if npaf(P, 0) + npaf(Q, 0) + npaf(X, 0) + npaf(Y, 0) != 2 * (m + n):
        return False
    for i in range(1, m):
        a = npaf(P, i) + npaf(Q, i)
        if i < n:
            a += npaf(X, i) + npaf(Y, i)
        if a:
            return False
    return True


def stabiliser(Tq, t):
    """The full multiplier group {k in Z_t^*: T_i[kq] = T_i[q] for all i,q}."""
    return [k for k in range(1, t)
            if all(Tq[i][(k * q) % t] == Tq[i][q]
                   for i in range(4) for q in range(t))]


def show(x):
    return "".join("+" if v > 0 else ("-" if v < 0 else "0") for v in x)


def unshow(s):
    return [1 if c == "+" else (-1 if c == "-" else 0) for c in s]


def divisors(m):
    return sorted(k for k in range(1, m + 1) if m % k == 0)


def tier_profiles(t, d, parity=True):
    """Every admissible row-sum profile of the d-tier at length t, from the
    orbit-parity row-sum law, recomputed from scratch.

    One class holds 0, so its row sum is eps + d*b (eps = +-1); the other
    three are d*b_i. The sum of squares must be t. With parity=True the
    refinement sum_i b_i = n (mod 2) is also imposed."""
    n = (t - 1) // d
    out = set()
    for eps in (1, -1):
        for b0 in range(-n, n + 1):
            r0 = eps + d * b0
            rem = t - r0 * r0
            if rem < 0 or rem % (d * d):
                continue
            q = rem // (d * d)
            for b1 in range(0, n + 1):
                if b1 * b1 > q:
                    break
                for b2 in range(b1, n + 1):
                    if b1 * b1 + b2 * b2 > q:
                        break
                    for b3 in range(b2, n + 1):
                        if b1 * b1 + b2 * b2 + b3 * b3 != q:
                            continue
                        if parity and (b0 + b1 + b2 + b3 - n) % 2:
                            continue
                        out.add((r0, d * b1, d * b2, d * b3))
    return sorted(out)


# ------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--quick", action="store_true",
                    help="skip every matrix of order > %d"
                         % asm.QUICK_MAX_ORDER)
    args = ap.parse_args()
    t0 = time.time()
    print("=" * 78)
    print("CERT 05 -- T-sequences of order 131 and five Hadamard matrices")
    print("           TT(44) -> BS(87,44) -> TS(131) -> H(4*131*w)")
    print("           trust chain: verify/verify.py")
    print("=" * 78)

    # ------------------------------------------------------------------ [A]
    print("\n[A] the Turyn-type payloads, decoded and checked firsthand")
    pay = json.load(open(PAYLOADS, encoding="ascii"))
    print("    provenance label: %s ..." % pay["provenance"]["label"][:72])
    quads = {}
    for key in sorted(pay["payloads"], key=int):
        n = int(key)
        h = pay["payloads"][key]
        check("TT(%d): payload is %d hex digits = 4n bits, the exact width of "
              "a 4-row array of length n with one pad bit" % (n, len(h)),
              len(h) == n)
        A, B, C, D = decode(h, n)
        r = turyn_residual(A, B, C, D, n)
        check("TT(%d): lengths (%d,%d,%d,%d) = (n,n,n,n-1), alphabet +-1"
              % (n, len(A), len(B), len(C), len(D)),
              [len(A), len(B), len(C), len(D)] == [n, n, n, n - 1]
              and all(v in (1, -1) for s in (A, B, C, D) for v in s))
        check("TT(%d): N_A+N_B+2N_C+2N_D = %d at lag 0 (= 6n-2 = %d) and 0 at "
              "every one of the %d positive lags -- the TURYN-TYPE axiom, own "
              "O(n^2) exact-integer loop" % (n, r[0], 6 * n - 2, n - 1),
              r[0] == 6 * n - 2 and not any(r[1:]))
        bad = decode(h, n, right_align_last=False)
        rb = turyn_residual(*bad, n=n)
        diff = sum(1 for s1, s2 in zip((A, B, C, D), bad)
                   for u, v in zip(s1, s2) if u != v)
        check("TT(%d): the left-aligned reading of the final digit FAILS "
              "(residual cost %d over the positive lags) and differs from the "
              "accepted decode in exactly %d signs, both in the last column"
              % (n, sum(v * v for v in rb[1:]), diff),
              any(rb[1:]) and diff == 2)
        quads[n] = (A, B, C, D)
    check("the banked decode of TT(44) reproduces the decoder's output "
          "character for character",
          [show(x) for x in quads[44]]
          == [pay["decoded_TT44"][k] for k in "ABCD"])

    # ----------------------------------------------------------------- [A']
    print("\n[A'] from-scratch control: the composition is not tuned to 44")
    for n in (6, 8):
        got = brute_TT(n)
        check("TT(%d) FOUND HERE by a complete search of the full 2^(4n-1) = "
              "2^%d space (no external input)" % (n, 4 * n - 1),
              got is not None)
        A, B, C, D = got
        r = turyn_residual(A, B, C, D, n)
        check("  the found TT(%d) satisfies the axiom by the same own loop "
              "(peak %d = 6n-2)" % (n, r[0]),
              r[0] == 6 * n - 2 and not any(r[1:]))
        P, Q, X, Y = TT_to_BS(A, B, C, D)
        check("  TT(%d) -> BS(%d,%d): base-sequence identity holds"
              % (n, 2 * n - 1, n), check_BS(P, Q, X, Y))
        Ts = BS_to_TS(P, Q, X, Y)
        tt = 3 * n - 1
        check("  BS(%d,%d) -> TS(%d): supports partition Z_%d and "
              "sum_i PAF(T_i) = %d*delta_0" % (2 * n - 1, n, tt, tt, tt),
              all(sum(1 for r in Ts if r[q]) == 1 for q in range(tt))
              and sum(paf(r, 0) for r in Ts) == tt
              and not any(sum(paf(r, i) for r in Ts) for i in range(1, tt)))

    # ------------------------------------------------------------------ [B]
    print("\n[B] arrow 1: TT(44) -> BS(87,44)")
    A, B, C, D = quads[44]
    P, Q, X, Y = TT_to_BS(A, B, C, D)
    check("(C;D, C;-D, A, B) has lengths (%d, %d, %d, %d) = (2n-1, 2n-1, n, n)"
          % (len(P), len(Q), len(X), len(Y)),
          [len(P), len(Q), len(X), len(Y)] == [87, 87, 44, 44])
    check("N_P+N_Q+N_X+N_Y = 2*(87+44) = %d at lag 0 and 0 at every positive "
          "lag -- the BASE-SEQUENCE axiom, own exact loop" % (2 * 131),
          check_BS(P, Q, X, Y))
    check("derived: N_P+N_Q = 2(N_C+N_D) identically, so arrow 1 is exactly "
          "the Turyn-type identity re-bracketed",
          all(npaf(P, i) + npaf(Q, i)
              == 2 * (npaf(C, i) + (npaf(D, i) if i < len(D) else 0))
              for i in range(87)))

    # ------------------------------------------------------------------ [C]
    print("\n[C] arrow 2: BS(87,44) -> TS(131), and the T-matrix axiom")
    Tq = BS_to_TS(P, Q, X, Y)
    check("four sequences of length %d over the alphabet {0,+-1}" % len(Tq[0]),
          all(len(r) == T for r in Tq)
          and all(v in (-1, 0, 1) for r in Tq for v in r))
    supp = [sum(1 for v in r if v) for r in Tq]
    rs = [sum(r) for r in Tq]
    check("exactly one T_i nonzero at each of the %d positions: supports %s "
          "partition Z_%d" % (T, supp, T),
          all(sum(1 for r in Tq if r[q]) == 1 for q in range(T))
          and sum(supp) == T and supp == EXPECT_SUPPORTS)
    check("row sums %s, and the derived identity sum_i r_i^2 = %d = t"
          % (rs, sum(v * v for v in rs)),
          rs == EXPECT_ROWSUMS and sum(v * v for v in rs) == T)
    an = asm.aggregate_npaf(Tq)
    check("sum_i NPAF(T_i) = %d*delta_0 -- the NONPERIODIC identity, i.e. "
          "these are T-SEQUENCES, strictly stronger than T-matrices" % an[0],
          an[0] == T and not any(an[1:]))
    ap_ = asm.aggregate_paf(Tq)
    check("sum_i PAF(T_i) = %d*delta_0 -- the T-MATRIX (periodic) axiom, "
          "which is what Cooper-Wallis needs" % ap_[0],
          ap_[0] == T and not any(ap_[1:]))
    K = stabiliser(Tq, T)
    check("full multiplier group = %s (TRIVIAL), so this object lies outside "
          "every K-invariant tier a multiplier search covers" % K, K == [1])
    check("the banked literature-route TS(131) reproduces this construction "
          "character for character",
          [show(r) for r in Tq]
          == json.load(open(T131_LIT, encoding="ascii"))["T"])

    # ----------------------------------------------------------------- [C']
    print("\n[C'] the SECOND, INDEPENDENT witness: a K-invariant T-matrix "
          "quadruple of order\n     131 from this laboratory's multiplier-tier "
          "search, and the tier lattice")
    print("     d    n    profiles (sum r^2 = 131 only) | + orbit-parity")
    alive = []
    for d in divisors(T - 1):
        n = (T - 1) // d
        raw = tier_profiles(T, d, parity=False)
        par = tier_profiles(T, d, parity=True)
        print("     %-4d %-4d %-28d | %d" % (d, n, len(raw), len(par)))
        if par:
            alive.append(d)
    dead = [d for d in divisors(T - 1) if not tier_profiles(T, d, parity=False)]
    check("the alive multiplier tiers at t=131 are exactly d in %s; every "
          "other divisor of 130 -- %s -- admits no profile even solving "
          "sum r_i^2 = 131 (d=2 dies on 131 - r0^2 = 2 mod 4 for odd r0)"
          % (alive, dead),
          alive == [1, 5] and dead == [2, 10, 13, 26, 65, 130])
    p5 = tier_profiles(T, 5, parity=True)
    norm5 = set(tuple(sorted(abs(x) for x in r)) for r in p5)
    check("the d=5 tier carries EXACTLY ONE admissible profile up to the "
          "per-class sign flip: |r| = %s (the two signed forms %s are the same "
          "profile)" % (sorted(norm5)[0], p5), norm5 == {(0, 5, 5, 9)})
    mul = json.load(open(T131_MUL, encoding="ascii"))
    Tm = [unshow(r) for r in mul["solutions"][0]["T"]]
    suppm = [sum(1 for v in r if v) for r in Tm]
    rsm = [sum(r) for r in Tm]
    check("the banked witness: supports %s partition Z_131, row sums %s"
          % (suppm, rsm),
          all(sum(1 for r in Tm if r[q]) == 1 for q in range(T))
          and suppm == [11, 10, 55, 55] and sorted(rsm) == [-9, 0, 5, 5])
    am = asm.aggregate_paf(Tm)
    check("sum_i PAF(T_i) = %d*delta_0 -- the T-MATRIX axiom, own exact loop"
          % am[0], am[0] == T and not any(am[1:]))
    Km = stabiliser(Tm, T)
    check("its full multiplier group has order %d, K = %s, and every k in K "
          "satisfies k^5 = 1 mod 131 -- the unique order-5 subgroup of Z_131^*"
          % (len(Km), Km),
          len(Km) == 5 and all(pow(k, 5, T) == 1 for k in Km))
    nm = asm.aggregate_npaf(Tm)
    nz = sum(1 for i in range(1, T) if nm[i])
    print("      [--] recorded, NOT a failure: this witness's aggregate "
          "NONPERIODIC autocorrelation\n           is nonzero at %d of the %d "
          "positive lags. It is a T-MATRIX quadruple,\n           exactly "
          "Cooper-Wallis's hypothesis, and not an aperiodic T-sequence,\n"
          "           unlike [C]'s." % (nz, T - 1))
    check("THE TWO WITNESSES ARE DIFFERENT OBJECTS: support multisets %s vs "
          "%s, row-sum multisets %s vs %s, multiplier groups of order %d vs "
          "%d, aggregate NPAF %s vs nonzero off-peak -- no equivalence can "
          "identify them"
          % (sorted(EXPECT_SUPPORTS), sorted(suppm), sorted(EXPECT_ROWSUMS),
             sorted(rsm), len(K), len(Km), "131*delta_0"),
          sorted(EXPECT_SUPPORTS) != sorted(suppm)
          and sorted(EXPECT_ROWSUMS) != sorted(rsm)
          and len(K) != len(Km)
          and not any(an[1:]) and any(nm[1:]))

    # ------------------------------------------------------------------ [D]
    print("\n[D] the Williamson quadruples the cash-in assumes")
    bank = asm.load_bank("williamson-9-13.json")
    for w in WS:
        W = asm.load_W(bank, w)
        a = asm.aggregate_paf(W)
        wr = [sum(r) for r in W]
        check("w=%-2d: four +-1 rows, symmetric circulant, sum_q PAF_q = "
              "%d*delta_0, sum (row sum)^2 = %d"
              % (w, a[0], sum(v * v for v in wr)),
              all(v in (1, -1) for r in W for v in r)
              and all(r[j] == r[(-j) % w] for r in W for j in range(w))
              and a[0] == 4 * w and not any(a[1:])
              and sum(v * v for v in wr) == 4 * w)

    # ------------------------------------------------------------------ [E]
    print("\n[E] the five target orders, and what the table says about them")
    print("    %-8s %-8s %s" % ("order", "n = 131w", "Cati-Pasechnik Table 4 "
                                "v2 entry"))
    for order in sorted(CP_TABLE4):
        n, m = CP_TABLE4[order]
        print("    %-8d %-8d %d(%d) -- their smallest known was 2^%d*%d = %d"
              % (order, n, n, m, m, n, (1 << m) * n))
    print("    Label: REPORTED-FROM-AUDITED-TABLE. The five entries are quoted "
          "from Table 4 of\n    arXiv:2411.18897v2 (2025-08-30), machine-"
          "diffed; re-check at release date.")
    cp, _ = asm.load_cp_table4()
    quoted = dict(CP_TABLE4.values())
    check("the five entries above agree with data/cp-table4-v2.json (%d "
          "entries), and t=131 is unresolved there at exactly w in %s: these "
          "five targets, exactly"
          % (len(cp), asm.cp_open_w(cp, 131)),
          asm.cp_agrees(cp, quoted)
          and asm.cp_open_w(cp, 131) == [7, 11, 13, 17, 19])

    # ------------------------------------------------------------------ [F]
    print("\n[F] calibration gates (known orders)")
    gates = sorted(GATE_PINS)
    if args.quick:
        skipped = [GATE_PINS[k][0] for k in gates
                   if GATE_PINS[k][0] > asm.QUICK_MAX_ORDER]
        gates = [k for k in gates if GATE_PINS[k][0] <= asm.QUICK_MAX_ORDER]
        print("    --quick: skipping orders > %d: %s"
              % (asm.QUICK_MAX_ORDER, skipped))
    for (tt, w) in gates:
        order, pin = GATE_PINS[(tt, w)]
        tag = "gate t=%d x w=%d -> H(%d)" % (tt, w, order)
        print("\n    %s" % tag)
        asm.build_and_pin(Tq, asm.load_W(bank, w), tt, w, order, pin, tag)

    # ------------------------------------------------------------------ [G]
    orders = sorted(TARGET_PINS)
    if args.quick:
        skipped = [o for o in orders if o > asm.QUICK_MAX_ORDER]
        orders = [o for o in orders if o <= asm.QUICK_MAX_ORDER]
        print("\n[G] the targets, witness 1  (--quick: skipping orders > %d: "
              "%s)" % (asm.QUICK_MAX_ORDER, skipped))
    else:
        print("\n[G] the targets, witness 1 (all five)")
    verdicts = []
    for order in orders:
        w, pin = TARGET_PINS[order]
        ok, sha = asm.build_and_pin(Tq, asm.load_W(bank, w), T, w, order, pin,
                                    "t=131 x w=%d -> H(%d)" % (w, order),
                                    verbose=False)
        verdicts.append((order, w, sha, ok))

    # ----------------------------------------------------------------- [G']
    print("\n[G'] the same orders, rebuilt from the second witness")
    w2 = sorted(WITNESS2_PINS)
    if args.quick:
        skipped = [o for o in w2 if o > asm.QUICK_MAX_ORDER]
        w2 = [o for o in w2 if o <= asm.QUICK_MAX_ORDER]
        print("     --quick: skipping orders > %d: %s"
              % (asm.QUICK_MAX_ORDER, skipped))
    for order in w2:
        w, pin = WITNESS2_PINS[order]
        ok, sha = asm.build_and_pin(Tm, asm.load_W(bank, w), T, w, order, pin,
                                    "witness 2, t=131 x w=%d -> H(%d)"
                                    % (w, order), verbose=False)
        other = dict((o, p) for o, (_, p) in TARGET_PINS.items())
        other[524] = GATE_PINS[(131, 1)][1]
        check("H(%d): the two witnesses give DIFFERENT matrices (%s... vs "
              "%s...)" % (order, sha[:12], other[order][:12]),
              sha != other[order])

    # ------------------------------------------------------------------ [H]
    print("\n[H] negative control")
    asm.tamper_reject(Tq, asm.load_W(bank, 1), T, 1, 524)

    # ------------------------------------------------------------------ [I]
    print("\n[I] verdicts")
    print("VERDICT TS(131): four (0,+-1) sequences, supports %s partitioning "
          "Z_131,\n        aggregate NPAF and PAF both = 131*delta_0  |  "
          "PROVEN-BY-CERTIFICATE" % EXPECT_SUPPORTS)
    print("VERDICT the second witness is a T-matrix quadruple of order 131 "
          "with multiplier\n        group of order 5, inequivalent to the "
          "first  |  PROVEN-BY-CERTIFICATE")
    for order, w, sha, ok in verdicts:
        n, m = CP_TABLE4[order]
        print("VERDICT H(%-5d) = 4*131*%-2d  Table 4 v2 entry %d(%d)  |  %s"
              "\n        sha=%s"
              % (order, w, n, m,
                 "PROVEN-BY-CERTIFICATE" if ok else "FAILED", sha))
    if not verdicts:
        print("        --quick: no target was built, so no verdict line on "
              "the five orders is\n        established by this run. Run "
              "without --quick for them.")
    print("\n        No novelty of existence is claimed for these five "
          "orders. A published\n        route closes them: London and "
          "Kotsireas 2025 supply TT(44) and both\n        compositions above "
          "are classical, as is Cooper-Wallis. What is claimed is\n        the "
          "first publicly accessible and independently reproducible "
          "construction\n        of H(3668), H(5764), H(6812), H(8908) and "
          "H(9956) located in our audit of\n        2026-08-29.")

    print("\n" + "=" * 78)
    print("-- elapsed %.1fs --" % (time.time() - t0))
    if asm.FAIL:
        print("CERT 05: FAIL (%d): %s" % (len(asm.FAIL), asm.FAIL))
        return 1
    print("CERT 05 GREEN%s" % (" (--quick)" if args.quick else ""))
    print("=" * 78)
    return 0


if __name__ == "__main__":
    sys.exit(main())
