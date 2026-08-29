#!/usr/bin/env python3
"""Cert 03 -- a T-matrix quadruple of order 163, and H(7172) + H(8476).

See NOTES.md. Exit 0 iff every check passes.

The order-163 T-matrix quadruple found by this laboratory is pushed through
the Cooper-Wallis product

    T-matrices of order t  x  Williamson-type of order w   ==>   H(4tw)

for w = 11 and w = 13, yielding

    H(7172) = 4*163*11        H(8476) = 4*163*13

  [A] the banked order-163 quadruple: alphabet {0,+-1}, supports partitioning
      Z_163, sum_i PAF(T_i) = 163*delta_0 by an own exact-integer double loop,
      sum_i r_i^2 = 163, and strict invariance under the unique order-9
      subgroup K of Z_163^* (re-derived here by closure, not read in);
  [B] the Williamson quadruples w = 1, 11, 13: +-1, symmetric circulant,
      sum_q PAF_q = 4w*delta_0 by the same own loop;
  [C] calibration gate t=163 x w=1 -> H(652), digest pinned three ways, plus a
      negative control (one flipped character must make verify.py exit 1);
  [D] the two targets, built one at a time (build -> write -> verify ->
      delete), each handed to verify/verify.py in a subprocess.

Standard library only. The trust chain is verify/verify.py alone; the
assembler is tools/assemble.py.

    python certs/03-t163-products/run.py           # ~1 min
    python certs/03-t163-products/run.py --quick   # gate only, no order > 1200
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
T = 163
D = 9
check = asm.check

# Canonical SHA-256 of verify.py's '+'/'-' serialisation.
PINS = {
    652: (1,
          "f488e58cad98e11624fba4c113fbd34b3a77bb1099f96f66772918db066b1c1a"),
    7172: (11,
           "8535de5debeae63611d22d39fc3a4b821404dd4320699202b8f17c6851e46aa2"),
    8476: (13,
           "25286dd5c3b1500ad6db0ebfdfb7aa26067dd3a4c779137dbdc18dab8a8dd77c"),
}
# Cati-Pasechnik arXiv:2411.18897v2 Table 4, entries n(m).
CP_TABLE4 = {7172: (1793, 4), 8476: (2119, 4)}


def order9_subgroup(p):
    """Re-derive the unique order-9 subgroup of Z_p^* by closure."""
    for g in range(2, p):
        K, x = [1], g
        while x != 1:
            K.append(x)
            x = x * g % p
        if len(K) == 9:
            K = sorted(K)
            if all(a * b % p in K for a in K for b in K):
                return K
    raise SystemExit("no order-9 subgroup")


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--quick", action="store_true",
                    help="skip every matrix of order > %d"
                         % asm.QUICK_MAX_ORDER)
    args = ap.parse_args()
    t0 = time.time()
    print("=" * 78)
    print("CERT 03 -- T(163) and the closures H(7172), H(8476)")
    print("           trust chain: verify/verify.py")
    print("=" * 78)

    # ------------------------------------------------------------------ [A]
    print("\n[A] the banked order-163 quadruple (data/T163.json)")
    d = json.load(open(os.path.join(DATA, "T163.json"), encoding="ascii"))
    if d["n"] != T:
        raise SystemExit("data/T163.json: n mismatch")
    cls, sgn = d["classes"], d["signs"]
    Tq = [[0] * T for _ in range(4)]
    for q in range(T):
        Tq[cls[q]][q] = sgn[q]
    check("alphabet {0,+-1}", all(x in (-1, 0, 1) for r in Tq for x in r))
    check("exactly one T_i nonzero at each of the 163 positions",
          all(sum(1 for r in Tq if r[q]) == 1 for q in range(T)))
    supp = [sum(1 for x in r if x) for r in Tq]
    rs = [sum(r) for r in Tq]
    print("      supports = %s (total %d) ; row sums = %s ; sum r^2 = %d"
          % (supp, sum(supp), rs, sum(x * x for x in rs)))
    check("supports match the banked profile %s" % d["support_sizes"],
          supp == list(d["support_sizes"]) and sum(supp) == T)
    check("sum_i r_i^2 = 163", sum(x * x for x in rs) == T)
    agg = asm.aggregate_paf(Tq)
    check("sum_i PAF(T_i) = 163*delta_0 -- the T-MATRIX (PERIODIC) axiom, "
          "which is exactly the Cooper-Wallis hypothesis",
          agg[0] == T and not any(agg[1:]))
    K = order9_subgroup(T)
    check("K re-derived by closure = the banked K (order 9): %s" % K,
          K == sorted(d["multiplier"]["K"]) and len(K) == D)
    check("T is strictly invariant under K",
          all(Tq[i][(k * q) % T] == Tq[i][q]
              for i in range(4) for k in K for q in range(T)))
    npaf = asm.aggregate_npaf(Tq)
    nz = sum(1 for s in range(1, T) if npaf[s])
    print("      [--] recorded, NOT a failure: the aggregate NON-PERIODIC "
          "autocorrelation is\n           nonzero at %d of the %d positive "
          "lags. This is a T-MATRIX quadruple,\n           not an aperiodic "
          "T-sequence." % (nz, T - 1))
    check("the NPAF peak still reads t (%d)" % npaf[0], npaf[0] == T)

    # ------------------------------------------------------------------ [B]
    print("\n[B] the Williamson inputs")
    bank = asm.load_bank("williamson-13.json")
    for w in (1, 11, 13):
        W = asm.load_W(bank, w)
        a = asm.aggregate_paf(W)
        check("w=%2d: +-1, symmetric circulant, sum_q PAF_q = %d*delta_0, "
              "sum (row sum)^2 = %d" % (w, 4 * w, 4 * w),
              all(x in (1, -1) for r in W for x in r)
              and all(r[k] == r[(-k) % w] for r in W for k in range(w))
              and a[0] == 4 * w and not any(a[1:])
              and sum(sum(r) ** 2 for r in W) == 4 * w)

    # ------------------------------ [B'] the table reading, checked
    print("\n[B'] the openness reading, checked against the transcription")
    cp, _ = asm.load_cp_table4()
    quoted = dict(CP_TABLE4.values())
    check("the %d entries quoted in the verdicts -- %s -- agree with "
          "data/cp-table4-v2.json (%d entries), and t=163 is unresolved there "
          "at exactly w in %s: these two targets, exactly"
          % (len(quoted),
             " ".join("%d(%d)" % (n, m) for n, m in sorted(quoted.items())),
             len(cp), asm.cp_open_w(cp, 163)),
          asm.cp_agrees(cp, quoted) and asm.cp_open_w(cp, 163) == [11, 13])

    # ---------------------------------------------------------------- [C][D]
    orders = sorted(PINS)
    if args.quick:
        skipped = [o for o in orders if o > asm.QUICK_MAX_ORDER]
        orders = [o for o in orders if o <= asm.QUICK_MAX_ORDER]
        print("\n    --quick: skipping orders > %d: %s"
              % (asm.QUICK_MAX_ORDER, skipped))
    results = []
    for order in orders:
        w, pin = PINS[order]
        tag = ("gate H(652)" if order == 652 else "H(%d)" % order)
        print("\n[%s] %s = 4*163*%d"
              % ("C" if order == 652 else "D", tag, w))
        ok, sha = asm.build_and_pin(Tq, asm.load_W(bank, w), T, w, order, pin,
                                    tag)
        results.append((order, w, sha, ok))

    # -------------------------------------------------- negative control
    print("\n[C'] negative control: one flipped character must be rejected")
    asm.tamper_reject(Tq, asm.load_W(bank, 1), T, 1, 652)

    # ------------------------------------------------------------ verdicts
    print("\n[E] verdicts")
    for order, w, sha, ok in results:
        label = "GATE (known order)" if order == 652 else "PROVEN-BY-CERTIFICATE"
        print("VERDICT H(%-5d) = 4*163*%-2d  |  %s  sha=%s"
              % (order, w, label if ok else "FAILED", sha))
    got = [o for o, _, _, _ in results if o in CP_TABLE4]
    if got:
        print("        Openness of %s, hedged and dated: listed at m = 4 in "
              "Cati-Pasechnik\n        Table 4 v2 (2025-08-30), machine-"
              "diffed, entries %s; re-check at release\n        date. m = 4 "
              "means neither 2^2*n nor 2^3*n carries an entry there, so the "
              "same\n        two matrices also give H(14344) and H(16952) by "
              "Sylvester doubling."
              % (", ".join("H(%d)" % o for o in got),
                 ", ".join("%d(%d)" % CP_TABLE4[o] for o in got)))

    print("\n" + "=" * 78)
    print("-- elapsed %.1fs --" % (time.time() - t0))
    if asm.FAIL:
        print("CERT 03: FAIL (%d): %s" % (len(asm.FAIL), asm.FAIL))
        return 1
    print("CERT 03 GREEN%s" % (" (--quick)" if args.quick else ""))
    print("=" * 78)
    return 0


if __name__ == "__main__":
    sys.exit(main())
