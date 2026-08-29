#!/usr/bin/env python3
"""Cert 02 -- six Hadamard matrices from the order-103 T-matrix quadruple.

See NOTES.md. Exit 0 iff every check passes.

The order-103 T-matrix quadruple certified by cert 01 and symmetric circulant
Williamson quadruples of order w are pushed through the Cooper-Wallis product

    T-matrices of order t  x  Williamson-type of order w   ==>   H(4tw)

for w in {7, 11, 17, 19, 23, 29}, yielding

    H(2884), H(4532), H(7004), H(7828), H(9476), H(11948)

  [A] the banked Williamson quadruples: +-1, symmetric circulant, and
      sum_q PAF_q = 4w*delta_0 by an own exact-integer double loop;
  [B] the banked order-103 quadruple: alphabet {0,+-1}, disjoint supports of
      sizes (33,30,19,21) partitioning Z_103, sum_i PAF(T_i) = 103*delta_0 by
      the same own loop. Also recorded: the aggregate NON-PERIODIC
      autocorrelation does not vanish, so the object is a T-MATRIX quadruple,
      which is exactly Cooper-Wallis's hypothesis, and not an aperiodic
      T-sequence;
  [C] calibration gates: t=3 (the brute-forced control of cert 01) x w=3 ->
      H(36), and t=103 x w=1 -> H(412), both handed to verify/verify.py with
      their canonical digests pinned; plus two negative controls (a perturbed
      Williamson row, and a one-character tamper of a green matrix);
  [D] the six targets, built one at a time, written, verified by
      verify/verify.py in a subprocess, digest pinned three ways, deleted;
  [E] one verdict line per order.

Standard library only. Exact integer arithmetic only. The trust chain is
verify/verify.py alone; the assembler is tools/assemble.py.

    python certs/02-t103-products/run.py           # all six (~3 min)
    python certs/02-t103-products/run.py --quick   # gates only, no order > 1200
"""

import argparse
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.join(ROOT, "tools"))
import assemble as asm  # noqa: E402

DATA = os.path.join(ROOT, "data")
T103_JSON = os.path.join(DATA, "T103.json")
T3_JSON = os.path.join(DATA, "T3-control.json")

T = 103
EXPECT_SUPPORTS = [33, 30, 19, 21]
check = asm.check

# Canonical SHA-256 of verify.py's '+'/'-' serialisation.
GATE_PINS = {
    36: ("t3", 3,
         "25f0a20be9d4ae4f29166241f2c324b772303740ff5c5538975bc17f1ec17fa9"),
    412: ("t103", 1,
          "fd7d2e174d3b52fecdc20cce963a3c278b00d07bc392018e0b8ded841991512e"),
}
TARGET_PINS = {
    2884: (7,
           "38c778c1de4d5a9af4e7a9391382ee8b98ca89143a4885efff0d013429ebff8e"),
    4532: (11,
           "3f00bb1954cc7fa77e2112ea147a76b59de18142fe17b28348d3c266f21a7348"),
    7004: (17,
           "50b9c71d6683a3bbf7bc1224243c8903932a1919d9683a91de49377433b707f0"),
    7828: (19,
           "6ea2946fe2ba73e5e06cf77fda33e99004a885228844e5da6deea96cb05286c5"),
    9476: (23,
           "ea4c04d1759029421df7572450365da4a16c466b9d435dfcc79930a43e548f96"),
    11948: (29,
            "8169db47b795b30cfc6036f4b090380902dd06f9e72726f0f9e2fe8120ed550b"),
}
# Cati-Pasechnik arXiv:2411.18897v2 Table 4, entries n(m): m is minimal such
# that H(2^m * n) is known. m = 3 means 2^2*n carries no entry there.
CP_TABLE4 = {2884: (721, 3), 4532: (1133, 3), 7004: (1751, 3),
             7828: (1957, 3), 9476: (2369, 3), 11948: (2987, 3)}


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--quick", action="store_true",
                    help="skip every matrix of order > %d"
                         % asm.QUICK_MAX_ORDER)
    args = ap.parse_args()
    t0 = time.time()
    print("=" * 78)
    print("CERT 02 -- Cooper-Wallis closures from the order-103 T-matrix "
          "quadruple")
    print("           trust chain: verify/verify.py")
    print("=" * 78)

    # ------------------------------------------- [A] banked Williamson bank
    print("\n[A] banked Williamson quadruples (data/williamson-quadruples.json)")
    bank = asm.load_bank()
    ws = sorted(int(k) for k in bank["quadruples"])
    check("bank holds exactly the gate + target orders "
          "{1, 3} u {7, 11, 17, 19, 23, 29}",
          ws == [1, 3, 7, 11, 17, 19, 23, 29])
    for w in ws:
        W = asm.load_W(bank, w)
        ent = bank["quadruples"][str(w)]
        pm1 = all(v in (1, -1) for r in W for v in r)
        widths = all(len(r) == w for r in W) and len(W) == 4
        sym = all(r[k] == r[(-k) % w] for r in W for k in range(w))
        agg = asm.aggregate_paf(W)
        rs = [sum(r) for r in W]
        ok = (pm1 and widths and sym and agg[0] == 4 * w
              and all(v == 0 for v in agg[1:]) and rs == ent["row_sums"]
              and sum(v * v for v in rs) == 4 * w)
        check("w=%-2d  4 rows of +-1, symmetric circulant, "
              "sum_q PAF_q = %d*delta_0, row sums %s with sum of squares %d"
              % (w, agg[0], rs, sum(v * v for v in rs)), ok)

    # ------------------------------------------ [B] the order-103 quadruple
    print("\n[B] the T-matrix quadruple of order 103 (data/T103.json)")
    cls, sgn, Tq, t, d103 = asm.load_T(T103_JSON)
    check("data/T103.json loads at n = 103", t == T)
    supports = [sum(1 for q in range(t) if cls[q] == i) for i in range(4)]
    check("alphabet {0, +-1} and exactly one T_i nonzero at each of the %d "
          "positions (supports %s partition Z_103)" % (t, supports),
          all(v in (0, 1, -1) for r in Tq for v in r)
          and all(sum(1 for i in range(4) if Tq[i][q]) == 1 for q in range(t))
          and sum(supports) == t)
    check("support sizes %s match the banked field %s"
          % (supports, d103["support_sizes"]),
          supports == EXPECT_SUPPORTS == d103["support_sizes"])
    agg = asm.aggregate_paf(Tq)
    check("sum_i PAF(T_i) = %d*delta_0 -- the T-MATRIX (periodic) axiom, "
          "own O(4n^2) loop" % agg[0],
          agg[0] == t and all(v == 0 for v in agg[1:]))
    rs = [sum(r) for r in Tq]
    check("derived identity sum_i (row sum)^2 = sum_s aggPAF(s) = %d = t"
          % sum(v * v for v in rs), sum(v * v for v in rs) == t)
    nagg = asm.aggregate_npaf(Tq)
    nz = sum(1 for s in range(1, t) if nagg[s])
    print("      [--] recorded, NOT a failure: the aggregate NON-PERIODIC "
          "autocorrelation is\n           nonzero at %d of the %d positive "
          "lags, so this quadruple is a T-MATRIX\n           quadruple, which "
          "is exactly Cooper-Wallis's hypothesis, and not an\n           "
          "aperiodic T-sequence." % (nz, t - 1))
    check("the NPAF peak still reads n (%d)" % nagg[0], nagg[0] == t)

    cls3, sgn3, T3, t3, _ = asm.load_T(T3_JSON)
    agg3 = asm.aggregate_paf(T3)
    check("the brute-forced length-3 control satisfies the same axiom "
          "(sum_i PAF = %d*delta_0)" % agg3[0],
          t3 == 3 and agg3[0] == 3 and all(v == 0 for v in agg3[1:]))

    # ------------------------------------------------- [C] calibration gates
    print("\n[C] calibration gates")
    for order, (tname, w, pin) in sorted(GATE_PINS.items()):
        Tsrc = T3 if tname == "t3" else Tq
        tt = 3 if tname == "t3" else T
        asm.build_and_pin(Tsrc, asm.load_W(bank, w), tt, w, order, pin,
                          "gate %s x w=%d -> H(%d)" % (tname, w, order))

    print("\n      negative controls")
    Wbad = asm.load_W(bank, 3)
    Wbad[1] = [-Wbad[1][0]] + Wbad[1][1:]        # break the w=3 quadruple
    rcb, vb, _, _ = asm.build(T3, Wbad, 3, 3, "perturbed w=3", verbose=False)
    check("perturbing one Williamson entry aborts before assembly "
          "(pre-check catches it: %r)" % vb, rcb != 0 and "PRE-CHECK" in vb)
    asm.tamper_reject(T3, asm.load_W(bank, 3), 3, 3, 36)

    # ---------------------------------------- [C'] the table reading, checked
    print("\n[C'] the openness reading, checked against the transcription")
    cp, _ = asm.load_cp_table4()
    quoted = dict(CP_TABLE4.values())
    check("the %d entries quoted in the verdicts -- %s -- agree with "
          "data/cp-table4-v2.json (%d entries), and t=103 is unresolved there "
          "at exactly w in %s: these six targets plus w=5, that is n=515 and "
          "H(2060), closed before this repository was cut"
          % (len(quoted),
             " ".join("%d(%d)" % (n, m) for n, m in sorted(quoted.items())),
             len(cp), asm.cp_open_w(cp, 103)),
          asm.cp_agrees(cp, quoted)
          and asm.cp_open_w(cp, 103) == [5, 7, 11, 17, 19, 23, 29])

    # ------------------------------------------------------- [D] the targets
    orders = sorted(TARGET_PINS)
    if args.quick:
        skipped = [o for o in orders if o > asm.QUICK_MAX_ORDER]
        orders = [o for o in orders if o <= asm.QUICK_MAX_ORDER]
        print("\n[D] targets  (--quick: skipping orders > %d: %s)"
              % (asm.QUICK_MAX_ORDER, skipped))
    else:
        print("\n[D] targets  (all six)")
    verdicts = []
    for order in orders:
        w, pin = TARGET_PINS[order]
        ok, sha = asm.build_and_pin(Tq, asm.load_W(bank, w), T, w, order, pin,
                                    "t=103 x w=%d -> H(%d)" % (w, order))
        n_cp, m_cp = CP_TABLE4[order]
        verdicts.append((order, w, n_cp, m_cp, sha, ok))

    # -------------------------------------------------------- [E] verdicts
    print("\n[E] verdicts")
    for order, w, n_cp, m_cp, sha, ok in verdicts:
        print("VERDICT H(%-5d) = 4*103*%-2d  |  %s  sha=%s"
              % (order, w, "PROVEN-BY-CERTIFICATE" if ok else "FAILED", sha))
    if verdicts:
        print("        Openness of these orders, hedged and dated: each is "
              "listed at m = 3 in\n        Cati-Pasechnik Table 4 v2 "
              "(2025-08-30), machine-diffed, entries %s;\n        re-check at "
              "release date. A database records what its authors knew, not\n"
              "        what exists."
              % ", ".join("%d(%d)" % (n, m) for _, _, n, m, _, _ in verdicts))
    else:
        print("        no target built under --quick; the gates above are "
              "known orders.")

    print("\n" + "=" * 78)
    print("-- elapsed %.1fs --" % (time.time() - t0))
    if asm.FAIL:
        print("CERT 02: FAIL (%d): %s" % (len(asm.FAIL), asm.FAIL))
        return 1
    print("CERT 02 GREEN%s" % (" (--quick)" if args.quick else ""))
    print("=" * 78)
    return 0


if __name__ == "__main__":
    sys.exit(main())
