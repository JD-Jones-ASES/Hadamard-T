#!/usr/bin/env python3
"""Cert 04 -- T-matrix quadruples of orders 101, 113 and 181.

See NOTES.md. Exit 0 iff every check passes.

Three K-invariant T-matrix quadruples found by this laboratory are re-verified
from scratch and banked. They close no order: for every t in {101, 113, 181}
and every Williamson order w this repository holds, the product order 4tw is
already known, so this certificate builds calibration gates only.

  [A] the three banked quadruples: alphabet {0,+-1}, disjoint supports
      partitioning Z_t, sum_i PAF(T_i) = t*delta_0 by an own exact-integer
      O(4t^2) double loop, sum_i r_i^2 = t, and the FULL multiplier group
      re-derived by direct stabiliser computation (not read in) and shown to
      be the unique subgroup of Z_t^* of its order;
  [B] the ten Williamson quadruples: +-1, symmetric circulant,
      sum_q PAF_q = 4w*delta_0, by the same own loop;
  [C] what the three quadruples buy: the Cati-Pasechnik Table 4 transcription
      in data/cp-table4-v2.json is PARSED at run time, two controls bind that
      parse to the two t this repository closes orders at, and the sweep of
      all three t against every odd w <= 63 is re-derived from it;
  [D] calibration gates -- t x w=1 and t x w=3 for each t, six known orders,
      each built, handed to verify/verify.py in a subprocess, digest pinned
      three ways, deleted. w=1 is degenerate in the group product, so the w=3
      gates are the ones that exercise Z_t x Z_w;
  [D'] a negative control: one flipped character must make verify.py exit 1.

Standard library only. The trust chain is verify/verify.py alone; the
assembler is tools/assemble.py.

    python certs/04-t-bank/run.py           # ~5 s
    python certs/04-t-bank/run.py --quick   # skips the three gates over 1200
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
TS = (101, 113, 181)
WS = (1, 3, 7, 9, 11, 13, 17, 19, 23, 29)
check = asm.check

# The two binding controls on the parse. Re-reading the same file at the two
# t this repository does close orders at must reproduce their unresolved sets
# exactly; a transcription that could not would not be trustworthy for the
# negative result in [C].
CONTROL_W = {103: [5, 7, 11, 17, 19, 23, 29], 163: [11, 13]}
CONTROL_NOTE = {
    103: "cert 02's six targets, plus w=5, that is n=515 and H(2060), closed "
         "before this repository was cut and so never a cert-02 target",
    163: "cert 03's two targets, exactly",
}

PROFILE = {
    101: {"K_order": 5, "supports": [21, 20, 20, 40],
          "row_sums": [-1, 0, -10, 0]},
    113: {"K_order": 7, "supports": [36, 35, 28, 14],
          "row_sums": [-8, 7, 0, 0]},
    181: {"K_order": 9, "supports": [82, 36, 36, 27],
          "row_sums": [10, 0, 0, -9]},
}
# Canonical SHA-256 of verify.py's '+'/'-' serialisation. (t, w) -> (order, sha)
GATE_PINS = {
    (101, 1): (404,
               "b9a1d492851c3b9e9a81c672647c89c6e4551f43fe9dc0a7ce9eecc5aad9f0ed"),
    (101, 3): (1212,
               "36ffec1af9be598bde2ec6682cfead01e833b3a69d970c7cff4ebd8b093637eb"),
    (113, 1): (452,
               "1ab70fee6221a39bf2fbf37e1cd10c7913153cbf00de8ee57e2a880c08ea167e"),
    (113, 3): (1356,
               "2aa915a52363e8c7bb8f21d819ff6dbc97f715e55e0226e6cf067cf66b019474"),
    (181, 1): (724,
               "a4d006afa086bf2a88e14f2bc13678b8cb85fefc058d74fbcd1581a27df96df2"),
    (181, 3): (2172,
               "7f4660c332ccaa25a805058b801c4291a0b33ae8e58655b623b7a81174d5a4e1"),
}


def stabiliser(Tq, t):
    """The FULL multiplier group {k in Z_t^*: T_i[kq] = T_i[q] for all i,q}."""
    return [k for k in range(1, t)
            if all(Tq[i][(k * q) % t] == Tq[i][q]
                   for i in range(4) for q in range(t))]


def unique_subgroup_of_order(p, d):
    """The unique order-d subgroup of the cyclic group Z_p^*, by closure."""
    for g in range(2, p):
        K, x = [1], g
        while x != 1:
            K.append(x)
            x = x * g % p
        if len(K) == d:
            K = sorted(K)
            if all(a * b % p in K for a in K for b in K):
                return K
    return None


def load_T(t):
    d = json.load(open(os.path.join(DATA, "T%d.json" % t), encoding="ascii"))
    if d["n"] != t:
        raise SystemExit("data/T%d.json: n mismatch" % t)
    Tq = [[0] * t for _ in range(4)]
    for q in range(t):
        Tq[d["classes"][q]][q] = d["signs"][q]
    return d, Tq


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--quick", action="store_true",
                    help="skip every matrix of order > %d"
                         % asm.QUICK_MAX_ORDER)
    args = ap.parse_args()
    t0 = time.time()
    print("=" * 78)
    print("CERT 04 -- T-matrix quadruples of orders 101, 113 and 181")
    print("           trust chain: verify/verify.py")
    print("=" * 78)

    # ------------------------------------------------------------------- [A]
    print("\n[A] the three banked T-matrix quadruples")
    quads = {}
    for t in TS:
        d, Tq = load_T(t)
        print("\n    --- t = %d ---" % t)
        check("t=%d: alphabet {0,+-1}" % t,
              all(v in (-1, 0, 1) for r in Tq for v in r))
        check("t=%d: exactly one T_i nonzero at each of the %d positions "
              "(supports are disjoint and cover Z_%d)" % (t, t, t),
              all(sum(1 for i in range(4) if Tq[i][q]) == 1
                  for q in range(t)))
        supp = [sum(1 for v in r if v) for r in Tq]
        rs = [sum(r) for r in Tq]
        check("t=%d: support sizes %s (total %d) match the banked field %s"
              % (t, supp, sum(supp), d["support_sizes"]),
              supp == list(d["support_sizes"]) and sum(supp) == t
              and supp == PROFILE[t]["supports"])
        check("t=%d: row sums %s match the banked field %s"
              % (t, rs, d["row_sums"]),
              rs == list(d["row_sums"]) == PROFILE[t]["row_sums"])
        check("t=%d: derived identity sum_i r_i^2 = %d = t"
              % (t, sum(v * v for v in rs)), sum(v * v for v in rs) == t)
        agg = asm.aggregate_paf(Tq)
        check("t=%d: sum_i PAF(T_i) = %d*delta_0 -- the T-MATRIX (PERIODIC) "
              "axiom, own O(4t^2) exact-integer loop" % (t, agg[0]),
              agg[0] == t and not any(agg[1:]))
        K = stabiliser(Tq, t)
        U = unique_subgroup_of_order(t, len(K))
        check("t=%d: full multiplier group has order %d, K = %s"
              % (t, len(K), K),
              len(K) == PROFILE[t]["K_order"] == d["multiplier"]["order"]
              and K == sorted(d["multiplier"]["K"]))
        check("t=%d: K is exactly the unique order-%d subgroup of the cyclic "
              "group Z_%d^* (re-derived by closure)" % (t, len(K), t),
              U is not None and U == K)
        check("t=%d: T is strictly K-invariant, and the %d K-orbits are the "
              "design's free cells" % (t, 1 + (t - 1) // len(K)),
              all(Tq[i][(k * q) % t] == Tq[i][q]
                  for i in range(4) for k in K for q in range(t)))
        npaf = asm.aggregate_npaf(Tq)
        nz = sum(1 for s in range(1, t) if npaf[s])
        check("t=%d: NPAF peak reads t (%d)" % (t, npaf[0]), npaf[0] == t)
        print("      [--] recorded, NOT a failure: the aggregate NON-PERIODIC "
              "autocorrelation is\n           nonzero at %d of the %d positive "
              "lags, so this is a T-MATRIX quadruple,\n           exactly "
              "Cooper-Wallis's hypothesis, and not an aperiodic T-sequence."
              % (nz, t - 1))
        quads[t] = Tq

    # ------------------------------------------------------------------- [B]
    print("\n[B] the Williamson quadruples")
    bank = asm.load_bank("williamson-9-13.json")
    check("bank covers exactly w in %s" % (list(WS),),
          sorted(int(k) for k in bank["quadruples"]) == sorted(WS))
    for w in WS:
        W = asm.load_W(bank, w)
        a = asm.aggregate_paf(W)
        rs = [sum(r) for r in W]
        check("w=%-2d: four +-1 rows, symmetric circulant, "
              "sum_q PAF_q = %d*delta_0, sum (row sum)^2 = %d, row sums %s"
              % (w, a[0], sum(v * v for v in rs), rs),
              all(v in (1, -1) for r in W for v in r)
              and all(r[j] == r[(-j) % w] for r in W for j in range(w))
              and a[0] == 4 * w and not any(a[1:])
              and sum(v * v for v in rs) == 4 * w
              and rs == list(bank["quadruples"][str(w)]["row_sums"]))

    # ------------------------------------------------------------------- [C]
    print("\n[C] openness sweep -- Cati-Pasechnik Table 4 v2 (2025-08-30),")
    print("    parsed at run time from data/cp-table4-v2.json")
    table, prov = asm.load_cp_table4()
    print("    source: %s" % prov["source"])
    check("the transcription loads and holds all %d entries, every one at an "
          "odd n <= %d" % (len(table), prov["range"]["n_max"]),
          len(table) == asm.CP_ENTRY_COUNT)
    check("reading convention: n(m) means 2^m*n is the smallest recorded known "
          "order, so H(4n) is unresolved there; an odd n <= %d ABSENT from the "
          "file has m(n) = 2, i.e. H(4n) was already known"
          % asm.CP_N_MAX, True)

    print("\n    two controls: the same file, parsed the same way, must "
          "reproduce the\n    unresolved sets of the two t this repository "
          "actually closes orders at.")
    for t_c, expect in sorted(CONTROL_W.items()):
        got = asm.cp_open_w(table, t_c)
        check("CONTROL t=%d: unresolved exactly at w in %s -- %s"
              % (t_c, got, CONTROL_NOTE[t_c]), got == expect)

    print("\n    the %d-cell sweep   (order = 4*t*w ; n = t*w)"
          % (len(TS) * len(WS)))
    print("    %-5s %-4s %-7s %-8s %s" % ("t", "w", "n", "order", "verdict"))
    print("    " + "-" * 70)
    reachable, unverifiable = [], []
    for t in TS:
        for w in WS:
            n, order = t * w, 4 * t * w
            if n > asm.CP_N_MAX:
                v = "UNVERIFIABLE -- n > %d, past the end of Table 4" \
                    % asm.CP_N_MAX
                unverifiable.append((t, w, order))
            elif n in table:
                m = table[n]
                v = ("UNRESOLVED -- Table 4 %d(%d), smallest recorded 2^%d*%d "
                     "= %d" % (n, m, m, n, (1 << m) * n))
                reachable.append((t, w, n, m, order))
            else:
                v = "known -- n absent from Table 4, so m(n) = 2"
            print("    %-5d %-4d %-7d %-8d %s" % (t, w, n, order, v))
    known = len(TS) * len(WS) - len(reachable) - len(unverifiable)
    check("every one of the %d cells is classified: %d known, %d unverifiable, "
          "%d reachable-unresolved"
          % (len(TS) * len(WS), known, len(unverifiable), len(reachable)),
          known + len(unverifiable) + len(reachable) == len(TS) * len(WS))
    check("REACHABLE UNRESOLVED ORDERS: %d -- these three quadruples open "
          "nothing with the Williamson orders this repository holds"
          % len(reachable), not reachable)
    print("    UNVERIFIABLE (beyond the end of the table): %s"
          % ["4*%d*%d = %d" % (t, w, o) for t, w, o in unverifiable])

    # Is the emptiness an artefact of this repository's small Williamson bank?
    # Sweep EVERY odd w <= 63 -- the range in which Williamson-type quadruples
    # are known for all w -- regardless of what is banked here.
    wide = [(t, w, t * w, table[t * w]) for t in TS for w in asm.cp_open_w(table, t)]
    check("no odd w <= 63 AT ALL puts t*w on a Table 4 entry, for any of "
          "t = %s -- so the empty result is not an artefact of this "
          "repository's Williamson bank (hits: %s)" % (list(TS), wide),
          not wide)
    print("    Label: REPORTED-FROM-AUDITED-TABLE. The sweep is re-derived "
          "here from the\n    transcription in data/, but a table records what "
          "its authors knew on its\n    date; this is a reading of that table, "
          "not a non-existence proof. The six\n    matrices built below are "
          "calibration gates at known orders, not results.")

    # ------------------------------------------------------------------- [D]
    print("\n[D] calibration gates")
    print("    w=1 is degenerate in the Z_t x Z_w product; the w=3 gates are "
          "the ones that\n    exercise the group structure.")
    keys = sorted(GATE_PINS)
    if args.quick:
        skipped = [GATE_PINS[k][0] for k in keys
                   if GATE_PINS[k][0] > asm.QUICK_MAX_ORDER]
        keys = [k for k in keys if GATE_PINS[k][0] <= asm.QUICK_MAX_ORDER]
        print("    --quick: skipping orders > %d: %s"
              % (asm.QUICK_MAX_ORDER, skipped))
    built = []
    for (t, w) in keys:
        order, pin = GATE_PINS[(t, w)]
        tag = "gate t=%d x w=%d -> H(%d)" % (t, w, order)
        print("\n    %s" % tag)
        ok, sha = asm.build_and_pin(quads[t], asm.load_W(bank, w), t, w,
                                    order, pin, tag)
        built.append((t, w, order, sha, ok))

    print("\n[D'] negative control")
    asm.tamper_reject(quads[101], asm.load_W(bank, 1), 101, 1, 404)

    # ------------------------------------------------------------------- [E]
    print("\n[E] verdicts")
    for t, w, order, sha, ok in built:
        print("VERDICT H(%-5d) = 4*%d*%-2d  GATE (known order)  |  %s  sha=%s"
              % (order, t, w, "verify.py-green" if ok else "FAILED", sha))
    for t in TS:
        print("VERDICT T-matrix quadruple of order %-3d  multiplier group "
              "order %d  |  PROVEN-BY-CERTIFICATE"
              % (t, PROFILE[t]["K_order"]))
    print("        T-sequences of lengths 101 and 113 are in the published "
          "literature, so\n        those two witnesses are independent "
          "re-derivations and no firstness is\n        claimed for them. For "
          "181 the claim ceiling is: the first publicly\n        accessible "
          "and independently reproducible order-181 T-matrix witness\n        "
          "located in our audit (audit closed 2026-08-29).")

    print("\n" + "=" * 78)
    print("-- elapsed %.1fs --" % (time.time() - t0))
    if asm.FAIL:
        print("CERT 04: FAIL (%d): %s" % (len(asm.FAIL), asm.FAIL))
        return 1
    print("CERT 04 GREEN%s" % (" (--quick)" if args.quick else ""))
    print("=" * 78)
    return 0


if __name__ == "__main__":
    sys.exit(main())
