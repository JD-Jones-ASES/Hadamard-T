# Cert 03 — a T-matrix quadruple of order 163, and H(7172) + H(8476)

## Statement

**1. The object.** `data/T163.json` holds four `(0, ±1)` sequences `T₁…T₄` on
`Z₁₆₃` with pairwise disjoint supports of sizes `(73, 36, 27, 27)`
partitioning `Z₁₆₃`, row sums `(1, 0, 9, 9)` with `Σ rᵢ² = 163`, satisfying
the **periodic** identity

```
Σᵢ PAF(Tᵢ)(s) = 163·δ₀(s),      PAF(x)(s) = Σ_q x[q]·x[(q+s) mod 163]
```

— a T-matrix quadruple of order 163. It is strictly invariant under
`K = {1, 38, 40, 53, 58, 85, 104, 133, 140}`, the unique subgroup of order 9
of `Z₁₆₃*`, which `run.py` re-derives by closure rather than reading in. Its
aggregate non-periodic autocorrelation is nonzero at 132 of the 162 positive
lags: a T-matrix quadruple, not an aperiodic T-sequence.

**2. Two Hadamard matrices,** through the Cooper–Wallis product with
symmetric circulant Williamson quadruples.

| order | `= 4·163·w` | canonical SHA-256 | file |
| --- | --- | --- | --- |
| 652 (gate) | `4·163·1` | `f488e58cad98e11624fba4c113fbd34b3a77bb1099f96f66772918db066b1c1a` | 0.4 MB |
| **7172** | `4·163·11` | `8535de5debeae63611d22d39fc3a4b821404dd4320699202b8f17c6851e46aa2` | 51.4 MB |
| **8476** | `4·163·13` | `25286dd5c3b1500ad6db0ebfdfb7aa26067dd3a4c779137dbdc18dab8a8dd77c` | 71.9 MB |

Each digest is pinned three ways: the run's own streaming digest,
`verify.py`'s `canonical_sha256=` verdict field, and the constant in `run.py`.
The matrices are not committed; they regenerate from `data/` in about a
minute.

## Honesty labels

| claim | label |
| --- | --- |
| the quadruple and the two matrices exist and are reproducible byte for byte | PROVEN-BY-CERTIFICATE |
| the openness of 7172 and 8476 | REPORTED-FROM-AUDITED-TABLE |
| the standing of the order-163 object | hedged, dated — see below |

The claim ceiling for the object, and nothing stronger:

> the first publicly accessible and independently reproducible order-163
> T-matrix witness located in our audit (audit closed 2026-08-29)

Openness: Cati–Pasechnik Table 4 v2 (arXiv:2411.18897v2, 2025-08-30) lists
`7172(1793,4)` and `8476(2119,4)` per this laboratory's machine-diffed
transcription; re-check at release date. `m = 4` means neither `2²·n` nor
`2³·n` carries an entry there, so the same two matrices also give **H(14344)**
and **H(16952)** by Sylvester doubling, under the same label. Cert 04 re-derives
the same unresolved set at `t = 163` — exactly `w ∈ {11, 13}` — from that
transcription at run time.

The Williamson quadruples of orders 1, 11 and 13 are classical, replay-grade
objects and carry no novelty claim. The `w = 13` quadruple used here is the
one in `data/williamson-13.json`; it differs from the `w = 13` quadruple in
`data/williamson-9-13.json`, and the pinned digests depend on that choice.

## How to run

```
python verify/verify.py --selftest              # the trust chain
python certs/03-t163-products/run.py            # ~1 min
python certs/03-t163-products/run.py --quick    # gate only, no order > 1200
```

Standard library only, from the repository root. One matrix is written to a
scratch file at a time and deleted after `verify.py` returns. Peak scratch
usage 72 MB, at H(8476).

## Expected verdict lines

```
VERDICT: HADAMARD order=652 all 212226 row pairs orthogonal
canonical_sha256=f488e58cad98e11624fba4c113fbd34b3a77bb1099f96f66772918db066b1c1a
VERDICT: HADAMARD order=7172 ... canonical_sha256=8535de5debeae636…51e46aa2
VERDICT: HADAMARD order=8476 ... canonical_sha256=25286dd5c3b1500a…8a8dd77c
CERT 03 GREEN
```

A negative control must also pass: flipping one character of the green H(652)
makes `verify.py` exit 1.

## Scope

This certificate binds the existence of the order-163 quadruple and of the two
matrices, and nothing more. It says nothing about other multiplier tiers at
163, which were not searched. On status: Table 4 v2 lists the entries
`1793(4)` and `2119(4)` unresolved; this certificate proves the two matrices
exist; we do not claim the table is a non-existence proof, and we do not need
it to be.
