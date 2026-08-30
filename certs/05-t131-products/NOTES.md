# Cert 05 — T-sequences of order 131, and five Hadamard matrices

## Statement

Two inequivalent T-matrix quadruples of order 131 are certified here, and each
builds, through Cooper–Wallis with the banked Williamson quadruples, all five
of

```
H(3668) = 4·131·7    H(5764) = 4·131·11    H(6812) = 4·131·13
H(8908) = 4·131·17   H(9956) = 4·131·19
```

— ten matrices in all, every one accepted by `verify/verify.py` at a digest
pinned three ways.

**Witness 1, the published route.** `TT(44) → BS(87,44) → TS(131)`. Both
arrows are re-derived here, not quoted. If `(A,B,C,D)` are Turyn-type of
lengths `n, n, n, n−1`, then `(C;D, C;−D, A, B)` are base sequences of lengths
`2n−1, 2n−1, n, n`; the cross terms of the two concatenations cancel, so
arrow 1 is the Turyn-type identity re-bracketed. From base sequences
`(P,Q,X,Y)` of lengths `m, m, n, n`,

```
T₁ = (P+Q)/2 ; 0_n    T₂ = (P−Q)/2 ; 0_n
T₃ = 0_m ; (X+Y)/2    T₄ = 0_m ; (X−Y)/2
```

are `(0, ±1)` with disjoint supports covering `Z_{m+n}`, and
`Σ N(Tᵢ) = ½ Σ N(base) = (m+n)·δ₀`. The result has supports `(44, 43, 23, 21)`,
row sums `(8, 3, −3, −7)`, **vanishing aggregate non-periodic**
autocorrelation — genuine T-sequences — and trivial multiplier group. It
equals `data/T131-literature-route.json` character for character.

**Witness 2, this laboratory's own.** A `K`-invariant quadruple from a
multiplier-tier search: supports `(11, 10, 55, 55)`, row sums `(−9, 0, 5, 5)`,
multiplier group of order 5, aggregate non-periodic autocorrelation nonzero
off-peak. The tier lattice is recomputed from scratch in `run.py`: of the
eight divisors of 130, six admit no row-sum profile even solving `Σ rᵢ² = 131`,
so only `d ∈ {1, 5}` survive, and `d = 5` carries exactly one admissible
profile up to the per-class sign flip, `|r| = (0, 5, 5, 9)`.

**The two are inequivalent** — the support multisets `(44, 43, 23, 21)` and
`(11, 10, 55, 55)` differ, and the support multiset is preserved by every
equivalence move on T-matrix quadruples (reorder, negate, shift, decimate);
`run.py` asserts that separation. The row-sum, strict-multiplier-group and
aggregate-NPAF differences are recorded as descriptions of the banked
representatives, not used as invariants — the strict stabiliser and the NPAF
are not preserved under cyclic shift. The ten matrices come in
five pairs with distinct canonical digests.

| order | `w` | witness 1 canonical SHA-256 | witness 2 canonical SHA-256 |
| --- | --- | --- | --- |
| 524 (gate) | 1 | `97aca328b00f3715fcce2a7e15c8637da09d279e383a5c3d41b6587be60d1205` | `003385ed9d3e528cea9fdfb85835d84b6cd92169f2062742547f43f48892851f` |
| 1572 (gate) | 3 | `4ce087b59bf79e0bfaf1d99c13c9f1a39c397ba24eb518b6e85a6fa9cb908196` | — |
| **3668** | 7 | `0dee08275687490a49df9c43bfc276ec93f95107365d05e3f3b48cfa31991ddf` | `defab1069bc0af50b5ec85dc913818f1d0763c69cf9c3eaefbedefb4a33cd90a` |
| **5764** | 11 | `b391b9bad54b235f6f1d9c778cd72256760a46f6e04cb26f2193e016d326ba02` | `3cbaf9fbcb28f601ef490e4b631961b4e5d60f463069188a0f830f265881d7cd` |
| **6812** | 13 | `eec46e29ad9b7933053c5e62912fb01c4113289e167a8b0ecfc9d9bb287a9d42` | `4fa5d4d6961eaa699590c5f45264cda2eaf6f42f1f8fd8c60e1ee899138f368b` |
| **8908** | 17 | `0107ba27c6c10f4ff65011f79e4e564a56a7ab4359535a634dc8c4b6bb2b9afe` | `88c6d5e9cd97eca3e3405779b68fb55f7c1d902d74719a4e8dfffc3979fdc0af` |
| **9956** | 19 | `a14f4b6991fbad498cebc568383119a0f4afc38eafc4be3669bf503619280496` | `2799b6154c4768b1b2df72b8ab472ddf6565f54c7d4d2d99fde83afbead41387` |

The `w = 13` quadruple used here is the one in `data/williamson-9-13.json`;
the pinned digests depend on that choice.

## Honesty labels

| claim | label |
| --- | --- |
| the two order-131 quadruples exist and satisfy the T-matrix axiom | PROVEN-BY-CERTIFICATE |
| the ten Hadamard matrices exist | PROVEN-BY-CERTIFICATE |
| the TT(40), TT(42), TT(44) payloads are London and Kotsireas's data | VERIFIED-PRIMARY (publisher PDF read firsthand 2026-08-29; the payloads are byte-identical to both printed copies) |
| the five orders are listed at `m = 3` in Cati–Pasechnik Table 4 v2 | REPORTED-FROM-AUDITED-TABLE (machine-diffed transcription of the 2025-08-30 e-print; re-check at release date) |
| "no one had these matrices before" | **not claimed** |

**No novelty of existence is claimed at 131.** These five orders are
deliberately separated from the eight orders of certs 02 and 03, and the
separation is the audit's own finding: at 131 the audit of 2026-08-29 located a
published route that closes them, and at 103 and 163 the same audit located no
such route.

> T-sequences of these lengths are in the published literature
> (Best–Đoković–Kharaghani–Ramp 2013 for 113; the Turyn-type ladder for
> 101/119; London–Kotsireas 2025 for 131). Our witnesses are independent
> re-derivations, not new objects.

A published route closes the five orders: London and Kotsireas 2025 supply
TT(44), and both compositions above are classical, as is Cooper–Wallis. Table 4
v2 simply does not record it. What is claimed is:

> the first publicly accessible and independently reproducible construction of
> H(3668), H(5764), H(6812), H(8908) and H(9956) located in our audit of
> 2026-08-29 — explicit matrices, rebuildable by anyone, with a trust chain
> of about 250 lines of standard-library Python.

The order-5 multiplier-tier witness is a separate matter: it is this
laboratory's own object. Whether a T-matrix quadruple of order 131 with a
five-element multiplier group is on record anywhere is **unaudited**. That
tier is **not** closed either: the search behind the banked snapshot was
partial, a full exhaustion is priced out, and a partial sweep is not a
nonexistence result.

`[A']` guards against the obvious failure mode, code tuned to make 44 work:
TT(6) and TT(8) are found *here*, by a complete search of the full `2^{4n−1}`
space with no external input, and pushed through the same two arrows to
TS(17) and TS(23). The composition code never sees the number 44.

## How to run

```
python verify/verify.py --selftest              # the trust chain
python certs/05-t131-products/run.py            # everything, ~4 min
python certs/05-t131-products/run.py --quick    # no order > 1200
```

Standard library only, from the repository root. One matrix is written to a
scratch file at a time and deleted after `verify.py` returns.

## Expected verdict lines

```
VERDICT: HADAMARD order=524 all 137026 row pairs orthogonal
canonical_sha256=97aca328b00f3715fcce2a7e15c8637da09d279e383a5c3d41b6587be60d1205
VERDICT: HADAMARD order=3668 ... canonical_sha256=0dee082756874904…31991ddf
VERDICT TS(131): four (0,+-1) sequences, supports [44, 43, 23, 21] partitioning Z_131,
        aggregate NPAF and PAF both = 131*delta_0  |  PROVEN-BY-CERTIFICATE
CERT 05 GREEN
```

A negative control must also pass: flipping one character of the green H(524)
makes `verify.py` exit 1.

## Scope

This certificate binds the existence of the two order-131 quadruples and of
the ten matrices, and nothing more; it makes no firstness claim at 131. On
status: Table 4 v2 lists the five entries unresolved; this certificate proves
the ten matrices exist; we do not claim the table is a non-existence proof, and
we do not need it to be.
