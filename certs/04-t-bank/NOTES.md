# Cert 04 — T-matrix quadruples of orders 101, 113 and 181

## Statement

**1. The objects.** `data/T101.json`, `data/T113.json` and `data/T181.json`
each hold four `(0, ±1)` sequences `T₁…T₄` on `Z_t` with pairwise disjoint
supports partitioning `Z_t`, satisfying the **periodic** identity
`Σᵢ PAF(Tᵢ)(s) = t·δ₀(s)`. `run.py` computes the **full stabiliser**
`{k ∈ Z_t* : Tᵢ[kq] = Tᵢ[q] ∀ i,q}` directly rather than reading `K` in, then
re-derives by closure that this stabiliser is exactly the unique subgroup of
`Z_t*` of that order.

| `t` | supports | row sums | `Σ rᵢ²` | `\|K\|` | `K` |
| --- | --- | --- | --- | --- | --- |
| **101** | `(21, 20, 20, 40)` | `(−1, 0, −10, 0)` | 101 | 5 | `{1, 36, 84, 87, 95}` |
| **113** | `(36, 35, 28, 14)` | `(−8, 7, 0, 0)` | 113 | 7 | `{1, 16, 28, 30, 49, 106, 109}` |
| **181** | `(82, 36, 36, 27)` | `(10, 0, 0, −9)` | 181 | 9 | `{1, 39, 43, 48, 62, 65, 73, 80, 132}` |

In all three the full stabiliser is exactly the design tier, with no
accidental enlargement. The aggregate non-periodic autocorrelation is nonzero
at 74/100, 90/112 and 150/180 positive lags: these are T-matrix quadruples,
not aperiodic T-sequences.

**2. They close no order, and `run.py` derives that rather than asserting it.**
The Cati–Pasechnik Table 4 transcription is published in this repository as
`data/cp-table4-v2.json` (195 entries, odd `n ≤ 2999`, the reading convention
in its own provenance block), and `run.py` **parses it at run time**. For every
`t ∈ {101, 113, 181}` and every odd `w ≤ 63`, the odd part `n = t·w` carries no
entry there, so `H(4tw)` was already known: **zero reachable unresolved
orders**. Five of the thirty banked cells — `4·113·29`, `4·181·17`, `4·181·19`,
`4·181·23`, `4·181·29` — have `n > 2999`, past the end of the table, and are
printed **UNVERIFIABLE**: not open, not closed. The sweep is run over every odd
`w ≤ 63`, not only the ten `w` this repository banks, so the empty result is
not an artefact of the Williamson bank.

**Two controls bind that parse.** A transcription that could not reproduce the
closures this repository already spends would not be trustworthy for a negative
result, so the same file is re-read at the two `t` where orders actually close:

| control | unresolved `w` found | expected |
| --- | --- | --- |
| `t = 103` | `5, 7, 11, 17, 19, 23, 29` | cert 02's **six** targets, **plus `w = 5`** — `n = 515`, `H(2060)`, closed before this repository was cut and so never a cert-02 target |
| `t = 163` | `11, 13` | cert 03's two targets, exactly |

`run.py` therefore builds **calibration gates only** — six known orders, run to
exercise the assembler end to end on each new `t`. `w = 1` is degenerate in the
`Z_t × Z_w` product, so a `w = 3` gate is run at each `t` as well; those are the
ones that exercise the group structure.

| order | `= 4·t·w` | canonical SHA-256 | file |
| --- | --- | --- | --- |
| 404 | `4·101·1` | `b9a1d492851c3b9e9a81c672647c89c6e4551f43fe9dc0a7ce9eecc5aad9f0ed` | 0.2 MB |
| 1212 | `4·101·3` | `36ffec1af9be598bde2ec6682cfead01e833b3a69d970c7cff4ebd8b093637eb` | 1.5 MB |
| 452 | `4·113·1` | `1ab70fee6221a39bf2fbf37e1cd10c7913153cbf00de8ee57e2a880c08ea167e` | 0.2 MB |
| 1356 | `4·113·3` | `2aa915a52363e8c7bb8f21d819ff6dbc97f715e55e0226e6cf067cf66b019474` | 1.8 MB |
| 724 | `4·181·1` | `a4d006afa086bf2a88e14f2bc13678b8cb85fefc058d74fbcd1581a27df96df2` | 0.5 MB |
| 2172 | `4·181·3` | `7f4660c332ccaa25a805058b801c4291a0b33ae8e58655b623b7a81174d5a4e1` | 4.7 MB |

These are gates, not results: all six orders were already known.

## Honesty labels

| claim | label |
| --- | --- |
| the three quadruples exist and satisfy the T-matrix axiom | PROVEN-BY-CERTIFICATE |
| the six gate matrices exist | PROVEN-BY-CERTIFICATE |
| firstness at 101 or 113 | **never claimed** |
| the standing of the order-181 object | hedged, dated — see below |
| "these three buy nothing" | REPORTED-FROM-AUDITED-TABLE — re-derived here from `data/cp-table4-v2.json`, but a table records what its authors knew on its date |

On 101 and 113:

> T-sequences of these lengths are in the published literature
> (Best–Đoković–Kharaghani–Ramp 2013 for 113; the Turyn-type ladder for
> 101/119; London–Kotsireas 2025 for 131). Our witnesses are independent
> re-derivations, not new objects.

On 181, the claim ceiling and nothing stronger:

> the first publicly accessible and independently reproducible order-181
> T-matrix witness located in our audit (audit closed 2026-08-29)

Three caveats on 181 are load-bearing. The map that lists 181 as unreachable
is **derived** from the published multiplication laws, not a published table.
The audit was **not** run against the direct T-matrix literature (Sawade 1985
and successors), which is where a prior order-181 T-matrix would sit. And the
Table 4 fingerprint that corroborates 103 and 163 is **absent** at 181: every
`181·w` is already known by other routes, so the table can neither corroborate
nor refute.

## How to run

```
python verify/verify.py --selftest        # the trust chain
python certs/04-t-bank/run.py             # ~5 s
python certs/04-t-bank/run.py --quick     # skips the three gates over 1200
```

Standard library only, from the repository root. One matrix is written to a
scratch file at a time and deleted after `verify.py` returns. Peak scratch
usage 4.7 MB, at H(2172).

## Expected verdict lines

```
VERDICT: HADAMARD order=404 all 81406 row pairs orthogonal
canonical_sha256=b9a1d492851c3b9e9a81c672647c89c6e4551f43fe9dc0a7ce9eecc5aad9f0ed
VERDICT H(404  ) = 4*101*1   GATE (known order)  |  verify.py-green  sha=b9a1d492…
VERDICT T-matrix quadruple of order 101  multiplier group order 5  |  PROVEN-BY-CERTIFICATE
CERT 04 GREEN
```

A negative control must also pass: flipping one character of the green H(404)
makes `verify.py` exit 1.

## Scope

This certificate binds the existence of the three quadruples and of the six
gate matrices, and nothing more. It makes no closure claim: no order it
touches was open. It says nothing about products other than Cooper–Wallis,
nothing about orders with `n > 2999`, and nothing about other multiplier tiers
at these three `t`.
