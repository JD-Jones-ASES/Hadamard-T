# Hadamard-T

T-matrix witnesses, and the Hadamard orders they close. This repository holds
the publicly posted Hadamard matrix of order 2060, re-verified and decoded here into
a periodic T-matrix quadruple of order 103; five further T-matrix quadruples, at
orders 101, 113, 163, 181, and 131 (the last obtained twice, by two independent
routes); and thirteen Hadamard orders built from those quadruples by the
Cooper–Wallis product, each accepted by the trust chain in exact integer
arithmetic. Every constructed matrix and every banked object replays from
`certs/`; search history is recorded and labelled in the note, not replayed.

## Definitions

**T-matrix quadruple of order t.** Write PAF(x)(s) = Σ_q x[q]·x[(q+s) mod t]
for x : Z_t → Z. A T-matrix quadruple of order t is four sequences
T₁,…,T₄ : Z_t → {0, +1, −1} with pairwise disjoint supports covering Z_t —
exactly one Tᵢ is nonzero at each of the t positions — satisfying
Σᵢ PAF(Tᵢ)(s) = t·δ₀(s). The four circulant matrices carrying these sequences
as first rows are T-matrices in the sense of Cooper–Wallis.

**The product** (Cooper–Wallis 1972). Such a quadruple together with a
Williamson-type quadruple of order w yields a Hadamard matrix of order 4tw.

**Periodic, not aperiodic.** The order-103 object is a periodic T-matrix
quadruple: an aperiodic T-sequence of length 103 is not constructed here and is
not claimed. The periodic identity is exactly the Cooper–Wallis hypothesis, and
it is what the products below need.

## Claims

Thirteen Hadamard orders — 2884, 3668, 4532, 5764, 6812, 7004, 7172, 7828,
8476, 8908, 9476, 9956 and 11948 — are constructed and machine-verified here,
and every one of them lands on an entry listed unresolved in Cati–Pasechnik
Table 4 v2 (arXiv:2411.18897v2, 2025-08-30), in a transcription taken from the
arXiv LaTeX e-print, machine-diffed 195/195 on 2026-08-29 and re-checked at
release date. Eight of the thirteen — 2884, 4532, 7004, 7172, 7828, 8476, 9476
and 11948, the products of the order-103 and order-163 quadruples — carry a
further result: the firsthand prior-art audit closed 2026-08-29 located no
published route to the order-103 or order-163 T-objects they are built from.

On the six orders built from the order-103 quadruple:

> Hadamard matrices of orders 2884, 4532, 7004, 7828, 9476 and 11948,
> constructed and machine-verified here. The Cati–Pasechnik database
> (arXiv:2411.18897v2, 2025-08-30) records no Hadamard matrix of order 2²·n for
> the corresponding n; on that basis these are the first publicly accessible and
> independently reproducible constructions of these orders located in our
> audit. A database records what its authors knew, not what exists.

On the two built from the order-163 quadruple: `H(7172)` and `H(8476)` land on
the entries `1793(4)` and `2119(4)`, and `m = 4` there records that neither
`2²·n` nor `2³·n` was known, so the same two matrices also give `H(14344)` and
`H(16952)` by Sylvester doubling under the same reading.

The discrimination between the eight and the other five is itself a result of
the audit, and it runs both ways. The five orders built from the order-131
quadruple were separated out precisely because a **published route closes
them**: London and Kotsireas 2025 supply `TT(44)`, and the classical chain
`TT(44) → BS(87, 44) → TS(131)` reaches order 131, so no novelty of existence
is claimed at any of the five. For the eight, the same audit — eleven primary
sources read firsthand on 2026-08-29 — located no such route. The arithmetic
reasons are recorded route by route in [`note/NOTE.md`](note/NOTE.md) §4 and
§6.1: base sequences at 103 would need `BS(52, 51)`, past the verified frontier
at `n ≤ 43`; the Turyn-type ladder gives lengths `3n − 1` and steps over 103;
103 and 163 are prime, so the composite routes need a factorisation they do not
have; and the `6m + 1` method family, in range at both orders, prints nothing
above order 79 directly, or 113 by one-off cyclotomic search.

Scope of the table reading. A database records what its authors knew, not what
exists. The status claims above are documentary — the named, dated table lists
these entries unresolved — and they are checkable twice over: `certs/04`
re-derives that reading from the shipped transcription at run time, and the
matrices themselves are proved by certificate independently of any table.

### Priority postures

On the order-103 object, we make no priority claim:

> The order-103 T-matrix quadruple used here was decoded from the publicly
> posted Hadamard matrix of order 2060 (gist, 2026-08-23). Our audit
> (2026-08-29) located no earlier explicit or independently reproducible
> order-103 T-matrix in the literature; we make no priority claim on the object
> itself.

At orders 101, 113, and 131 we claim no firstness of any kind:

> T-sequences of these lengths are in the published literature
> (Best–Đoković–Kharaghani–Ramp 2013 and Zuo–Xia–Xia 2013 for 113; the
> Turyn-type ladder for 101/119; London–Kotsireas 2025 for 131). Our witnesses
> are independent re-derivations, not new objects.

At orders 163 and 181, the ceiling wording is:

> the first publicly accessible and independently reproducible order-163
> [resp. order-181] T-matrix witness located in our audit (audit closed
> 2026-08-29)

Each posture is fixed by the audit of 2026-08-29 and is carried unchanged in
the note, in the `data/` records, and in the certificate that replays the
object.

## The orders

An entry n(m) means m is minimal with a Hadamard matrix of order 2ᵐ·n known in
Table 4 v2, so an entry with m ≥ 3 is an order 4n the table does not record.
Matrix counts are a remark, not the scoreboard: the order-131 row is built
twice, both builds in cert 05.

| t | orders 4tw, with w in parentheses | Table 4 entries n(m) |
| --- | --- | --- |
| 103 | 2884 (7), 4532 (11), 7004 (17), 7828 (19), 9476 (23), 11948 (29) | 721(3), 1133(3), 1751(3), 1957(3), 2369(3), 2987(3) |
| 163 | 7172 (11), 8476 (13) | 1793(4), 2119(4) |
| 131 | 3668 (7), 5764 (11), 6812 (13), 8908 (17), 9956 (19) | 917(3), 1441(3), 1703(3), 2227(3), 2489(3) |

## Replay

Everything runs from the repository root on bare Python 3.9 or newer. Standard
library only, no network.

```
python verify/verify.py --selftest
python certs/01-h2060-artifact/run.py
python certs/02-t103-products/run.py
python certs/03-t163-products/run.py
python certs/04-t-bank/run.py
python certs/05-t131-products/run.py
```

`verify/verify.py` is the trust chain. It accepts a matrix file only if the
matrix is square, has every entry in {+1, −1}, and satisfies H·Hᵀ = n·I, by
exact integer arithmetic on packed rows with no floating point anywhere. Each
cert rebuilds its matrices from the banked data, hands them to the trust chain,
compares the canonical SHA-256 in the verdict against the digest pinned in its
`NOTES.md`, and deletes them. Generated matrices go to `out/` or to a temporary
file and are never committed.

## Layout

| path | contents |
| --- | --- |
| `note/NOTE.md` | the note: definitions, statements, and each claim's cert |
| `verify/verify.py` | the trust chain: Hadamard check, exact arithmetic, stdlib only |
| `tools/assemble.py` | the shared Goethals–Seidel and Cooper–Wallis assembler |
| `data/` | banked quadruples, the Williamson bank, the order-2060 artifact, provenance pins |
| `certs/` | one directory per claim, each with `run.py` and `NOTES.md` |

## Credits, license, provenance

The order-2060 payload is Schneider's, posted publicly on 2026-08-23; it was
verified and decoded here. The Turyn-type sequences behind one of the two order-131 routes
are London and Kotsireas 2025; prior T-sequences at these lengths are due to
Đoković and to London. The machinery is classical: Williamson; Baumert and
Hall; Goethals and Seidel; Turyn; Cooper and Wallis. The audited status table
is Cati and Pasechnik, arXiv:2411.18897v2 (2025-08-30).

Code is MIT; see [LICENSE](LICENSE). The note (`note/`) and the prose
documentation in `data/` are CC BY-SA 4.0; see [LICENSE-DOCS](LICENSE-DOCS.md)
for the exact boundary. The banked matrices and sequences are mathematical
data: no license is claimed over them, and the redistributed public artifact
keeps its provenance pins. Authorship is in [DISCLOSURE.md](DISCLOSURE.md).

Companion repository: [Hadamard-M](https://github.com/JD-Jones-ASES/Hadamard-M),
the erratum on Miyamoto's 1991 construction at order 515 and an explicit
Hadamard matrix of order 7796. Formalization:
[Hadamard-formal](https://github.com/JD-Jones-ASES/Hadamard-formal) — a Lean 4
/ Mathlib development proving the Cooper–Wallis construction and the existence
theorems at orders 2884, 4532, 7004, 7172, 7828, 8476, 9476 and 11948 from
kernel-checked copies of this repository's banked witnesses, together with the
parity obstruction, its empty-class corollary, and the C₂-display
nonexistence theorem behind Hadamard-M's erratum. That development is
registered in the Palomar Registry as PALOMAR-2026-08-31-000001.
