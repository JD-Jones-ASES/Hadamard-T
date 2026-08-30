# Cert 02 — six Hadamard matrices from the order-103 T-matrix quadruple

## Statement

The order-103 T-matrix quadruple certified in
[cert 01](../01-h2060-artifact/NOTES.md), together with symmetric circulant
Williamson quadruples of order `w` held in `data/`, yields through the
**Cooper–Wallis product** explicit Hadamard matrices of six orders `4·103·w`.

| order | `w` | canonical SHA-256 | file |
| --- | --- | --- | --- |
| **2884** | 7 | `38c778c1de4d5a9af4e7a9391382ee8b98ca89143a4885efff0d013429ebff8e` | 8.3 MB |
| **4532** | 11 | `3f00bb1954cc7fa77e2112ea147a76b59de18142fe17b28348d3c266f21a7348` | 20.5 MB |
| **7004** | 17 | `50b9c71d6683a3bbf7bc1224243c8903932a1919d9683a91de49377433b707f0` | 49.1 MB |
| **7828** | 19 | `6ea2946fe2ba73e5e06cf77fda33e99004a885228844e5da6deea96cb05286c5` | 61.3 MB |
| **9476** | 23 | `ea4c04d1759029421df7572450365da4a16c466b9d435dfcc79930a43e548f96` | 89.8 MB |
| **11948** | 29 | `8169db47b795b30cfc6036f4b090380902dd06f9e72726f0f9e2fe8120ed550b` | 142.8 MB |

Two calibration gates run the same pipeline end to end at small orders:

| gate | `t × w` | canonical SHA-256 |
| --- | --- | --- |
| H(36) | 3 × 3 | `25f0a20be9d4ae4f29166241f2c324b772303740ff5c5538975bc17f1ec17fa9` |
| H(412) | 103 × 1 | `fd7d2e174d3b52fecdc20cce963a3c278b00d07bc392018e0b8ded841991512e` |

The matrices are not committed. They are regenerated deterministically from
`data/` in about three minutes, and the digests above are what binds them.
Each digest is pinned three ways: the run's own streaming digest,
`verify.py`'s `canonical_sha256=` verdict field, and the constant in `run.py`.

The product itself, with the sign table and the Goethals–Seidel orientation
written out, is documented at the top of
[`tools/assemble.py`](../../tools/assemble.py). The order-103 object is a
T-matrix quadruple in the periodic sense; the terminology paragraph in cert 01
applies here unchanged.

## Honesty labels

| claim | label |
| --- | --- |
| the six matrices exist and are reproducible byte for byte | PROVEN-BY-CERTIFICATE |
| the openness of the six orders | REPORTED-FROM-AUDITED-TABLE |

> Hadamard matrices of orders 2884, 4532, 7004, 7828, 9476 and 11948,
> constructed and machine-verified here. The Cati–Pasechnik database
> (arXiv:2411.18897v2, 2025-08-30) records no Hadamard matrix of order 2²·n
> for the corresponding n; on that basis these are the first publicly
> accessible and independently reproducible constructions of these orders
> located in our audit. A database records what its authors knew, not what
> exists.

The Table 4 entries the six orders land on are `721(3)`, `1133(3)`,
`1751(3)`, `1957(3)`, `2369(3)`, `2987(3)`, quoted from the machine-diffed
transcription of the v2 e-print (2025-08-30) and re-checked at release date;
cert 04 re-derives the same unresolved set at `t = 103` from that transcription
at run time. The Williamson quadruples are classical, replay-grade objects and
carry no novelty claim.

## How to run

```
python verify/verify.py --selftest              # the trust chain
python certs/02-t103-products/run.py            # all six, ~3 min
python certs/02-t103-products/run.py --quick    # gates only, no order > 1200
```

Standard library only, from the repository root. One matrix is written to a
scratch file at a time and deleted after `verify.py` returns; the row cache is
freed before the subprocess starts, so builder and verifier never hold a large
object at once. Peak scratch usage 143 MB, at H(11948).

## Expected verdict lines

```
VERDICT: HADAMARD order=2884 ... canonical_sha256=38c778c1de4d5a9a…29ebff8e
VERDICT H(2884 ) = 4*103*7   |  PROVEN-BY-CERTIFICATE  sha=38c778c1…
...
CERT 02 GREEN
```

Two negative controls must also pass: perturbing one entry of the `w = 3`
Williamson quadruple aborts the build before assembly, and flipping one
character of the green H(36) makes `verify.py` exit 1.

## Scope

This certificate binds the existence of the six matrices and nothing more; the
status of the six orders is reported from an audited table, under the label
above. Table 4 v2 lists the six entries unresolved; this certificate proves the
six matrices exist; we do not claim the table is a non-existence proof, and we
do not need it to be.
