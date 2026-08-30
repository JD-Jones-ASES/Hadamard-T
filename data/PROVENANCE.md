# data/ — what each file is and where it came from

Every file here is an input to a certificate in `certs/`. No generated matrix
is stored: the certificates rebuild them and delete them.

## The order-2060 artifact

`H2060.txt.gz` is a byte-exact copy of the file `H2060.txt` from the public
gist

> https://gist.github.com/schneiderlo/b866a2ff2fcd93934f0db54cfa4069d0

gist id `b866a2ff2fcd93934f0db54cfa4069d0`, user `schneiderlo`, description
"Hadamard matrix of order 2060", created **2026-08-23**, retrieved
**2026-08-29**. No gist revision sha was recorded at retrieval; the payload
digest below is what binds the file.

| quantity | value |
| --- | --- |
| payload SHA-256 (decompressed `H2060.txt`) | `c7a145d86210740dd3f8ea21ca896a54d6916007a042638f17c8c47f097200f7` |
| SHA-256 of `H2060.txt.gz` as stored here | `499bdd1fdad3cc67282bc87737c41782ee97a05956617d682a915ba09c5d31dd` |

The file is already in `verify/verify.py`'s canonical serialisation, so its
payload digest and its canonical digest coincide. Cert `01-h2060-artifact`
re-checks both.

**Announcement, cited separately from the artifact.** The gist is paired with
the post https://x.com/modkin_mp/status/2091352950039785791 (Loïc Schneider,
`@modkin_mp`, 2026-08-23 02:32 UTC), which states that a Hadamard matrix of
order 2060 exists and shows a 2060 x 2060 pixel image. Neither the gist nor
the post discloses a construction method. The artifact was located and
downloaded directly from the gist URL and verified firsthand; nothing in this
repository depends on the post.

Credit for the order-2060 matrix is Schneider's. This laboratory decoded its
structure; it did not find it.

## The T-objects decoded from that artifact

| file | contents |
| --- | --- |
| `seed515.json` | the four length-515 `+/-` sequences that generate the gist matrix, the normalising shifts `(0, 206, 206, 206)`, the row sums, and the two pinned canonical digests |
| `T103.json` | the order-103 T-matrix quadruple, as `classes[q]` and `signs[q]` with `T_i[q] = signs[q]` if `classes[q] = i` else `0` |
| `T3-control.json` | the brute-forced order-3 control quadruple, same encoding, with the pinned H(60) digest |

Decoded firsthand from `H2060.txt.gz` alone. No table, paper or third-party
method disclosure was used, and no priority is claimed on the order-103
object (see `certs/01-h2060-artifact/NOTES.md`).

## Turyn-type payloads and the order-131 witnesses

`turyn-type-payloads.json` holds the hexadecimal payloads for TT(40), TT(42)
and TT(44) printed in

> S. London and I. S. Kotsireas, *New Turyn-type sequences*, Cryptography and
> Communications **17**(5) (2025) 1601–1610,
> doi:[10.1007/s12095-025-00829-z](https://doi.org/10.1007/s12095-025-00829-z),
> published 2025-08-11.

The paper was read firsthand from the publisher PDF on 2026-08-29. It prints
each payload twice, in §3 and in Appendix A, and the two copies agree; the
strings banked here are byte-identical to both. The decoding convention
recorded in the file — final hex digit right-aligned, carrying `X, Y, Z` only
— is the paper's own printed Appendix-A formula, and was also determined
independently here from the payloads alone. Cert `05-t131-products` decodes
and re-checks all three against the Turyn-type identity rather than trusting
the transcription.

`T131-literature-route.json` is the TS(131) that the two published
compositions produce from TT(44). `T131-multiplier-tier.json` is a second,
inequivalent order-131 witness found by this laboratory's own multiplier-tier
search; its `note` field states what that search did and did not settle.

## Objects found by this laboratory

| file | contents |
| --- | --- |
| `T163.json` | K-invariant T-matrix quadruple of order 163, `K` the unique order-9 subgroup of `Z_163*` |
| `T101.json`, `T113.json`, `T181.json` | K-invariant T-matrix quadruples of orders 101, 113, 181 |

Each carries its own `derivation` and `literature_status` fields. T-sequences
of lengths 101 and 113 are in the published literature, so those two witnesses
are independent re-derivations and no firstness is claimed for them. The
claim ceilings for 163 and 181 are quoted in the certificates.

## The Cati–Pasechnik Table 4 transcription

`cp-table4-v2.json` holds all **195** entries of Table 4 of

> M. Cati and D. V. Pasechnik, *A database of constructions of Hadamard
> matrices*, arXiv:[2411.18897](https://arxiv.org/abs/2411.18897)**v2**
> (2025-08-30).

transcribed from the arXiv LaTeX e-print and machine-diffed against it,
195/195 agreeing, on 2026-08-29. An entry `n(m)` means `m` is the least
exponent for which a Hadamard matrix of order `2ᵐ·n` is recorded known, so
`H(4n)` is not recorded known while the entry stands; the table covers odd
`n ≤ 2999`, and entries with `m = 2` are omitted because those orders are
already known. An odd `n ≤ 2999` **absent** from the file therefore has
`m(n) = 2`. The file carries that convention in its own `provenance` block.

Cert `04-t-bank` parses this file at run time to derive its openness sweep,
and re-reads it at `t = 103` and `t = 163` as two controls that must reproduce
the unresolved sets certs 02 and 03 spend. Certs 02, 03 and 05 check the
entries they quote against it. Nothing here treats the table as a
non-existence proof: it is a dated report of what its authors knew, and every
certificate that reads it says so.

## Williamson quadruples

| file | orders `w` |
| --- | --- |
| `williamson-quadruples.json` | 1, 3, 7, 11, 17, 19, 23, 29 |
| `williamson-13.json` | 13 |
| `williamson-9-13.json` | 9, 13 (the 13 here is a **different** quadruple from `williamson-13.json`) |

Symmetric circulant Williamson-type quadruples, brute-forced by this
laboratory. They are classical, replay-grade objects — Williamson type is
known in the literature for every odd `w <= 63` except 35, 47, 53 and 59 — and
carry no novelty claim. Every one is re-verified exactly by the certificate
that uses it. The two order-13 quadruples are kept separate because the pinned
digests of the products depend on which one is used.

## File digests

| file | bytes | SHA-256 |
| --- | --- | --- |
| `H2060.txt.gz` | 54711 | `499bdd1fdad3cc67282bc87737c41782ee97a05956617d682a915ba09c5d31dd` |
| `cp-table4-v2.json` | 3620 | `74bdb560ccfb63bafbd717270c43b030a5db9a8a42f1c4131142aa3d2b151362` |
| `T101.json` | 2055 | `4d40d2456e7bd6b759b39ca57aa19e64db2b37d0a7dc1bb1fa5bb4a6c6e9e3d1` |
| `T103.json` | 1989 | `cf59cb40fcdc74f9309161ad2d50ddd7ad462446da365efb527de0dd4e295a06` |
| `T113.json` | 2207 | `521ae6053b9fc20b99267ad8f2e94b36bb594910537d1f87a56ee19e8c3528cb` |
| `T131-literature-route.json` | 572 | `a5f290fabd8d4e426010c1a26527cf846bbda1281444eb384fea2d743b473a2d` |
| `T131-multiplier-tier.json` | 13903 | `4c91eba6b5c8edc1535760aac7b59cd4aa775608d3fbbeaeb6b7789151297e50` |
| `T163.json` | 2376 | `285557b7fc0f846af4fcec2e83d5d97eabcff3b82b2171bd0a6bdf32f8d27a4d` |
| `T181.json` | 3156 | `b1a737ddd4477b733dbd99a49aba22058f71df9437706efa37e854fcf98e6a75` |
| `T3-control.json` | 423 | `e803340908031150490d2515080e39263da93d5d97ecdd2344dcc12084de71f4` |
| `seed515.json` | 2923 | `9be8a0f0b637aa93a3071dd4bf4ece11767b236cc86022506567cc928467068d` |
| `turyn-type-payloads.json` | 2049 | `752f78c4b3ba188ac0031c98bb990e2f020f7f174a306590a1c4f421b265872b` |
| `williamson-13.json` | 770 | `249ba0919c5f0b82faab5d6c248a1a2ab18107d5a6ee2175195236fc3509935f` |
| `williamson-9-13.json` | 876 | `4621ea45635f7895beabbe01a692e76989f22b1abf5cb7b2141c784685e07810` |
| `williamson-quadruples.json` | 2390 | `515dac963269e31694445c5ad7635b001e16b3ea3cdb58fa8ee12eaae47e0446` |

## Redistribution

`H2060.txt.gz` is a redistributed copy of a publicly posted matrix, kept
byte-identical to the source payload and pinned by the digests above; the
matrix is mathematical data and no license is claimed over it here. The
Table 4 transcription is a transcription of a published table, credited to
its authors. License boundaries for this repository are in LICENSE-DOCS.md.
