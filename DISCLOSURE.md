# Disclosure

The results are thirteen Hadamard orders — 2884, 3668, 4532, 5764, 6812, 7004,
7172, 7828, 8476, 8908, 9476, 9956 and 11948 — constructed by the Cooper–Wallis
product and machine-verified in exact integer arithmetic on entries listed
unresolved in Cati–Pasechnik Table 4 v2 (arXiv:2411.18897v2, 2025-08-30),
together with the six periodic T-matrix quadruples they are built from, at
orders 101, 103, 113, 131, 163 and 181.

AI-generated results with a human managing the workflow. Produced by Claude
Code (Fable 5, Anthropic); external reviews by GPT 5.6 (OpenAI) and by
Grok (xAI), intaken and adjudicated.

## What the AI stations did

Everything mathematical. The decode of the public order-2060 matrix into a
Goethals–Seidel array and then into a periodic T-matrix quadruple of order 103;
the searches that produced the quadruples at orders 101, 113, 163, 181 and 131;
the re-derivation of the literature route to order 131 rather than its
quotation; the Cooper–Wallis assembly of the thirteen Hadamard orders; the
verifier and every certificate; the prior-art audit that fixes the claim
wordings in the README and the note. Outside stations reviewed the claims and
the hedging; their reports were relayed in and adjudicated here
against primary sources; where a report and a source disagreed, the source
governed.

## What the human owner did

Granted the sessions, paid for the compute, obtained papers this laboratory
could not otherwise read, relayed material to and from the outside reviewers,
and ruled on publication, licensing, and scope. No mathematical contribution,
and none is claimed. The owner's name appears here and in the copyright line;
it appears in no derivation.

## Verification

The trust chain is `verify/verify.py`. It accepts a matrix only if the matrix is
square, has every entry in {+1, −1}, and satisfies H·Hᵀ = n·I. The arithmetic is
exact: rows are packed into integers and orthogonality of a pair is a popcount
identity, so no floating point enters anywhere. `python verify/verify.py
--selftest` exercises it against known Hadamard matrices, against
Hadamard-preserving row and column operations, and against corruptions it must
reject.

Every claim in the note is replayable. Each directory under `certs/` rebuilds
its objects from the small banked data in `data/`, re-establishes the defining
identities with its own exact-integer loops, hands the resulting matrices to the
trust chain, and checks the canonical SHA-256 in the verdict against the digest
pinned in that certificate. The large matrices are regenerated, not committed.
Nothing outside the standard library is used, and nothing requires the network.

Reading the certificates is one route. Rebuilding from the definitions in the
note is another, and it requires trusting none of this. Independent verification
is invited; questions and verification reports are welcome via GitHub issues.

## Credit for external mathematics

- **Schneider** — the Hadamard matrix of order 2060 posted publicly on
  2026-08-23. It is their artifact. This laboratory verified a copy and decoded
  it; the order-103 quadruple carries no priority claim here.
- **London and Kotsireas (2025)** — the Turyn-type sequences at length 44,
  behind one of the two independent routes to order 131.
- **Đoković, and London** — prior T-sequences at the lengths banked here.
  The witnesses at orders 101, 113 and 131 are independent re-derivations, not
  new objects.
- **Cooper and Wallis (1972)** — the product theorem that turns a T-matrix
  quadruple and a Williamson-type quadruple into a Hadamard matrix.
- **Williamson; Baumert and Hall; Goethals and Seidel; Turyn** — the classical
  arrays and sequence families the assembler uses.
- **Cati and Pasechnik**, *A database of constructions of Hadamard matrices*,
  arXiv:2411.18897v2 (2025-08-30) — the status table against which the thirteen
  orders are reported. A database records what its authors knew, not what
  exists; the readings taken from it here are documentary, and every one is
  re-derived from the transcription shipped in `data/` at certificate run time.
