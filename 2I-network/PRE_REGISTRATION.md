# Pre-registration — 2I-network simulation

**Date**: 25 May 2026
**Author**: Mr Code
**Brief**: see [BRIEF.md](BRIEF.md) (CinC, 25 May 2026)
**Binding**: this document is committed to git **before any computation has been run**.

## Construction choice

The brief delegates one ambiguity to me: "4D embedding or careful combinatorial
bookkeeping. Use whichever is cleaner." I commit to:

**Geometric face-sharing growth in ℝ³ with vertex-merge consolidation.**

Concrete rules:

1. **Reference icosahedron.** 12 vertices at the standard positions
   `(0, ±1, ±φ), (±1, ±φ, 0), (±φ, 0, ±1)` divided by `√(1+φ²)` so the
   circumradius is 1. The 20 triangular faces are enumerated by finding all
   triples of vertices at the icosahedron's edge length.
2. **Placement.** Icosahedron 0 sits at the origin. To attach a new
   icosahedron N to face F of an existing icosahedron M, N is obtained by
   reflecting M across the plane of F. N's three F-vertices coincide
   with M's by construction; N's other nine vertices are placed in ℝ³.
3. **Vertex pool.** Every placed vertex is added to a global vertex pool.
   When a new icosahedron is placed, each of its 12 vertex positions is
   compared against the pool with tolerance ε = 1×10⁻⁶ (relative to unit
   circumradius). Any coincidence merges the vertices via union-find.
4. **Consolidation = vertex-pool sharing.** After placement, any face of
   N whose 3 vertices all coincide with 3 vertices of some other existing
   icosahedron's face is, by definition, the same face — those two
   icosahedra share that face. Edges in the network are added accordingly.
5. **Attachment order.** First-in-first-out across open faces. Within the
   pool of open faces, the oldest is picked. This is deterministic given
   the initial face ordering. RNG seed for any tiebreak: `seed = 20260525`.
6. **Edges in the network graph.** Two icosahedra share an edge in G
   iff they share **at least 3 vertices** (i.e., share a face). This is
   the face-sharing rule; it's what defines "neighbour" for Observable 1
   and 2.

Falsification of this construction choice is allowed by the brief's
stop-on-fail clause: if at N=120 the spectral gap is wildly off from any
reasonable interpretation of 168, I switch to **axis-sharing** (the
secondary construction) and report.

## Implementation choices that are NOT pre-registered

These are small choices, free to revise as needed, but reported in the
final write-up:

- BFS-order vs DFS-order for face selection — chosen FIFO above; will
  note if I revisit.
- Number of starting nodes for averaging `N_u(R)` in the dimension fit:
  `min(200, N)`.
- For the null model: 20 random regular graphs per N (as the brief
  specifies), generated with `networkx.random_regular_graph` and a
  reproducible seed.
- For the random regular null at non-integer average degree, round
  average degree to the nearest *even* integer (random regular graphs
  require `n*d` even).
- For sparse eigenvalue computation: `scipy.sparse.linalg.eigsh` with
  `sigma=0`, `which='LM'`, requesting 8 eigenvalues.

## What I will NOT do

- I will not look at Paper 203, Paper 189, or any paper that may bias the
  construction choice.
- I will not adjust the thresholds after seeing data.
- I will not silently revise the brief; if the construction smells wrong
  before running, I stop and flag (the brief's hygiene rule).
- I will not "fix" disagreements between the 2I-network's λ₁ and 168 by
  re-normalising or re-selecting nodes; the three normalisations listed
  in the brief are the only normalisations I report.
