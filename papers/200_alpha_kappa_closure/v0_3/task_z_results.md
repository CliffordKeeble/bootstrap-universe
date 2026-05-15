# Paper 200 v0.3 Task Z - results (HALT)

Pre-registered in `task_brief.md`. Implementation in
`task_z_within_regime_null.py` (locked at pre-reg commit `61aff82`,
unmodified after run).

## Outcome: HALT on structural anchor failure

The brief's method step 4 (verify `S_5(31) = 0`) failed. Per the brief:

> If it doesn't reproduce, halt and report - something is wrong with the
> convention or the primes.

And per the brief's "What Mr Code should NOT decide" section:

> If during implementation any of these need to change, halt and flag to
> CinC. Do not adjust mid-run.

The locked script halted at the anchor check and did NOT proceed to compute
regime fractions or apply the bin decision. **No bin decision was
produced.** See `findings.md` for the diagnostic narrative and the flag
to CinC.

## What the script computed before halting

- Generated first 200 primes by sieve.
- Pre-flight cross-check (brief method step 1): `p_1 = 2`, `p_5 = 11`,
  `p_31 = 127` -> all match standard prime table. **OK.**
- Computed `chi_5(p_n)` under Definition A (per brief: ramified prime
  `p = 5` contributes 0; `chi_5(p) = +1` if `p % 5 in {1, 4}`,
  `-1` if `p % 5 in {2, 3}`).
- Computed running sum `S_5(n)` for `n = 1, ..., 200`.
- Anchor check: `S_5(31) = -6` (expected `0`). **HALT.**

## Anchor verification (hand-checked)

The first 31 primes, their residues mod 5, and `chi_5` under Definition A:

| n | p_n | p_n mod 5 | chi_5 | S_5(n) |
|---:|---:|---:|---:|---:|
| 1 | 2 | 2 | -1 | -1 |
| 2 | 3 | 3 | -1 | -2 |
| 3 | 5 | 0 |  0 | -2 |
| 4 | 7 | 2 | -1 | -3 |
| 5 | 11 | 1 | +1 | -2 |
| 6 | 13 | 3 | -1 | -3 |
| 7 | 17 | 2 | -1 | -4 |
| 8 | 19 | 4 | +1 | -3 |
| 9 | 23 | 3 | -1 | -4 |
| 10 | 29 | 4 | +1 | -3 |
| 11 | 31 | 1 | +1 | -2 |
| 12 | 37 | 2 | -1 | -3 |
| 13 | 41 | 1 | +1 | -2 |
| 14 | 43 | 3 | -1 | -3 |
| 15 | 47 | 2 | -1 | -4 |
| 16 | 53 | 3 | -1 | -5 |
| 17 | 59 | 4 | +1 | -4 |
| 18 | 61 | 1 | +1 | -3 |
| 19 | 67 | 2 | -1 | -4 |
| 20 | 71 | 1 | +1 | -3 |
| 21 | 73 | 3 | -1 | -4 |
| 22 | 79 | 4 | +1 | -3 |
| 23 | 83 | 3 | -1 | -4 |
| 24 | 89 | 4 | +1 | -3 |
| 25 | 97 | 2 | -1 | -4 |
| 26 | 101 | 1 | +1 | -3 |
| 27 | 103 | 3 | -1 | -4 |
| 28 | 107 | 2 | -1 | -5 |
| 29 | 109 | 4 | +1 | -4 |
| 30 | 113 | 3 | -1 | -5 |
| 31 | 127 | 2 | -1 | **-6** |

Counts in `[1, 31]`: 12 positive contributions, 18 negative, 1 zero -> sum = -6.

`S_5(31) = -6 != 0`. The anchor stated in the brief does not reproduce
under Definition A as defined.

## Definition B comparison (diagnostic)

Skipping the ramified prime `p = 5` (Definition B per brief: shifted
indexing; "out of scope for this task" but useful as a diagnostic) gives
`S_5(31) = -5`. Still not zero. See `findings.md` for the full Def B
column.

## Tie count in [1, 33] under Definition A

Across `n in [1, 33]`, `S_5(n)` takes the values:

```
-1, -2, -2, -3, -2, -3, -4, -3, -4, -3, -2, -3, -2, -3, -4, -5, -4, -3,
-4, -3, -4, -3, -4, -3, -4, -3, -4, -5, -4, -5, -6, -5, -6
```

`S_5(n) = 0` count: **0 of 33** (`f_tight = 0.0000` if we forced a count
diagnostically). Brief's Paper 199 sec.6.1 cross-check expected 14.
Discrepancy: **14**. This was the secondary flag the brief anticipated
("Paper 199 sec.6.1 may need its own audit"); it is now subsumed by the
primary anchor-failure flag.

## No bin decision

The brief's decision rule applies to `f_tight`, but `f_tight` is computed
*after* the anchor check passes. The anchor check failed; no bin decision
was reached.

## Files produced

- `findings.md` - integrated narrative, flags, candidate causes for CinC.
- `task_z_results.json` - machine-readable HALT record.
- (script output `task_z_results.json` was NOT written - the locked script
  halted before reaching the JSON-writing block; the JSON in this folder
  is hand-built from the diagnostic side-computation and clearly marked
  as a HALT record, not a results payload.)
