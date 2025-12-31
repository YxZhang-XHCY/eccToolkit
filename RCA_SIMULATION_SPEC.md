# RCA Simulation Spec (CTMC Mode)

This document summarizes the current RCA simulation logic and nails down
the review-critical edge cases.

## Core Flow

1) Generate compartments (Poisson size), assign eccDNA instances by weight.
2) For each eccDNA instance, simulate RCA to get repeat count `R`.
3) Build an RCA molecule graph: trunk length `T = R * L`, plus branches.
4) Inject inter-molecule chimeras inside each compartment.
5) Debranch + nick/break to produce linear molecules.
6) Sample platform-specific reads (NGS/HiFi/ONT) from the linear pool.
7) Emit FASTQ + truth fields (segments, repeats, breakpoints, etc.).

## CTMC Repeat Model (Gillespie)

States: `E_free -> E_bound -> E_elong <-> E_pause -> E_off` with irreversible
`E_* -> E_inact`. Product accumulates only in `E_elong`.

Repeat count:
- `P_nt` accumulates as `P_nt += v_eff * dt` in elongation.
- `R_raw = P_nt / L`.
- Integer `R` uses unbiased rounding:
  `R = floor(R_raw) + Bernoulli(R_raw - floor(R_raw))`.

## Length Penalty (No Double Counting)

Length scale (only for `L > L_ref`):
- `scale = 1 + a * (1 - exp(-(L/L_ref - 1) / tau))`
- This is monotone (longer rings are penalized more), with saturation.

Penalty modes:
- `pause_off` (default):
  - `k_pause *= scale`, `k_off *= scale`.
  - `k_on` and `v_eff` do NOT change with length.
  - **Length factor inside `v_eff` is forced to 1 in this mode.**
- `full`:
  - `k_pause`, `k_off` increase by `scale`.
  - `k_on` and `v_eff` are multiplied by `1/scale`.
- `speed_only`:
  - only `v_eff` uses `1/scale`.
- `off_only`:
  - only `k_off` uses `scale`.

## Resource Competition

CTMC already produces a per-molecule yield. To avoid double penalties:
- `resource_mode=noise` (default) applies mild multiplicative noise only.
  - `R_eff = R_ctmc * lognormal(0, sigma)`, then **unbiased rounding**.
- `resource_mode=total_pool` is optional and stronger (shared pool).

## Rebinding Semantics

`rebinding_target` is explicit:
- `same_molecule` (implemented): an enzyme that goes off can rebind to the
  same template instance.
- `random_in_compartment` (not implemented yet): would allow rebinding to
  any template in the same compartment.

## HiFi Fragmentation (No Phase Artifact)

Long molecules are fragmented to 15-25 kb (mean 20 kb, std 1.2 kb).
To avoid fixed phase relative to molecule start:
- sample a random `start_offset` and begin cutting from there,
  then sequentially slice.
- if the prefix (0..start_offset) is long enough, it is emitted as a fragment.

## Repeat Truth for Chimeras

Truth now includes a per-source mapping:
- `repeat_count_by_source = {ecc_id: total_len / ecc_len}`
  (computed from TRUNK + CHIMERA segments; falls back to BRANCH if needed).
- `repeat_count_truth` is the **max integer** across sources for compatibility.
 - Background reads have `repeat_count_truth = 0` and empty `repeat_count_by_source`.

## Branch Resource Assumption

Branch generation does **not** reduce trunk length. This is a deliberate
model simplification: branch synthesis is treated as an additive process
that does not deplete the trunk yield.

## Recommended Sanity Checks (Optional)

1) CTMC identity:
   - With `k_pause = k_off = 0`, expect `P_nt ~= v * reaction_time` and
     `R_raw ~= v * T / L` (only discretization causes deviation).
2) Length penalty monotonicity:
   - In `pause_off`, mean `R_raw` decreases with `L` and saturates for large `L`.
3) HiFi read-length histogram:
   - Reads concentrated in 15-25 kb, mean ~20 kb, std ~1.2 kb.
   - Fragment breakpoints should look uniform (no phase artifact).
4) ctcR frequency:
   - For each `L` bin, check `P(repeat>=2)` vs. `read_len / L`.
5) Determinism:
   - With a fixed seed, `truth` should be bitwise reproducible.
