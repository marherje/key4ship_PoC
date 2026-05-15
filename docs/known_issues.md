# Known Issues — MTC Tracking and Event Display

Diagnosed from `gaudi_jobs/muon_pipeline/debug_tracking_20.txt` (50-event run, single 5 GeV µ⁻ at xyz=(5, −82.5, 750) mm, dir (0,0,1)).

---

## Fixed

### U–V grouping straddle at 10 mm boundary (FIXED)

Adjacent U and V SciFi planes in the same MTC layer are separated by 2.35 mm in z.
The seeding code grouped them by `round(center.x() / 10.0)`, which placed ~2–3 layer pairs per station into different groups (U in bin N, V in bin N+1), preventing U×V crossings from forming for those layers.

**Fix applied:** changed divisor from `10.0` to `5.0` in `ACTSProtoTracker.cpp` (Source C grouping key). A 2.35 mm gap always fits within a 5 mm bin.

---

## Open Issues

### 1. Hough seeder does not handle curved tracks — wrong reconstructed position

**Symptom:** Reconstructed tracks appear far from SciFi hits in the event display. The dominant Hough peak frequently corresponds to a spurious crossing rather than the true muon position. Example from evt=0:
- True muon: (x=5, y=−82.5) mm
- Winning Hough seed: (x=12.0, y=−3.0) mm (6 votes)
- Correct-region peak: (x=6.5, y=−85.7) mm (4 votes, loses, flagged duplicate)

**Root cause:** The Hough accumulator bins in constant (x, y). This implicitly assumes the track is straight — i.e., has the same transverse position at every z. The MTC muon is **not** straight:

- B_y = 1.7 T in iron absorbers; R = p/(qB) = 5 GeV / (0.3 × 1.7) ≈ 9.8 m
- 45 iron slabs × 50 mm = 2250 mm of iron path
- Per-slab deflection angle ≈ 50/9800 ≈ 0.005 rad; over 45 slabs ≈ 0.23 rad total
- Cumulative transverse displacement over the full MTC (z = 800–4090 mm): hundreds of mm in x

Because the muon curves, the U×V crossings from different MTC layers land at different (x, y) positions following the curved trajectory. With HoughBinSize = 5 mm, the votes are spread across many bins. No single bin accumulates enough votes to dominate. Low-multiplicity spurious peaks (from secondary hits or random U–V pairings) win by default.

**Evidence:** chi2/nMeas ≈ 120 for most events (expected ≈ 1 for a good fit). The KF starts from the wrong seed position and fails to converge, since the 1D strip measurements do not give enough pull in both transverse coordinates simultaneously to recover from a large seed error.

**Required fix:** The straight-line Hough is the wrong tool for this geometry. Options:
- **Per-station local seeding:** fit a straight-line seed within each MTC station independently (curvature within a single 1m station is small), then link station seeds by extrapolating with a bending-angle hypothesis.
- **Hough in (x₀, p_x/p_z) space:** parameterise the track by its transverse intercept and slope; each measurement votes for a line in this space. Requires knowing the approximate momentum to constrain the Hough range.
- **Cellular automaton / track-following:** build track candidates by chaining compatible hits layer-by-layer, allowing a configurable bending between layers.

Until the seeder is fixed, tracks reconstructed in MTC-only mode (`muon_pipeline`) should not be trusted for position or momentum.

---

### 2. Event display draws a vertical straight line instead of the fitted trajectory

**Symptom:** The reconstructed track is rendered as a vertical line at constant (seed_x, seed_y) passing through every detector layer — it does not follow the actual track.

**Root cause — code:** In `event_display/event_display_eve.py`, `read_track_points()` builds the display points as:

```python
sx, sy = seed['seed_x'], seed['seed_y']
points = [(sx*mm2cm, sy*mm2cm, z*mm2cm, z) for z in layer_zs_mm]
```

This creates a column at fixed (x, y) for all z. `draw_tracks()` then connects these with straight line segments.

**Root cause — data:** `EDM4HEP2RNTuple.cpp` writes only `seed_x` and `seed_y` to the `ACTSTracks` RNTuple (extracted from the `AtIP` track state's `D0`/`Z0` fields). The Kalman filter output — fitted (x, y, dx/dz, dy/dz) at each measurement surface — is never written.

**Two separate errors:**

| Error | Effect |
|-------|--------|
| Track direction ignored | Even for a straight track, the particle enters at a slope (dx/dz, dy/dz); the display ignores this and shows a vertical column. |
| Magnetic curvature ignored | In the MTC iron (B_y = 1.7 T) the track bends in x by hundreds of mm over 3.3 m. A straight line in any direction deviates arbitrarily from the true path. |

**Required fix:**

1. **Export per-surface fitted track states from `ACTSProtoTracker`.** For each accepted ACTS track, write the fitted position (global x, y, z) at each measurement surface into the RNTuple (or into the `edm4hep::TrackState` collection).

2. **Update `event_display_eve.py` to use per-surface positions.** Connect the fitted positions with line segments rather than projecting a seed point through all layers.

As a minimal interim step, export the seed **direction vector** (fitted px/pz, py/pz at AtIP) so the display can draw a sloped straight line rather than a vertical one. This does not fix the curvature issue but at least gives the correct first-order direction.

---

## Interaction between the two issues

Issue 1 (wrong seed) feeds directly into Issue 2 (wrong display). Even if the display code were fixed to draw the actual fitted trajectory, the fitted trajectory would still be wrong because the Hough seeder gives the KF a bad starting position, causing poor convergence (chi2/nMeas ≫ 1).

The correct order of fixes:
1. Fix the seeder to find the correct track position → gives the KF a good seed.
2. Export KF output states to the RNTuple.
3. Use those states in the event display.
