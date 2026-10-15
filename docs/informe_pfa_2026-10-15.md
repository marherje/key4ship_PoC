# SND-PFA — closing report (2026-10-15)

Particle flow for the SHiP SND in key4hep: what was built, how well it measures,
what was chosen and why, what did not work, and what is left. The step-by-step
record with every intermediate number is `docs/snd_pfa_progress.md`; the plan it
followed is `docs/snd_pfa_plan.md`. A shorter "how it is implemented" guide lives
in the Claude Docs page *SND-PFA in key4hep*.

All numbers are on the reference sample `pfa_ref_a` (2000 GENIE interactions, all
flavours, generated 2026-09-22 with the corrected flux frame, current geometry),
reconstructed with the sequential tracker and the job6 defaults of today.
Outputs on EOS: `/eos/experiment/ship/user/marquezh/snd_pfa/genie/`
(`pfa_pfa_ref_a_v2.root`, `ShipPFA_pfa_ref_a.root`, `final_*_ref_a.txt`).

---

## 1. Summary

- **What exists.** A Gaudi module `SND_pfa` with six reconstruction stages
  (calibrated hits → layer clusters → truth clusters → muon → CLUE3D 3D clusters →
  PFOs/event) in `gaudi_jobs/genie_pipeline/job6_pfa.py`, an RNTuple writer
  (`job7_pfa_rntuple.py` → `ShipPFA.root`), and the MOBO integration
  (`pfa_metrics`, experiment `snd_genie_pfa`). A baseline trial with its own GENIE
  generation (100 interactions, HTCondor) completes in 5 minutes; the PFA takes 12 s of it.
- **What works.** The hadronic energy sits on its physical ceiling: E_had/E_had,true
  = 0.952 overall against a ceiling of 0.933, NC −0.001 and ν_μ CC with the muon
  found +0.003 from it. The muon's own deposit is removed without bias (1.05 of the
  truth). When the muon is found and its momentum measured by curvature, E_ν/E_ν,true
  = 0.950 with σ68 0.205.
- **What does not meet the plan.** E_ν on all ν_μ CC: bias −0.193 against the
  criterion |bias| < 0.10. Muon identification: 0.647 of ν_μ CC (0.895 of those with a
  GlobalTrack, criterion 0.97). Both are limited by the same thing — **the muon is
  not tracked, or tracked too short** — not by the PFA's energy reconstruction.
  Muon fakes in NC meet the criterion (0.049 ≤ 0.05); in ν_e CC they do not (0.069).

## 2. Final metrics (job6 defaults, 2000 events)

| Metric | Value | Truth ceiling / reference | Plan criterion | |
|---|---|---|---|---|
| Muon ID efficiency, ν_μ CC | 0.647 | tracker gives the muon a GlobalTrack in 0.705 | — | |
| … among muons with a GlobalTrack | 0.895 | 1 | ≥ 0.97 | ✗ |
| Muon fakes, NC | 0.049 | — | ≤ 0.05 | ✓ |
| Muon fakes, ν_e CC | 0.069 | — | ≤ 0.05 | ✗ |
| Charge sign correct (p > 5 GeV) | 0.84–0.86 | — | > 0.80 | ✓ |
| p_μ bias / σ68, 10–20 GeV | −0.12 / 0.55 | — | — | |
| E_had / E_had,true (ν_μ CC + NC) | 0.952, σ68 0.381 | ceiling 0.933, σ68 0.345 | — | |
| E_had gap to ceiling: NC / ν_μ CC μ found / ν_μ CC no μ | −0.001 / +0.003 / +0.154 | 0 | — | |
| E_μ,subtracted / E_μ,deposited | 1.055, σ68 0.218 | 1 | — | |
| E_ν / E_ν,true, all ν_μ CC | 0.807, σ68 0.407 | — | \|bias\| < 0.10 | ✗ |
| … muon found | 0.896, σ68 0.333 | — | | |
| … p_μ by curvature (582 ev) | 0.950, σ68 0.205 | — | | |
| … p_μ only a range lower bound (126 ev) | 0.293 | — | | |
| … no muon found (357 ev) | 0.433 | — | | |
| Vertex found (cluster mode) | 0.94; Δz median −14.7 mm, σ68 65 mm; 10.4 % > 3 layers | — | — | |

The `pfa_metrics` dict of the baseline (`ShipPFA_pfa_ref_a.root`): eff_muon_id 0.646,
fake_muon_nc 0.049, fake_muon_nue 0.070, res_pmu 0.580, bias_ehad −0.048,
res_ehad 0.381, bias_enu −0.193, res_enu 0.407, frac_leak_mtc 0.753 (the 1-event
difference in the ν_μ CC count comes from classifying by neutrino flavour + CC
instead of by the primary lepton).

The E_had ceiling is the calibrated energy of every hit the primary muon does not
dominate; below it is longitudinal leakage (0.95 for vertices in the SiTarget, 0.64
at the back of the MTC — phase 1.5) and neutral energy that never deposits. No
algorithm recovers either.

## 3. Where E_ν loses its 19 %

| ν_μ CC events | Share | E_ν ratio | Why |
|---|---|---|---|
| muon found, p by curvature | 55 % | 0.950 | leakage and p_μ resolution only |
| muon found, p only bounded by its range in iron | 12 % | 0.293 | the fit never measured q/p; the iron crossed is a lower bound |
| no muon found | 34 % | 0.433 | E_μ is missing entirely |

The 45 % of events in the last two rows carry the bias. Phase 3's loss analysis
(`muon_loss.py`) puts the missing muons at: **32 % crossing fewer than 8 MTC layers**
(geometry: vertices deep in the MTC, low p, large angle — only 76 of those 355 have
any track), **3.6 % with enough layers but no track**, ~1 % lost in the finder. The
iron-bound row is also a tracking statement: those muons have a track whose
curvature was not measured. **The next lever on E_ν is the MTC tracking**, which is
the work that follows this report.

## 4. Defaults chosen

| What | Default | Why (phase) |
|---|---|---|
| Hit threshold | 0.5 MIP | 1.2 |
| 1 MIP | Landau MPV of a 10 GeV muon per plane type (`MIPEnergyOverride`) | mean dE/dx put the MIP peak at 0.73–0.88 (1.4) |
| Sampling weights | from the compact's material, per layer | follows MOBO variants (1.1) |
| `GlobalScaleECAL` / `GlobalScaleMTC` | 1.51601 / 1.31372 | e⁻ 1–20 GeV from SiTarget; π⁺ from inside MTC40 (1.4) |
| e/h | silicon 1.07546, MTC 1.0 | MTC scale already fitted on pions (4.1) |
| Layer clustering | connectivity + split at local maxima > 2 MIP | 2.1–2.3: purity median 1.00 |
| Muon momentum | curvature + ionisation loss; range only as lower bound | the range measured tracking length, not p (3.2) |
| Muon ID | 8 MTC layers (6 when the corridor reaches the last layer), 60 % of spanned, mean ≤ 2 MIP; SciFi half-width 2 mm with the stereo term | 3.4 / efficiency study, 2026-10-14 |
| Muon subtraction | `trajectory`: mean deposit per plane vs p_μ (guns 3/10/30/100 GeV), remainder from the event total | 4.2b |
| E_had | sum over layer clusters (CLUE3D optional: `PFA_HADRONIC_CLUSTERS=HadronicClusters`) | CLUE3D loses 6–7 % in outliers (4.3) |
| Vertex | first layer with a > 2 MIP cluster followed by activity | found in 94 % vs 14 % with tracks (4.2) |
| CLUE3D | dc 10/15/40 mm, ρc 0.02/0.04/0.2 GeV (SiTarget/SiPad/MTC), transverse distance, ±2/±3 layer windows | 4.3 |

## 5. What did not work, and why it is worth remembering

- **Muon momentum from range** (decided 2026-09-22, reverted 2026-09-27). The MTC has
  2.35 m of iron, worth 2.7 GeV; SHiP muons are 2–100+ GeV and none stop. The "range"
  was the reconstructed track length, biasing p by −55 to −85 %.
- **One `ELossScale` for the muon.** It overcorrects below 5 GeV and undercorrects above
  10 (radiative losses); left at the physical 1.0.
- **One-MIP corridor subtraction (4.1).** One MIP is the Landau peak; a muon deposits a
  mean of 1.35–2.9 of those per plane, and the corridor also refused the > 3 MIP
  clusters that hold that tail. E_had came out 16–20 % high in ν_μ CC.
- **Trajectory subtraction capped per plane (first 4.2b).** min(mean, cluster) takes
  only 56–83 % of a Landau-distributed deposit, whatever the mean. Making the mean
  depend on p_μ changed nothing until the remainder went to the event total.
- **Taking isolated muon clusters whole.** 2D connectivity in the 50 mm scintillator
  merges the shower and the muon into one cluster with no neighbours: 4.5 GeV taken
  for 2.75 deposited.
- **CLUE3D for E_had.** Purity 0.72 and a useful 3D picture, but its outliers only ever
  remove energy; E_ν stays better with the layer sum.
- **The 30 mm cluster-axis-to-vertex link.** Straight fits through pieces of hadron
  showers in iron point back badly (47 % within 100 mm); it threw away 80 % of the 3D
  energy. Disabled until a sample with backgrounds gives it something to reject.
- **The plan's single e/h = 1.4.** Predated the calibration; on the MTC it would have
  pushed E_had 40 % high.

Bugs found on the way, all fixed: the MTC cell decoder sent every scintillator tile
to channel 0 (0.4); the q/p sign of the persisted track state (0.3); objects keyed by
(collectionID, index) collided because every collection made in the job has
collectionID 0xffffffff until written (4.2b — it had affected 4.1–4.3); the muon
corridor compared stereo U/V clusters with the track x without ±y·tanα (3.4); the
shared-sample GENIE chain of the MOBO (`snd_genie`) replays the August sample with the
misplaced flux (5.3) — the PFA experiment is built on `snd_genie_gen`, which generates
GENIE on each variant and does not use it.

## 6. Pending

1. **MTC tracking** (agreed with the user as the next step): muons with 4–7 MTC
   layers that get no track or no curvature measurement. Improve `ACTSSequentialTracker`
   itself, no extra post-algorithms.
2. ν_e CC muon fakes (0.069) — pions crossing the MTC as MIPs; beyond the cuts, needs
   more information (e.g. the vertex topology).
3. A z-dependent leakage correction for E_had (phase 1.5 measured the ceiling vs vertex
   z); the plan deliberately did not fit it.
4. p_μ resolution above 5 GeV (σ68 0.50–0.55) is the MTC curvature measurement.
5. Phase 6 (EM showers, ν_e CC), only if asked.
6. `snd_genie` (shared sample) still replays the August GENIE file with the misplaced
   flux. Not changed (outside the PFA plan); `snd_genie_gen` and `snd_genie_pfa` are
   unaffected.
7. The calibration constants were fitted on the baseline geometry; for far MOBO
   variants a recalibration step (or an `E/E_true` monitor) would make `bias_*` metrics
   trustworthy. `res_enu` (a width) is the safer objective meanwhile.
