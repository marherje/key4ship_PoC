# SND-PFA — registro de progreso

Guión maestro: `docs/snd_pfa_plan.md`. Este fichero lo edita la sesión que implementa: una
fila por paso, con fecha y los números que pide el bloque **Anota** del guión. Estado:
`☐` pendiente · `◐` en curso · `☑` cerrado · `✗` bloqueado (explicar en la fila).

## Estado por paso

| Paso | Estado | Fecha | Anotaciones (números, rutas, defaults elegidos) |
|---|---|---|---|
| 0.1 rama y registro | ☑ | 2026-09-22 | rama `feature/snd-pfa` desde `feature/genie-events`; enlace añadido al final de `docs/gaudi_pipeline.md` |
| 0.2 `MCParticles` en overlay | ☑ | 2026-09-22 | 2000 ev: 0 sin MCParticles, 0 índices fuera de rango, PDG(MCParticle)==ContribPDGs 3729389/3729389 (100 %). **El neutrino SÍ está**: generatorStatus=4 (beam) en los 2000, E de 1.00 a 231 GeV (mediana 5.65) → `E_ν,true` se lee directamente del neutrino. 1481/2000 (74 %) con leptón cargado primario (CC-like). El vértice de interacción es el de las hijas, no el del neutrino (z=0). |
| 0.3 q/p persistido (`AtLastHit`/`AtFirstHit`) | ☑ | 2026-09-22 | Escrito y utilizable: 687 trazas MTC dominadas por el μ primario, 0 sin estado, relVar~0.026. **Signo corregido** (ver desviaciones): carga correcta 0.90/0.90/0.97/0.88 por bins (antes 0.10). dp/p vs endpoint +0.42/+0.16/+0.08/−0.02 → criterio p>10 GeV **pasa**. σ68 vs endpoint 0.39–0.98: la resolución de curvatura del MTC está dominada por dispersión múltiple. **Decisión del usuario (2026-09-22): el rango pasa a ser la medida principal de momento en la Fase 3**; la curvatura queda para la carga y para los muones que salen. Criterio final: signo > 0.80 y |dp/p vs endpoint| < 0.15 (p>10 GeV) → **OK**. |
| 0.4 esqueleto `SND_pfa` + `SNDCellDecoder` | ☑ | 2026-09-22 | módulo `SND_pfa` compila e instala; `job6_pfa.py` corre. `check_decoder.py`: **59439 celdas comparadas de las 4 colecciones, 0 discrepancias**. El cross-check cazó un bug real: el descriptor del MTC declara `strip`, `x` e `y` en todos los planos, así que decidir por "existe el campo" mandaba cada tile del scint al canal 0; ahora se decide por la vista. |
| 1.1 pesos de muestreo | ☑ | 2026-09-23 | 388 planos (240 SiTarget + 22 SiPad + 126 MTC), 0 materiales desconocidos. **Peso por capa** (ver desviaciones): SiTarget pasivo 7.656 / activo 0.2326 → **w = 33.910** (X y Y); SiPad 22.50 / 0.2520 → **w = 92.252**; MTC 65.26 / 4.582 → **w = 14.963** (U, V y scint). Capa 0 distinta por detector (90.27, 15.24): tiene menos material delante. **Periodicidad verificada**: un único w por detector. **Variante** (`variants/pfa_wtest.xml`, SiTarget W 3.5→7.0, SiPad 10→5.0): SiTarget ×2.00, SiPad ×0.52, MTC sin cambio. |
| 1.2 `SNDCaloHitBuilder` completo | ☑ | 2026-09-23 | 50 eventos GENIE: 532.5 hits/evento → **~250 escritos** (umbral 0.5 MIP), 0 sin decodificar. `check_calohits.py` **OK**: hits = links, 0 discrepancias de cellID, 0 bajo umbral, **0/5910** hits del MTC con la `y` del sim hit. E_reco/E_dep idéntico dentro de cada detector (33.910 / 92.252 / 14.963). **Asimetría entre vistas ortogonales: X/Y 1.5 %, U/V 0.6 %** — las dos vistas miden el mismo paso, como debe. `energyError` transporta el valor en MIP. |
| 1.3 guns de calibración (Condor) | ☑ | 2026-09-23 | Lanzados 11 jobs (clusters 20712106–20712118) con `simulation/run_script/launch_calib_guns.sh` (nuevo): μ⁻ 10 GeV; e⁻ 1/2/5/10/20; π⁺ 2/5/10/20/40, 500 eventos cada uno, vértice (1,1,−1640) mm = 50 mm dentro del SiTarget. Salidas en **EOS**: `/eos/experiment/ship/user/marquezh/snd_pfa/calib/output_calib_<part>_<E>GeV.edm4hep.root`. Todas terminadas (367 MB en EOS). Reconstruidas con `gaudi_jobs/calib_pipeline/` (nuevo: `job1_overlay.py`, `job6_pfa.py`, `run_calib.sh`), salidas en `.../calib/reco/`. **Añadidas 8 guns más con vértice z=+200 mm (dentro de MTC40)**, clusters 20712931–38: sin ellas no hay forma de calibrar el MTC (ver 1.4). |
| 1.4 calibración | ☑ | 2026-09-26 | **MIP**: la primera pasada dio modas 0.73–0.88, no ~1, porque el peso usa el dE/dx **medio** mínimo mientras una capa fina registra el **valor más probable** (Landau MPV, los δ escapan); el efecto disminuye con el espesor (0.73 en 300 µm de Si, 0.88 en 20 mm de centelleador), que es la firma esperada. Fijado con `MIPEnergyOverride` medido del μ de 10 GeV → segunda pasada: **modas 0.98–1.03** en los seis planos ✓. Validación cruzada: el MIP del scint del MTC sale **3.53 MeV** y el informe de tracking midió ~3.5 MeV/MIP por otro camino. **EM**: E_vis/E = 0.660 lineal, residuo máximo **1.3 %** (criterio <5 %) ✓; σ/E ≈ **21 %/√E** (0.215 a 1 GeV → 0.046 a 20 GeV), dentro del 15–25 % esperado. → **GlobalScaleECAL = 1.51601**. **Hadrones**: E_vis/E crece de 0.521 (2 GeV) a 0.620 (40 GeV), residuo **15.5 %** — **no cumple el criterio <10 %**, pero es no-linealidad intrínseca de un calorímetro no compensado, no un fallo del código: un solo factor no puede corregirla (ver preguntas abiertas). e/π medio = 1.070. **Fuga al MTC**: media 0.1 %→2.6 % y fracción de eventos que lo alcanzan 0.8 %→38 % entre 2 y 40 GeV. **MTC cerrado con las guns de vértice z=+200 mm** (dentro de MTC40, donde de hecho interacciona la mayoría de los neutrinos): `GlobalScaleMTC = 1.31372`, y con ella E_vis/E = 0.949 / 0.955 / 0.974 / 0.995 / 1.005 para π⁺ de 2/5/10/20/40 GeV, **residuo 5.1 %** — el MTC es **más lineal que el silicio** (15.1 %) para hadrones, lo que conviene tener en cuenta al construir E_had. Con la escala EM aplicada la respuesta a electrones es exacta (E_vis = 0.99/1.97/5.01/10.00/20.00 GeV). |
| 1.5 cota E_had en GENIE | ☑ | 2026-09-26 | 400 ev. **El criterio del guión (corr>0.9, sesgo<10 %) no se cumple, y la causa es física, no del código: fuga longitudinal que depende de dónde interaccione el neutrino.** Ratio E_vis,nomu/E_had,sum por z del vértice: **0.949** (−2000..−500, SiTarget) → 0.888 → 0.790 → 0.728 → **0.638** (2500..3400, final de MTC60). Una interacción en el blanco tiene 1.8 m de detector delante; una en la última estación pierde el 36 %. De ahí el sesgo global −19.4 % y la correlación 0.781 (la fuga varía evento a evento). Desglose de `E_had,sum`: 41.9 % π±, 23.1 % kaones/charm/γ/Λ, 20.8 % π⁰, 8.1 % p, 6.0 % n — parte nunca se deposita (neutrones, K⁰ PDG 311 que Geant4 no propaga, ν de decaimientos). La definición alternativa `E_ν − E_leptón` sale **peor** (corr 0.590, sesgo −32.9 %) porque añade la energía del blanco nuclear y los fragmentos. **Cota para la Fase 4: 0.95 si el vértice está en el SiTarget, 0.64 al final del MTC**; la corrección tiene que depender de z, no ser un factor. |
| 2.1 `SNDLayerClusterer` | ☑ | 2026-09-26 | 2000 ev de `pfa_ref_a`: clusters/evento 237.9 (SiTarget) / 15.0 (SiPad) / 101.0 (MTC SciFi) / 20.6 (MTC scint). Media de hits por cluster 1.6–2.1, con 19–32 % de más de un canal. Anchura media 0.96 mm (SiTarget, strips de 75 µm), 2.43 mm (SiPad, pads de 5.53), 0.81 mm (SciFi, strips de 1 mm), 35.6 mm (scint, celdas de 50 mm) — coherente con cada granularidad. Nota: la **mediana** de clusters/evento es 0 en SiTarget y SiPad porque, con el flujo bien colocado, más de la mitad de los neutrinos interacciona aguas abajo del silicio. |
| 2.2 `SNDTruthClusterer` + `SNDTruthLink` | ☑ | 2026-09-26 | 15105 clusters de verdad en 1980 de 2000 eventos (20 sin truth utilizable: no dejan hits sobre umbral). `SNDTruthLink.h` es el único sitio que une `*_Contributions` con `*ContribParticleIDs` + offsets. El índice de la MCParticle va en `shapeParameters[0]`, su PDG en `[2]`. Con IOSvc el event store **no** es un `PodioDataSvc`, así que el algoritmo reabre el fichero con su propio `podio::ROOTReader`, igual que `EDM4HEP2RNTuple`. |
| 2.3 validación capa (pureza, aislamiento del muón) | ☑ | 2026-09-26 | 300 ev. **Pureza mediana 1.00 en los cuatro** (criterio ≥0.90 en SiPad y scint) ✓; fracción con pureza exactamente 1: 0.89 SiTarget / 0.91 SiPad / 0.87 SciFi / **0.76 scint** — el scint es el más mezclado, como corresponde a celdas de 50 mm. Clusters por partícula y capa: 2.39 / 1.00 / 2.00 / 1.00 — los ~2 de SiTarget y SciFi son las **dos vistas** (X/Y, U/V) de la misma partícula, no sobre-división. **Aislamiento del muón: 87.2 % dentro de la mancha hadrónica, 93.8 % fuera** (cluster propio de ≤3 MIP, ≤3 canales, pureza ≥0.8). Es una buena noticia para la Fase 3. |
| 3.1 corredor traza–calo | ☑ | 2026-09-27 | `SNDMuonFinder` construye el corredor por la polilínea de los estados `AtOther` con extrapolación recta (`ExtrapolateMm` 400). Perfil de hierro derivado del mismo barrido de material que los pesos: **42 capas MTC, 2352.75 mm de Fe** (= 42 × 56 mm, exacto). Un cluster sólo puede pertenecer a un corredor; las trazas se ordenan por nº de hits. 842 candidatos en 2000 eventos (0.42/evento). |
| 3.2 ID + momento | ☑ | 2026-09-27 | ID: `type&4`, MinMTCLayers 8, MinFracLayers 0.6, MaxMeanMIP 2.0. **Momento = curvatura, con el rango sólo como cota inferior** (decisión del usuario del 2026-09-27, revirtiendo la del 22-sep; motivo en desviaciones). Añadida la **corrección de pérdida de energía** con la misma tabla de material de la Fase 1 (el fit mide q/p al final del corredor, tras 2–4 GeV de ionización). Evolución del sesgo por bins de p_true: rango **−0.55…−0.85** → curvatura **−0.08…−0.38** → curvatura + dE/dx **+0.31 / +0.19 / −0.11 / −0.14 / −0.12**. σ68 0.16–0.50. Sobrecorrige por debajo de 5 GeV y subcorrige por encima de 10 (el muón radia, lo que la estimación MIP no modela); un solo `ELossScale` no arregla ambos extremos, así que queda en 1.0 (la estimación física). |
| 3.3 validación muón | ☑ | 2026-09-27 | 500 ev. **Eficiencia** ν_μ CC: 0.608 (todos), 0.781 (μ con ≥20 hits). Referencia del tracker: el μ tuvo GlobalTrack en **0.703**; de ésos el finder identifica **0.849** → cuesta **15.1 puntos**, y el criterio del guión era ≤3. **Fakes**: NC **0.111**, ν_e CC **0.096** (criterio ≤0.05); los corredores falsos los construyen π± (17 de 23), protones (3). **Carga**: 0.84–1.00 correcta por bins ✓. **Momento**: sesgo −0.26 a −0.85 — inaceptable, causa identificada (ver preguntas). **Barrido 3×3 hecho** (`Optimization/Analysis/pfa/scan_muon_cuts.sh`, 500 ev por punto, salidas en `.../scratch/muon_scan/`): MinMTCLayers × CorridorHalfWidthMTCSciFi. Eficiencia relativa / fakes NC: L6W25 **0.905/0.137**, L6W15 0.899/0.128, L6W10 0.879/0.111, L8W25 0.859/0.111, L8W15 0.849/0.111, **L8W10 0.824/0.085** (elegido), L10W15 0.754/0.060, L10W10 0.744/**0.051**. **Ningún punto cumple los dos criterios del guión** (≥0.97 relativa y ≤0.05 fakes): son un compromiso continuo. Lo que lo limita no son los cortes sino cuántas capas del MTC reconstruye el tracking (mediana 8 capas a 2 GeV, 34 por encima de 20); mejorar ambos a la vez pasa por el tracking del MTC, no por el finder. Los fakes son π± (17 de 23) que cruzan el MTC con perfil MIP: físicamente indistinguibles de un muón sin más información. |
| 4.1 `SNDPFOBuilder` (E_had, E_ν, vértice) | ☑ | 2026-10-13 | `gaudi_source/pfa/SNDPFOBuilder.cpp`, añadido a `job6_pfa.py`. Cadena de 2000 ev en 30 s: 750 con muón, vértice por clusters en 1872, por trazas en 279. 184 capas en la secuencia del haz (120+22+42). Los `pfa_*` van en la colección **`PFASummary`** (un `ParticleID` por evento) con los nombres en el metadato `PFASummaryNames`, no como parámetros de frame (ver desviaciones). Se guardan **los dos vértices** (`pfa_vtx_clu_*`, `pfa_vtx_trk_*`) además del elegido, para comparar sin re-ejecutar. **E/H medido**: `EoverHECAL = 1.075` (pendiente e 1.000 / π 0.930 con las escalas actuales, `calib.py`), `EoverHMTC = 1.0`. |
| 4.2 validación E_had/E_ν | ✗ | 2026-10-13 | `Optimization/Analysis/pfa/enu_perf.py`, 2000 ev (salida en `.../snd_pfa/genie/enu_perf_ref_a.txt`). **bias_ehad = +0.023, res_ehad = 0.454** (ν_μ CC + NC; cota medida en los mismos eventos 0.933, σ68 0.345). **bias_enu = −0.172, res_enu = 0.426** → **no cumple** |sesgo| < 10 %. **VertexMode = clusters**: lo encuentra en el 94 %, Δz mediana −14.7 mm, σ68 65 mm, 10.4 % a más de 3 capas; `track` sólo en el 14 % y con el 41 % a más de 3 capas. Tablas en la sección "Fase 4.2" más abajo. |
| 4.3 clusters 3D (decisión del usuario) | ☑ | 2026-10-13 | **decisión: CLUE3D en lugar del cono** (usuario, 2026-10-13). `SNDHadronicClusterer.cpp` + `SNDClue3D.h` (algoritmo) + `SNDMuonChoice.h` (regla común de qué candidato es el muón). 2000 ev: 201 puntos de capa → **17.6 clusters 3D/evento** (mediana 4), 15 % de los puntos sin agrupar (5.8 % de la energía). **Pureza ponderada en energía 0.64**; el 21 % de la energía está en clusters de pureza > 0.9. **Para E_had no mejora la suma por capas** (tabla en "Fase 4.3"): E_ν 0.792 frente a 0.828. Queda como **opcional** (`PFA_HADRONIC_CLUSTERS=HadronicClusters`); el default sigue siendo el de 4.1. La colección `HadronicClusters` se escribe siempre. |
| 4.2b resta del muón por trayectoria (opción (a), usuario) | ☑ | 2026-10-14 | `MuonSubtraction = trajectory` en `SNDPFOBuilder` (default). Exceso de E_had en ν_μ CC con muón **+0.199 → +0.099**, σ68 0.642 → 0.432; NC intacto. Corrige además un bug de claves (collectionID sin asignar en memoria) que afectaba a 4.1–4.3. **Closed (English from here):** `MuonMeanMIP` tabulated vs p_μ from 3/10/30/100 GeV guns; the per-plane cap was the real bias (takes 56–83 % of a Landau), so the remainder is now taken from the event total → excess **+0.002**, σ68 0.372, E_mu,sub/E_mu,true 1.05. E_ν criterion still fails (−0.190): no-muon events and leakage. See "4.2b, closed". |
| 3.4 corridor stereo fix (user, 2026-10-14) | ☑ | 2026-10-14 | `SNDMuonFinder` compared U/V SciFi clusters with the track x without ±y·tanα (`MTCStereoAngleDeg`, 5°). ν_μ CC efficiency 0.615 → **0.636**; identified among muons with a GlobalTrack 0.850 → 0.877; fakes NC 0.066 → 0.076, ν_e CC 0.083 → 0.093. See "Muon efficiency study". |
| 5.1 `PFASummaryWriter` | ☑ | 2026-10-15 | `gaudi_source/pfa/PFASummaryWriter.cpp` + `gaudi_jobs/genie_pipeline/job7_pfa_rntuple.py`, appended (with job6) to `genie_pipeline_seq.sh`. Reads the PFA file with its own podio reader (truth in frame parameters, column names in run metadata). `ShipPFA_pfa_ref_a.root`: **EventPFA 2000 rows** ✓, PFOs 2765, LayerClusters 749049 (`PFA_WRITE_CLUSTERS=1`). Cross-check with the 6/2 defaults: muon efficiency 689/1066 = 0.646, NC fakes 0.049, ν_e CC fakes 0.070 — muon_perf.py gives 689/1065: it classifies by the primary lepton, EventPFA by the neutrino flavour + `is_cc`, which differ on one event. Extra columns beyond the plan: `E_vis_nomu` (truth ceiling of 1.5), `E_mu_dep_true`, `E_mu_sub`, `E_had_layers`, `vtx_source`. |
| 5.2 `pfa_metrics` + payload + test | ☑ | 2026-10-15 | `mobo.ship.metrics.pfa_metrics` (+ `pfa_metrics_from_columns`, testable without ROOT); `Optimization/Analysis/job{6_pfa,7_pfa_rntuple}.py` generated from the pipeline ones ($SND_COMPACT, $SND_OUTPUT_LEVEL, bare trial filenames); evaluator option `pfa_metrics` appends them after job5 and merges the dict (missing EventPFA → `StepFailed`, not a zero); `CONSUMED_BY` drops tracks_seq after job6 and pfa_seq after job7. New `conf/evaluator/snd_genie_pfa.yaml` + `conf/experiment/snd_genie_pfa.yaml`, copies of **`snd_genie_gen`** (GENIE generated per variant; the plan's `snd_genie` predates it) with `pfa_metrics: true` and objectives `eff_muon_id` max, `res_enu` min, `cost_proxy` min; `snd_genie*.yaml` untouched. **`tests/ship/test_pfa_metrics.py`: 3 passed** (definitions by hand, missing file, a 20-row EventPFA written with RNTupleWriter). The 22 old MOBO failures are unchanged. **Baseline dict** (`ShipPFA_pfa_ref_a.root`, 2000 events, job6 defaults): eff_muon_id 0.646, fake_muon_nc 0.049, fake_muon_nue 0.070, res_pmu 0.580, bias_ehad −0.048, res_ehad 0.381, bias_enu −0.193, res_enu 0.407, frac_leak_mtc 0.753, nevents_pfa 2000 — same as muon_perf.py / enu_perf.py. |
| 5.3 trial local de humo | ☑ | 2026-10-15 | **First version, shared sample** (`executor=local`, 200 events): the chain works end to end in 111 s wall, the PFA adds 6.9 s (job6 3.9, job7 2.9), `metrics.json` has every PFA key, no ROOT file left behind. But it replayed `run_sim_genie.py`'s default `genie_allflavours_0_2000.hepmc` — the August sample with the misplaced flux — and gave bias_ehad +0.30, res_ehad 0.96. On `genie_pfa_ref_a.hepmc` the same trial gave eff_muon_id 0.670, fake_muon_nc 0.061, res_pmu 0.552, bias_ehad −0.013, res_ehad 0.378, bias_enu −0.167, res_enu 0.406 — the 2000-event reference within statistics. **The user pointed out the MOBO generates GENIE per variant**: `snd_genie_pfa` was rebuilt on `snd_genie_gen` (the temporary `hepmc` evaluator option was removed again) and its smoke trial ran on HTCondor (`executor=htcondor_genie`, 6 GB slot; driver local under `nice -n 19 timeout 6h`, run_dir on AFS, trials on EOS `.../snd_pfa/scratch/mobo_smoke_pfa_gen/`): **COMPLETED, 298 s wall** for 100 generated interactions (GENIE + ddsim 68 s, job6 6.7 s, job7 4.9 s). eff_muon_id 0.622, fake_muon_nc 0.034, fake_muon_nue 0.0, res_pmu 0.501, bias_ehad −0.073, res_ehad 0.487, bias_enu −0.315, res_enu 0.334, frac_leak_mtc 0.716, eff_track 0.49, n_int_per_1e18pot 85901 — consistent with the reference within the statistics of ~45 ν_μ CC (a campaign wants more than 100 interactions per trial for these objectives). |
| 6 EM (si se pide) | ☐ | | |
| Cierre: `docs/informe_pfa_<fecha>.md` | ☑ | 2026-10-15 | `docs/informe_pfa_2026-10-15.md`. End-to-end on `pfa_ref_a` with today's defaults (`final_*_ref_a.txt` on EOS): check_truth OK, check_qop OK, muon_perf 0.647 / NC 0.049 / ν_e 0.069, enu_perf bias_ehad −0.048, bias_enu −0.193 (FAIL, as known). |

## Fase 4.2 — tablas (2026-10-13, `pfa_ref_a`, 2000 ev)

E_had,reco / E_had,true (Σ E_cin de los hadrones primarios; ν_μ CC + NC, 1548 ev). "Cota" =
energía calibrada de los hits que no domina el muón primario, con el mismo e/h, en los mismos
eventos (la del paso 1.5 se midió antes de `GlobalScaleMTC` y sin e/h, por eso cambia).

| E_had,true [GeV] | N | ratio | σ68/E | cota | diferencia |
|---|---|---|---|---|---|
| 0–2 | 391 | 1.495 | 1.595 | 0.823 | +0.672 |
| 2–5 | 227 | 1.030 | 0.510 | 0.881 | +0.149 |
| 5–10 | 200 | 0.971 | 0.351 | 0.878 | +0.093 |
| 10–20 | 299 | 0.966 | 0.238 | 0.923 | +0.043 |
| >20 | 431 | 0.992 | 0.189 | 0.974 | +0.018 |

| z_vtx,true [mm] | N | ratio | σ68/E | cota | diferencia | cota 1.5 |
|---|---|---|---|---|---|---|
| −2000..−500 | 233 | 1.105 | 0.386 | 0.962 | +0.143 | 0.949 |
| −500..500 | 416 | 1.075 | 0.521 | 0.974 | +0.102 | 0.888 |
| 500..1500 | 323 | 1.037 | 0.391 | 0.951 | +0.086 | 0.790 |
| 1500..2500 | 308 | 1.009 | 0.315 | 0.929 | +0.081 | 0.728 |
| 2500..3400 | 258 | 0.814 | 0.502 | 0.727 | +0.088 | 0.638 |

| clase | N | ratio | σ68/E | cota | diferencia |
|---|---|---|---|---|---|
| NC | 485 | 0.914 | 0.332 | 0.915 | **−0.001** |
| ν_μ CC, muón encontrado | 676 | 1.130 | 0.530 | 0.969 | **+0.161** |
| ν_μ CC, muón no encontrado | 387 | 0.990 | 0.484 | 0.832 | **+0.158** |

E_ν,reco / E_ν,true en ν_μ CC (1065 ev):

| E_ν,true [GeV] | N | ratio | σ68/E | muón encontrado |
|---|---|---|---|---|
| 0–2 | 45 | 0.812 | 0.621 | 0.22 |
| 2–5 | 83 | 1.024 | 0.572 | 0.57 |
| 5–10 | 104 | 0.654 | 0.445 | 0.41 |
| 10–20 | 180 | 0.813 | 0.424 | 0.59 |
| >20 | 653 | 0.831 | 0.409 | 0.72 |

| `pfa_p_mu_source` | N | ratio | σ68/E |
|---|---|---|---|
| curvatura | 576 | **0.992** | 0.231 |
| cota de hierro | 102 | 0.291 | 0.263 |
| sin muón | 387 | 0.432 | 0.378 |

Vértice (reco − MC, mm):

| modo | encontrado | med Δx | med Δy | med Δz | σ68 Δx | σ68 Δy | σ68 Δz | >3 capas |
|---|---|---|---|---|---|---|---|---|
| clusters | 0.94 | −0.1 | 0.0 | −14.7 | 48.8 | 77.5 | 65.2 | 0.104 |
| track | 0.14 | 0.0 | 0.2 | 1.2 | 17.9 | 17.7 | 99.1 | 0.409 |

Fuga (E_had,MTC/E_had > 0.3): 0.856 / 0.836 / 0.780 / 0.736 / 0.763 por bins de E_had — no es
fuga sino que la mayoría de los neutrinos interacciona dentro del MTC.

**Diagnóstico.** El sesgo de E_ν (−17 %) tiene dos causas, y ninguna es la suma de energía:
1. **El 36 % de los ν_μ CC no tiene muón** (387/1065) y entonces E_ν = E_had: ratio 0.43. Es la
   eficiencia del `SNDMuonFinder` (fase 3.3), no algo del PFO builder. Con el muón medido por
   curvatura E_ν sale **sin sesgo (0.992)**, y con muón encontrado el sesgo es **−5.1 %**, dentro
   del criterio.
2. **La cota de hierro como momento** (102 ev) da 0.29: cuando el ajuste no midió q/p, el mínimo
   compatible con el rango se queda muy corto, por lo explicado en 3.2.
Además, en ν_μ CC E_had queda **+16 % por encima de la cota** tanto si el muón se encuentra como
si no, mientras que en NC coincide con ella (−0.001): parte del depósito del muón acaba en E_had.
Un muón que cruza todo el detector deja, calibrado como si fuera cascada, ~1 GeV en el SiTarget
y ~3 GeV en el MTC. Con muón encontrado se escapa lo que el corredor no recoge: el semiancho de
10 mm en el SciFi (elegido en 3.3 para bajar los fakes), los clusters de > 3 MIP que el corredor
rechaza por parecer cascada, y el tramo antes de que empiece la traza.

## Fase 4.3 — CLUE3D frente a la suma por capas (2026-10-13, `pfa_ref_a`, 2000 ev)

Mismos eventos, misma escala; la diferencia es contra la cota de cada clase (E_had) o contra 1 (E_ν).
CLUE3D con dc = 10/15/40 mm y ρc = 0.02/0.04/0.2 GeV (SiTarget/SiPad/MTC), sin enlace al vértice.

| clase | 4.1 ratio | σ68 | diferencia | 4.3 ratio | σ68 | diferencia |
|---|---|---|---|---|---|---|
| E_had NC | 0.914 | 0.332 | −0.001 | 0.855 | 0.362 | −0.060 |
| E_had ν_μ CC, muón encontrado | 1.130 | 0.530 | +0.161 | 1.067 | 0.517 | +0.098 |
| E_had ν_μ CC, sin muón | 0.990 | 0.484 | +0.158 | 0.901 | 0.456 | +0.070 |
| E_had todos | 1.023 | 0.454 | +0.091 | 0.965 | 0.423 | +0.033 |
| **E_ν ν_μ CC** | **0.828** | 0.426 | −0.172 | **0.792** | 0.425 | −0.208 |

Que CLUE3D quede más cerca de la cota en ν_μ CC no es separación del muón: es que los outliers
de CLUE quitan un ~6 % de energía, y en NC (sin muón) eso mismo deja la E_had un 6 % por debajo de
la cota. El sesgo de E_ν empeora 3.6 puntos.

Barrido de parámetros (300 ev por punto, `.../scratch/clue_scan/`):

| dc [mm] | ρc [GeV] | 3D/ev | pureza | sin agrupar (E) | E_had/E_true | E_ν/E_ν,true |
|---|---|---|---|---|---|---|
| 10/15/40 | 0.02/0.04/0.2 | 17 | 0.641 | 5.2 % | 0.954 | 0.778 |
| 5/7.5/20 | 0.02/0.04/0.2 | 28 | 0.674 | 8.1 % | 0.924 | 0.758 |
| 2.5/5/10 | 0.02/0.04/0.2 | 39 | 0.709 | 11.5 % | 0.899 | 0.743 |
| 5/7.5/20 | 0.01/0.02/0.1 | 45 | 0.680 | 4.3 % | 0.965 | 0.784 |

La pureza sube poco al cerrar dc y se paga en energía sin agrupar. En una interacción de neutrino
en el SND los hadrones salen en un cono estrecho y se solapan: 0.64–0.71 parece el techo práctico
de la separación por densidad con esta granularidad.

**El enlace al vértice del guión (eje a < 30 mm) no funciona con CLUE3D**: sólo el 20 % de la
energía de los clusters 3D lo cumple (47 % a < 100 mm, 82 % a < 300 mm). El guión lo pensó para el
cono, cuyo eje pasa por el vértice por construcción; el eje ajustado a un trozo de cascada en el
hierro no apunta tan bien. Además, en una muestra sólo GENIE toda la energía es de la interacción,
así que el enlace no tiene nada que rechazar. Queda desactivado (`VertexLinkMM = 1e9` en job6)
hasta tener una muestra con fondo.

**De dónde no viene el +16 % de ν_μ CC**: los clusters 3D dominados por el muón primario sólo
tienen el 1.4 % de E_had (mediana, 300 ev). El depósito del muón está mezclado dentro de clusters
hadrónicos, así que no se puede quitar eliminando clusters 3D enteros; lo que lo ataca es restar
el depósito MIP esperado a lo largo de la trayectoria (opción (a) de las preguntas abiertas).

## Fase 4.2b — resta del muón a lo largo de su trayectoria (2026-10-14, `pfa_ref_a`, 2000 ev)

Opción (a) de las preguntas abiertas, pedida por el usuario el 2026-10-14, junto con la comparación
capas vs CLUE3D con la corrección aplicada en los dos.

**Bug encontrado y corregido antes de comparar** (afecta a las tablas de 4.2 y 4.3): `SNDPFOBuilder`
y `SNDHadronicClusterer` identificaban clusters e hits por `(collectionID, índice)`, pero una
colección creada en el mismo job no tiene collectionID hasta que se escribe (todas valen
0xffffffff). El índice 3 de los clusters del SiTarget era así "el mismo" que el índice 3 de los del
MTC: el corredor marcaba clusters de otras colecciones, CLUE3D excluía los equivocados y los
enlaces de verdad de los PFOs mezclaban hits. Ahora la clave es `podio::detail::getOrderKey`
(la dirección del objeto, lo que compara el propio `operator<` de podio).

**Método** (`MuonSubtraction = "trajectory"`, default de job6; `PFA_MUON_SUBTRACTION=corridor`
recupera el 4.1). Se recorre la traza ajustada del muón (estados AtOther, interpolados y
extrapolados recto) por todos los planos sensibles desde el vértice; en cada uno, el cluster más
cercano a la trayectoria (ventana 5/15/20/50 mm en SiTarget/SiPad/SciFi/tiles) cede
`MuonMeanMIP × (E/MIP del cluster) × 1/cosθ`, como mucho lo que tiene. `MuonMeanMIP` es la media de
MIPs por plano del muón de 10 GeV de la gun (1.78/1.81/1.55/1.67/1.73/1.36 para SiT X/Y, SiPad,
U/V, tiles). El 4.1 restaba 1 MIP, pero 1 MIP es el **pico** de Landau (`MIPEnergyOverride`)
y la media por plano es de 1.35 a 1.8 (el 4–10 % de los planos se llevan un δ o una pérdida
radiativa), y además rechazaba los clusters de > 3 MIP, justo esa cola. CLUE3D recibe ahora
también los clusters del muón (`ExcludeMuonCorridor = False`) y la resta sigue a cada cluster
de capa al cluster 3D en que acabe. Estéreo del MTC comprobado con muones de verdad:
s_U = x − y·tanα, s_V = x + y·tanα.

E_had: diferencia con la cota por clase; E_ν: ratio en ν_μ CC (con el bug arreglado en todas):

| | corredor + capas | corredor + CLUE3D | **trayectoria + capas** | trayectoria + CLUE3D |
|---|---|---|---|---|
| E_had NC | −0.001 (σ 0.332) | −0.060 (σ 0.364) | **−0.001** (σ 0.330) | −0.061 (σ 0.362) |
| E_had ν_μ CC, muón encontrado | +0.199 (σ 0.642) | +0.118 (σ 0.548) | **+0.099** (σ 0.432) | +0.040 (σ 0.441) |
| E_had ν_μ CC, sin muón | +0.158 | +0.070 | +0.158 | +0.070 |
| sesgo E_had (todos) | +0.031 | −0.029 | **+0.004** | −0.059 |
| E_ν ν_μ CC, muón encontrado | 0.960 (σ 0.339) | 0.937 (σ 0.336) | **0.933 (σ 0.330)** | 0.913 (σ 0.327) |
| E_ν ν_μ CC, todos | 0.832 | 0.795 | 0.817 | 0.782 |
| E_ν, p_μ por curvatura | 1.002 (σ 0.230) | 0.979 | 0.977 (σ 0.212) | 0.954 |

CLUE3D: 19.4 clusters 3D/evento, pureza ponderada en energía 0.72 (era 0.64 con el bug), 7.4 % de
la energía sin agrupar.

- La resta por trayectoria **reduce el exceso a la mitad (+0.199 → +0.099) y mejora la resolución
  de E_had en ν_μ CC un 33 % (σ68 0.642 → 0.432)**, sin tocar NC. Con muón por curvatura, σ68 de
  E_ν baja de 0.230 a 0.212.
- Se resta el **84 %** del depósito verdadero del muón (mediana 2.26 de 2.75 GeV; σ68 0.15); el
  cluster de la trayectoria se encuentra en el 90 % de los planos. En los 32 NC con muón falso se
  quitan 0.96 GeV.
- Lo que falta es dependencia en energía: por plano, los muones de GENIE dejan de media 1.5–1.8 MIP
  entre 10 y 30 GeV y **2.2–2.3 MIP por encima de 30 GeV** (SiTarget, SciFi), un 25–30 % más que la
  gun de 10 GeV, por las pérdidas radiativas. **Guns de muón de 3, 30 y 100 GeV enviadas a
  HTCondor** (clusters 20913929–31) para hacer `MuonMeanMIP` función de p_μ.
- CLUE3D vs capas con la misma corrección: CLUE3D queda más cerca de la cota en ν_μ CC (+0.040), pero
  por la misma razón que en el 4.3: sus outliers quitan un ~6–7 % de toda la energía, y en NC
  eso lo deja en −0.061. Con la resta por trayectoria la ventaja que se le suponía en 4.3
  (separar el muón) ya la tiene la suma por capas. **E_ν sigue mejor con capas** (0.933 vs 0.913
  con muón encontrado). El default sigue siendo capas.
- El criterio del 4.2 (|sesgo E_ν| < 10 % en todos los ν_μ CC) **sigue sin cumplirse** (−0.183):
  el 36 % de ν_μ CC sin muón (ratio 0.43) y la cota de hierro (0.27) dominan, y eso no lo cambia
  ninguna resta.
- Hallazgo lateral, **no corregido** (es de la fase 3): el corredor de `SNDMuonFinder` compara la s
  de las vistas U/V del SciFi con la x de la traza sin el término ±y·tanα; a |y| de 100–200 mm eso
  es un error de 9–17 mm con una ventana de 10 mm. Puede costar eficiencia de muón (opción (b)).

Ficheros: `.../snd_pfa/genie/pfa_pfa_ref_a_{corr,corr_clue,traj,traj_clue}.root` con su
`enu_perf_ref_a_*.txt` y `job6_pfa_ref_a_*.log`.

### 4.2b, closed — momentum dependence and the per-plane cap (2026-10-14)

*(Entries from here on are in English.)*

**Muon guns at 3/30/100 GeV** (500 events each, reconstructed with `run_calib.sh`). Mean MIPs per
crossed plane; the median barely moves (1.09–1.19) while the radiative tail grows (planes above
3 MIP: 6 % at 3 GeV, 15 % at 100):

| p (GeV) | SiT X | SiT Y | SiPad | MTC U | MTC V | tiles |
|---|---|---|---|---|---|---|
| 3 | 1.529 | 1.559 | 1.397 | 1.526 | 1.563 | 1.261 |
| 10 | 1.782 | 1.810 | 1.548 | 1.665 | 1.723 | 1.350 |
| 30 | 1.956 | 1.998 | 1.916 | 1.882 | 1.929 | 1.508 |
| 100 | 2.875 | 2.919 | 2.738 | 2.490 | 2.571 | 1.928 |

`MuonMeanMIP` is now this table (`MuonMeanMIPMomenta` = 3, 10, 30, 100), interpolated linearly in
log p at the reconstructed p_μ and held flat outside it. **Alone it changed nothing**: E_mu,sub /
E_mu,true 0.837 → 0.849, the E_had excess +0.099 → +0.101 (`pfa_pfa_ref_a_trajp`).

**The real cause is the per-plane cap.** Taking min(mean, cluster) from each plane is biased low
whatever the mean: the deposit is a Landau, so the typical plane holds less than the mean and a few
hold far more. On the guns, min(mean, X) takes only **56–83 %** of the deposit (`mu_cap.py`).
Two fixes were tried:

- *Take isolated clusters whole* (no other cluster on the plane within 20/40/60/150 mm) — **rejected**:
  it took 4.5 GeV for a true 2.75 (σ68 1.5). In the scintillator 93 % of the muon clusters looked
  "isolated" because 2D connectivity merges the whole shower, muon included, into one cluster with no
  neighbours. Even restricted to small clusters it cannot reach the tail: at 100 GeV, "whole if
  ≤ 20 MIP, else mean" still takes only 70–86 %, because a radiative loss is a real shower.
- *Take the remainder from the event* — **adopted** (`MuonSubtraction = "trajectory"`, the default;
  the first version is kept as `"trajectory-capped"`). Each cluster still gives at most what it holds
  (so the per-cluster split for CLUE3D is unchanged), and what it could not give is subtracted from the
  event's silicon / MTC totals. Summed over a whole track this is the unbiased mean deposit.

Results (layer sum, 2000 events; `pfa_pfa_ref_a_{unc,uncst,uncst_clue}`):

| | capped (`trajp`) | **uncapped (`unc`)** | uncapped + corridor stereo fix (`uncst`) | same, CLUE3D |
|---|---|---|---|---|
| E_mu,sub / E_mu,true (median, σ68) | 0.849, 0.150 | **1.054, 0.208** | 1.051, 0.212 | 1.051, 0.212 |
| E_had gap, NC | −0.001 | −0.001 | −0.007 | −0.063 |
| E_had gap, ν_μ CC muon found (σ68) | +0.101 (0.441) | **+0.001 (0.385)** | +0.002 (0.372) | −0.061 (0.380) |
| E_had gap, all | +0.068 | +0.020 | +0.018 | −0.046 |
| E_ν ratio, ν_μ CC muon found (σ68) | 0.933 (0.332) | 0.903 (0.311) | 0.902 (0.314) | 0.886 (0.318) |
| E_ν ratio, all ν_μ CC | 0.817 | 0.807 | 0.810 | 0.777 |

- **The muon excess in E_had is gone (+0.101 → +0.002) and its σ68 improves another 16 %** (0.441 →
  0.372); NC is untouched. In NC events with a fake muon (a pion corridor) 0.92 GeV is taken.
- E_ν with the muon found drops to 0.90: the E_had excess used to hide part of the iron leakage.
  The remaining E_ν bias is now leakage (bound 0.933 on E_had) and p_μ, not the muon deposit.
- The layer sum stays the default; CLUE3D is still behind for E_ν.
- The 4.2 criterion (|E_ν bias| < 10 % on all ν_μ CC) **still fails (−0.190)**, dominated by the ν_μ CC
  without a found muon and by the leakage. It is no longer a muon-subtraction problem.

## Muon efficiency study (2026-10-14, `pfa_ref_a`, 2000 events, 1065 ν_μ CC)

Asked by the user after the corridor stereo fix (row 3.4). New script
`Optimization/Analysis/pfa/muon_loss.py`: puts every ν_μ CC at the first stage it fails, using a second
PFA file made with the finder's cuts opened up (`MinMTCLayers 0`, `MinFracLayers 0`, `MaxMeanMIP 1e9`)
to read which cut rejects a muon that has a track.

**Where the muon is lost** (with the stereo fix, before any retuning):

| stage | all | p < 5 GeV | 5–20 | > 20 |
|---|---|---|---|---|
| acceptance: the muon lights < 8 MTC layers (truth) | 0.321 | 0.578 | 0.274 | 0.182 |
| tracking: ≥ 8 layers but no track dominated by it | 0.036 | 0.024 | 0.041 | 0.039 |
| finder: has a track, no candidate | 0.008 | 0.007 | 0.014 | 0.002 |
| **found** | **0.636** | 0.391 | 0.672 | 0.776 |

**Given acceptance the finder already finds 94 %.** The loss is geometry: of the 33 % of muons crossing
fewer than 8 MTC layers, 8.6 % cross none (40 % of those from SiTarget vertices, mostly p < 3 GeV or
θ > 30°), 11.8 % cross 1–3 and 12.9 % 4–7, half of these from vertices in MTC60 — the muon simply
leaves through the back. Only 76 of those 355 have a track at all.

**Two finder changes** (`SNDMuonFinder`):

- `MinMTCLayersAtExit` (new, off by default): a corridor that reaches the last MTC layer
  (within `ExitToleranceLayers` = 1) needs only that many MTC layers. It uses the `m_lastMTCLayer`
  that was computed but never used.
- `CorridorHalfWidthMTCSciFi` re-scanned now that U/V are compared at the right s.

| MinMTCLayersAtExit / SciFi half-width | eff. ν_μ CC | among muons with a GlobalTrack | fakes NC | fakes ν_e CC |
|---|---|---|---|---|
| off / 10 mm (current job6) | 0.636 | 0.877 | 0.076 | 0.093 |
| 6 / 15 | 0.655 | 0.905 | 0.088 | 0.104 |
| 6 / 10 | 0.656 | 0.907 | 0.082 | 0.095 |
| 4 / 6 | 0.661 | — | 0.076 | 0.095 |
| 6 / 6 | 0.654 | 0.904 | 0.070 | 0.090 |
| 6 / 4 | 0.653 | 0.903 | 0.056 | 0.088 |
| **6 / 3** | **0.653** | **0.903** | **0.056** | **0.081** |
| 6 / 2 | 0.647 | 0.895 | **0.049** | 0.069 |
| 4 / 10 | 0.664 | 0.917 | 0.088 | 0.100 |
| 3 / 10 | 0.665 | 0.919 | 0.095 | 0.102 |

- Exit rule at 6 + a narrow SciFi window **improves efficiency and fakes at the same time**: before the
  stereo fix the window had to be wide to absorb the missing ±y·tanα, and it let pions in.
- 6 / 2 is the first point that meets the plan's NC fake criterion (≤ 0.05) — with efficiency still
  above today's. Neither point meets "≥ 0.97 of the muons with a GlobalTrack" (0.90).
- E_ν does not move (ν_μ CC −0.191 at 6 / 3): the muons recovered are few and soft.
- **Set as the job6 default on 2026-10-14 (user): 6 / 2** (`MinMTCLayersAtExit = 6`, `CorridorHalfWidthMTCSciFi = 2.0`).

**What is left, by size:** acceptance (32 %, physics and geometry: a short muon with no MTC lever arm
looks like a pion); tracking with acceptance (3.6 %: a calorimeter-seeded muon finder in the MTC could
recover part of it); finder (~1 %). Against the tracker the finder still costs ~10 points (plan
budget: 3), but those are muons that DO have a GlobalTrack while crossing fewer than 8 MTC layers —
the same acceptance term, where a looser cut buys efficiency only with pion fakes (rows 4/10, 3/10).

Files: `.../snd_pfa/genie/pfa_pfa_ref_a_{uncst,loosest,exit*}.root`, `muon_perf_ref_a_*.txt`,
`muon_loss_ref_a_isost.txt`.

## Muestras y ficheros producidos

| Qué | Ruta | Fecha |
|---|---|---|
| GENIE 2000 ev, overlay con MCParticles | `gaudi_jobs/genie_pipeline/events.edm4hep.root` | |
| GENIE tracking con q/p | `gaudi_jobs/genie_pipeline/tracks_seq.edm4hep.root` | |
| GENIE PFA (con PFOs) | `/eos/experiment/ship/user/marquezh/snd_pfa/genie/pfa_pfa_ref_a.root`, 2000 ev; log `job6_pfa_ref_a.log`. Con E_had por CLUE3D: `pfa_pfa_ref_a_clue.root` | 2026-10-13 |
| Guns de calibración | `/eos/experiment/ship/user/marquezh/snd_pfa/calib/` — 11 jobs lanzados | 2026-09-23 |
| GENIE de referencia (flujo correcto) | `simulation/run_script/data/output_genie_pfa_ref_{a,b,c,d}.edm4hep.root`, 4×2000 ev, 651 MB | 2026-09-22 |
| Área EOS para pruebas | `/eos/experiment/ship/user/marquezh/snd_pfa/{genie,calib,scratch}/` | 2026-09-23 |
| Scripts de análisis PFA | `Optimization/Analysis/pfa/` (`pfa_truth.py`, `check_truth.py`, `check_qop.py`, `check_decoder.py`, `check_calohits.py`, `enu_perf.py`) | 2026-09-23 |
| Variante de geometría para el test de pesos | `Optimization/Simulation/geometry/variants/pfa_wtest.xml` (SiTarget W 7.0, SiPad W 5.0) | 2026-09-23 |

## Desviaciones del guión (fecha · qué · por qué)

- 2026-10-13 · Paso 4.3: **CLUE3D en lugar del enlace por cono** (decisión del usuario). Se
  implementa el algoritmo en `SNDClue3D.h` en vez de usar la librería `cluestering` del stack:
  ésta es genérica en N dimensiones sobre alpaka, con una sola distancia crítica para todos los
  puntos y sin ventana de capas, y el SND necesita valores por detector. Cambios respecto a CLUE3D
  de TICL: dc y ρc por detector; distancia sólo transversal (la ventana de capas hace el trabajo
  longitudinal, porque cada detector tiene un paso de capa distinto); sin tiles (O(N²), con ~200
  puntos por evento). Las vistas X/Y (SiTarget) y U/V (MTC) se emparejan en puntos (x, y) con el
  criterio de asimetría del guión (< 0.5); una vista sin pareja toma la otra coordenada del
  baricentro de la otra vista de la capa, o del evento.
- 2026-10-13 · Paso 4.3: **el enlace al vértice queda desactivado** (`VertexLinkMM = 1e9` en
  job6) y **la E_had por clusters 3D es opcional**, no la de por defecto: ver la sección
  "Fase 4.3". Con 3D activado se escribe un PFO por cluster (`goodnessOfPID` = 1 si enlazado al
  vértice), y cada uno se enlaza a su MCParticle dominante con la pureza como peso, en lugar de
  enlazar el sistema hadrónico al neutrino.
- 2026-10-13 · Paso 4.3: `SNDMuonChoice.h` reúne la regla "el muón del evento es el candidato con
  más capas del MTC", que usan tanto `SNDHadronicClusterer` como `SNDPFOBuilder`. Si cada uno la
  tuviera por su cuenta, uno podría restar un corredor que el otro cuenta.
- 2026-10-13 · Paso 4.1: **`EoverH` se separa en `EoverHECAL` = 1.075 y `EoverHMTC` = 1.0**, no
  un único 1.4. El 1.4 del guión es anterior a la calibración: los hits del silicio llegan a
  escala EM (`GlobalScaleECAL`, ajustada con electrones) y los del MTC ya a escala hadrónica
  (`GlobalScaleMTC` se ajustó con piones). Medido con `calib.py` y las escalas actuales:
  pendiente e 1.000, π 0.930 → e/h = 1.075. Con 1.4 en ambos el MTC habría salido un 40 % alto.
- 2026-10-13 · Paso 4.1: **los `pfa_*` no van como parámetros de frame** sino en la colección
  `PFASummary` (un `ParticleID` por evento, `PDG` = `pfa_flavour`, `particle` = el PFO
  hadrónico), con el orden de los valores en el metadato `PFASummaryNames`. Con IOSvc un
  algoritmo sólo puede escribir metadatos de run (`k4FWCore::putParameter`), no parámetros de
  evento.
- 2026-10-13 · Paso 4.1: el PFO hadrónico no lleva `type = 1`: EDM4hep 1.x quitó
  `ReconstructedParticle::type`. Se distingue por `PDG = 0` (el muón es ±13).
- 2026-10-13 · Paso 4.1: **sólo el corredor del muón elegido** (el de más capas del MTC) se
  excluye de E_had. Los clusters de los demás candidatos vuelven al sistema hadrónico: en 3.3 esos
  corredores eran π± que cruzaban el MTC, y quitarlos hacía desaparecer su energía.
- 2026-10-13 · Paso 4.1: la corrección del corredor `(E_MIP − 1) × E_activo × w` se calcula como
  `E_cluster × (MIP − 1)/MIP`, que es lo mismo (la energía de un cluster es lineal en sus MIP) y
  ya incluye la escala global.
- 2026-10-13 · Paso 4.1: el vértice `track` usa las trazas con el bit 1 (SiTarget) de
  `SiTargetTracks` **y** `GlobalTracks`, porque las trazas empalmadas sólo están en la segunda.
  Si `VertexMode = track` no encuentra traza, cae al de clusters. En el MTC la `y` del vértice
  por clusters sale del par U/V con el ángulo estéreo (`Y = (s_V − s_U)/(2 tan α)`, igual que
  `MTCSciFiMeasConverter`), o del centelleador si falta una vista.

- 2026-09-26 · **Bug propio, encontrado por la comprobación de anchura de la Fase 2.3**: en
  `SNDLayerClusterer` y `SNDTruthClusterer` guardaba `&hit` de un objeto entregado por el
  iterador de una colección podio. Ese objeto es un **handle temporal**, no un elemento
  persistente: el puntero apunta a un buffer que la iteración siguiente reutiliza, así que
  todos los hits guardados se leían como el mismo (clusters con 3–5 hits de un único
  cellID y posición idéntica). Los objetos podio se copian **por valor**. Se detectó porque
  la anchura de los clusters salía 0.00 mm en los cuatro detectores teniendo el 30 % más de
  un canal. Corregido; las cifras de 2.1–2.3 son posteriores al arreglo.
- 2026-09-26 · La anchura de un cluster de lectura **2D** se mide en las dos coordenadas
  (`hypot`): dos tiles contiguos en y comparten x, así que una extensión sólo en x daba 0.

- 2026-09-23 · **La muestra GENIE de referencia del guión no servía, por una razón nueva
  además del flujo.** `output_genie_allflavours_0_2000.edm4hep.root` (31-ago) se simuló
  con el compact de **166 capas de SiTarget**, el que el informe del 11-sep identificó
  como erróneo: sus hits llegan a `layer=165` y la geometría actual tiene 120 (0–119).
  Resultado: el 12.8 % de los hits del SiTarget apuntan a volúmenes que ya no existen y
  `CellIDPositionConverter` no los puede colocar (1929 de 15071 hits en 50 eventos). Se
  detectó porque, al bajar el umbral efectivo con los MIP medidos, esos hits dejaron de
  caer por debajo del corte y el `VolumeManager` empezó a quejarse. **El pipeline GENIE
  apunta ahora a `output_genie_pfa_ref_a.edm4hep.root`** (2000 ev, geometría y flujo
  correctos: layer máximo 119, 9.78 M contribuciones frente a 3.73 M — el flujo bien
  colocado deposita ~3× más). Los números de los pasos 0.2/0.3 y de la Fase 1 medidos
  sobre la muestra vieja **hay que rehacerlos**; los cualitativos (signo de q/p invertido,
  pesos de muestreo, escala EM de las guns) no dependen de ella.
- 2026-09-23 · `SNDCaloHitBuilder` captura el fallo de `CellIDPositionConverter` y cae a
  la posición simulada, contándolo y avisando en `finalize()`. Con la muestra correcta el
  contador debe quedar **en 0**; si no, es señal de que la muestra y la geometría no
  casan.
- 2026-09-23 · Las salidas de la cadena GENIE pasan a EOS
  (`/eos/.../snd_pfa/genie/{events,tracks,pfa}_pfa_ref_a.root`), con `INPUT_FILE` /
  `OUTPUT_FILE` como variables de entorno en job1/job4/job5/job6.

- 2026-09-23 · Paso 1.1/1.2: **el peso de muestreo es POR CAPA, no por plano** (decisión
  del usuario tras preguntar por qué MTC V no aparecía en las comprobaciones). Con el
  esquema del guión (todo el absorbedor al primer plano sensible que lo sigue) la segunda
  vista quedaba con peso ~1: medido en 50 eventos, MTC V tenía **más** hits que U (2307 vs
  2251) y aportaba el **0.45 %** de la energía frente al 55 % de U; SiTarget Y, el 0.17 %.
  Eso tira la información de la mitad de los canales (dos vistas ortogonales muestrean el
  mismo paso: promediarlas vale √2 en la fluctuación) y hace que una capa pierda el 99.5 %
  de su energía si el plano pesado cae bajo umbral. Ahora `w = (pasivo_capa + Σactivo)/
  Σactivo` para todos los planos de la capa. Efecto: X/Y coinciden al 1.5 % y U/V al 0.6 %.
- 2026-09-23 · Paso 1.1: **bug encontrado al comprobar la periodicidad**. El barrido TGeo
  rompía cuando `FindNextBoundary` devolvía un paso nulo, lo que ocurre justo cuando el
  punto de partida está sobre una frontera — y los extremos del intervalo son, por
  construcción, las caras de un plano sensible. Resultado: el SiTarget salía con **dos**
  pesos (33.91 en 63 capas, 47.06 en 54) siendo una geometría periódica, y el MTC con
  tres. Arreglado arrancando 1 nm dentro y cruzando la frontera nula en vez de parar.

- 2026-09-23 · Paso 1.1: **no se usa `dd4hep::rec::MaterialManager::materialsBetween`**,
  como decía el guión, sino un barrido TGeo explícito (`InitTrack`/`FindNextBoundary`/
  `Step`) en `SNDSamplingWeights.h`. Motivo: sobre esta geometría `materialsBetween`
  devolvió sólo los últimos 3 mm de los 53 mm de hierro entre dos capas del MTC (tramo
  z = 81..134 mm), mientras acertaba en otros tramos; el barrido TGeo da 53.000 mm ahí y
  reproduce el compact en los tres detectores. Además `Vector3D` va en cm en esta build
  (`dd4hep::mm = 0.1`), lo que hacía que un valor en mm cayera fuera del mundo de 20 m.
- 2026-09-23 · Paso 1.1: los planos sensibles **no** se reconocen por nombre de slice sino
  por `dd4hep::Volume::isSensitive()` + deduplicación por (det, view, z). La sensibilidad
  vive en sitios distintos según el detector: el MTC la marca en el slice (42/estación),
  el SiTarget en cada tile `_colN_rowM` (1920 volúmenes → 240 planos) y el SiPad en cada
  `_wafer_pads` (792 → 22).

- 2026-09-22 · Paso 0.3: **el criterio de aceptación del guión estaba mal planteado**.
  Pedía |sesgo| < 15 % de `p_reco` frente al momento del muón **en el vértice**, pero
  `AtLastHit` está dentro del MTC, después de ~4 GeV de pérdida por ionización (42 cm de
  W en el SiTarget + 22 cm en el SiPad + hasta 2.2 m de Fe). Medido: dp/p vs vértice
  −0.29/−0.29/−0.24/−0.14 por bins, y vs `momentumAtEndpoint` +0.42/+0.16/+0.08/−0.02 —
  es decir `p_reco` cae **entre** ambos, como debe. Criterio nuevo: signo de carga > 80 %,
  `p_reco` entre endpoint y vértice en > 80 % de los casos, y |dp/p vs endpoint| < 15 %
  para p > 10 GeV. La corrección de dE/dx es trabajo de `SNDMuonFinder` (Fase 3) con la
  tabla de material de la Fase 1.
- 2026-09-22 · Paso 0.3: añadida la relación **Track → TrackerHit3D** en
  `ACTSSequentialTracker::writeTrack` (aditiva; no cambia la reconstrucción). Sin ella no
  hay forma de saber qué partícula MC construyó una traza salvo parseando las líneas
  `[SegDump]` del log, y la Fase 3 necesita ese enlace en tiempo de ejecución.
  `SeqMeasurement` gana un `srcIndex` para ello.
- 2026-09-22 · Paso 0.3: **corregido el signo de q/p al escribirlo**. El marco de tracking
  mapea DD4hep (x,y,z) a ACTS (z,y,x), que es una permutación impar (reflexión); B es un
  pseudovector, así que el fit, alimentado con el By de DD4hep y con medidas igualmente
  reflejadas (loc0 = −ddX), es autoconsistente en trayectoria y |q/p| pero deduce la
  carga espejada. Medido antes del cambio: q/p < 0 en el 99 % de los μ+ verdaderos y en
  sólo el 13 % de los μ−. El fit **no se toca**; sólo se invierte el signo al escribir.

- 2026-09-22 · Paso 0.2: las relaciones padre/hija **sí** se preservan también con N
  fuentes (remapeadas con los offsets, en dos pasadas), no sólo con una fuente como decía
  el guión. Cuesta poco y hace falta para reconocer al leptón primario por su madre.
- 2026-09-22 · Paso 0.2: añadido un parámetro de frame extra `MCParticleSourceIDs`
  (source_id de cada offset, en el mismo orden que `MCParticleSourceOffsets`). Sin él un
  consumidor no puede mapear el source_id de una contribución (que EventOverlay guarda en
  `CaloHitContribution::PDG`) a su offset.

## Preguntas abiertas para el usuario

- **Paso 4.2 no cumple |sesgo E_ν| < 10 % (−17 %) — 2026-10-13.** Diagnóstico en la sección
  "Fase 4.2". El PFO builder no es el cuello de botella: con muón por curvatura E_ν sale sin
  sesgo. Opciones: (a) **restar el depósito MIP esperado del muón a lo largo de su trayectoria**
  (1 MIP × w por plano cruzado, con la tabla de material de la Fase 1) en vez de fiarse de los
  clusters que recoge el corredor; ataca el +16 % de E_had en ν_μ CC y no depende de los cortes
  del corredor; (b) mejorar la eficiencia del muón (36 % de ν_μ CC sin muón), que según 3.3 pasa
  por el tracking del MTC; (c) aceptar y reportar el sesgo separado por `pfa_p_mu_source`.
  El paso 4.3 (CLUE3D, 2026-10-13) no lo resuelve: el depósito del muón está mezclado dentro de
  clusters hadrónicos, no en clusters propios. **(a) sigue siendo la recomendación.** → Implementada el 2026-10-14 (sección "Fase 4.2b"): reduce el exceso de E_had a la mitad pero no cambia el criterio de E_ν, que dominan (b) y la cota de hierro.
- **CLUE** — CLUE3D implementado en el 4.3 (ver arriba). Lo que sigue se refiere a sustituir
  también el clustering 2D de capa (comentado con el usuario el 2026-10-13). El stack
  trae `k4clue` 01-00-09 (envoltorio Gaudi para barril + endcaps, en principio no aplicable tal
  cual) y la librería `cluestering` 2.7.0, que sólo trae headers. Si en 4.3 o en la Fase 6 hace
  falta separar cascadas solapadas, un `SNDClueClusterer` con la misma salida que
  `SNDLayerClusterer` podría sustituirlo sin cambiar nada aguas abajo; evaluarlo contra
  `SNDTruthClusterer` con `layer_perf.py`/`muon_perf.py`. CLUE3D sería la alternativa al
  enlace por cono del 4.3.

- ~~El momento del muón por rango no es viable en el SND~~ — **resuelto 2026-09-27**: el
  usuario eligió **curvatura + rango como cota inferior**, ya implementado. Se conserva el
  diagnóstico porque explica por qué. Números: el MTC entero tiene
  **2352.75 mm de hierro**, cuyo poder de frenado sólo detiene muones de hasta
  **2.69 GeV**, mientras los muones de GENIE van de 2 a >100 GeV. **Ningún muón se para
  dentro**: lo que mide `ironMm` no es el rango del muón sino cuántas capas reconstruyó el
  tracking (mediana 476 mm a 2 GeV → 1960 mm a >20 GeV), es decir una medida de la calidad
  del tracking disfrazada de momento. De ahí el sesgo de −0.55 a −0.85. Y sólo el 32 % de
  los muones de >20 GeV llegan a la capa 42, así que el criterio "sale del detector" que
  conmuta a curvatura falla en dos tercios de los casos. Opciones: (a) **curvatura como
  medida principal** (σ68 40–50 %, sin sesgo grande) con el rango como **cota inferior**
  p > 0.0011427 × ironMm, que sí es información válida; (b) curvatura sola, sin rango;
  (c) ajustar una función p(ironMm, nLayers) empírica sobre las guns de muón, aceptando
  que mide tracking y no rango. **Elegida (a).**
- **Limitación de diseño del MTC, anotada** (a petición del usuario, sin proponer cambios
  al espacio de búsqueda del MOBO): con 2352.75 mm de hierro el MTC (a) no puede medir
  momento por rango para el espectro de muones que le llega, y (b) no contiene las cascadas
  hadrónicas de las interacciones que ocurren en su mitad trasera — fuga del 36 % medida en
  el paso 1.5. Ambas cosas dependen de parámetros del MTC que hoy están fijos en el MOBO.

- **Respuesta hadrónica no lineal (paso 1.4).** E_vis/E_true crece monótonamente de 0.521
  (π de 2 GeV) a 0.620 (40 GeV): +19 %. Es lo normal en un calorímetro no compensado (a
  más energía, más fracción EM en la cascada), y el residuo del 15.5 % frente a un ajuste
  de un solo parámetro **no se puede quitar con `EoverH`**. Opciones para la Fase 4:
  (a) corrección dependiente de la energía `E_had = f(E_vis)` ajustada con estas guns;
  (b) *software compensation* por densidad de energía por hit (pesar más los hits de baja
  energía, que son la componente hadrónica); (c) aceptar el sesgo y reportarlo. Recomiendo
  (a) como primer paso: es barato y las guns ya están.

- **Muestra de referencia GENIE** — *en curso*: decidido el 2026-09-22 **generar una
  muestra correcta antes de seguir con la Fase 1**. `output_genie_allflavours_0_2000.edm4hep.root` (la que
  usa `job1_overlay.py` y el informe del 2026-09-11) tiene los vértices mal colocados: de
  300 eventos, mediana z = +625 mm (dentro del MTC) y **sólo el 22 % dentro del SiTarget**
  (z −1690..−370). Es el bug de flujo documentado en la memoria `genie-flux-misplaced`.
  `output_genie_gen_test_1000.edm4hep.root` (17-sep, 999 eventos) tiene el mismo problema:
  mediana +906 mm, 19 % en el SiTarget. Toda la física de las Fases 1.5, 3 y 4 (E_had,
  E_ν, eficiencia de muón) se mediría sobre interacciones que ocurren en el sitio
  equivocado.
- ~~Criterio de momento del muón (paso 0.3)~~ — **resuelto 2026-09-22**: rango como
  medida principal en la Fase 3, curvatura para carga y muones que escapan. El guión
  (§7, paso 3.2) está actualizado.
