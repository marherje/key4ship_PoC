# SND-PFA — guión de implementación por órdenes

Particle flow propio para SiTarget + SiPad + MTC, con la arquitectura por etapas de CyberPFA
(CEPCSW, `Reconstruction/RecPFACyber`). Escrito el 2026-09-21 para que Opus lo implemente
en varias sesiones. **Este fichero es el maestro y no se edita** (salvo errores de hecho);
el estado se lleva en `docs/snd_pfa_progress.md`.

## Cómo usar este guión

1. Al empezar una sesión: lee este fichero entero y `docs/snd_pfa_progress.md`.
2. Continúa por el **primer paso sin marcar**. Los pasos son secuenciales dentro de cada
   fase; no empieces una fase sin cerrar la anterior.
3. Cada paso tiene tres bloques: **Haz** (órdenes), **Comprueba** (qué ejecutar y qué debe
   salir), **Anota** (qué escribir en el registro). Un paso está cerrado sólo cuando los
   tres están hechos.
4. Si un criterio de aceptación no se cumple: anótalo en el registro, **no bajes el
   listón**, y pregunta al usuario antes de seguir.
5. Cualquier desviación del guión va a la sección "Desviaciones" del registro, con fecha y
   motivo.
6. Al cerrar un paso propón un commit al usuario (no commitees sin preguntar).

## A. Por qué este diseño (contexto, no órdenes)

- CyberPFA no se porta: está atado a la geometría CEPC (`CyberPFAlg.cpp:98-111`,
  `CrystalEcalSvc`, `GearSvc`, LCIO, EDM4CEPC, helices desde un IP). Pandora
  (`k4gaudipandora` + `lccontent`, en el stack) tampoco: `LCContent` asume barril+endcaps y
  trazas desde un IP; adaptarlo al SND (blanco fijo, calorímetro plano) es escribir una
  *content library* nueva. Queda como "camino B" documentado.
- Se copia de CyberPFA: (1) la cadena hits calibrados → clusters 1D por vista (*half
  clusters*, que aquí son los planos X/Y del SiTarget y U/V del MTC) → clusters 3D →
  traza–cluster → PFOs; (2) una versión **Truth** de cada etapa como cota; (3) un ntuple de
  depuración por etapa.
- Decisiones del usuario (2026-09-21): un `Gaudi::Algorithm` por etapa con colecciones
  EDM4hep entre ellas; entrada `*HitsWindowed` + umbral; prioridad (1) muón ID + momento,
  (2) E_had y E_ν. Cascadas EM: fase opcional al final.
- Físicamente, SiTarget+SiPad son ≈180 X₀ / ≈6 λ de W/Si: ahí se queda casi toda la
  energía EM y hadrónica; el MTC (≈13 λ Fe, campo en el hierro) ve muones y colas. Por eso
  "ECAL" = SiTarget+SiPad y "HCAL/muon" = MTC en todo el guión.
- Hechos del repo en los que se apoya cada orden: **Anexo 1**.

## B. Reglas fijas (cúmplelas en todos los pasos)

- Cumple `CLAUDE.md`. Los pasos de PFA leen `gaudi_jobs/genie_pipeline/tracks_seq.edm4hep.root`
  (2000 eventos GENIE; contiene hits `*Windowed`, medidas y trazas) y tardan
  segundos–minutos: **un** `k4run` local cada vez. Todo `ddsim` va a HTCondor
  (`simulation/run_script/launch_multiplePG.sh` → `generic_condor_multiplePG.sh`).
  Reejecuta job4 (tracking, ~21 min) sólo cuando el guión lo diga.
- Compila con `cmake --build build --target install -j2` (o `make -j2 install` en `build/`);
  **no** uses `build.sh` tal cual (lleva `-j16`).
- `source init_key4ship.sh` antes de nada; `k4run` desde el directorio de la pipeline.
- Tras añadir o cambiar una `Gaudi::Property`, reinstala (`make install`) o Python no la ve.
  Los booleanos no se pueden poner a `False` desde la línea de `k4run`.
- Scripts de validación: Python con `podio.root_io.Reader`, en `Optimization/Analysis/pfa/`.
  Cada script imprime una tabla y termina con `OK`/`FAIL` según los criterios del paso.
- **Las salidas de las pruebas van a EOS, no a AFS** (decidido el 2026-09-23; la cuota AFS
  es de 10 GB y estaba al 73 %). Base:
  `/eos/experiment/ship/user/marquezh/snd_pfa/`, hermana de la de las campañas MOBO, con
  `genie/` (muestras de referencia y las cadenas job1→job6 sobre ellas), `calib/` (guns de
  la Fase 1) y `scratch/` (pruebas de usar y tirar). En AFS se quedan sólo el código, los
  docs y los ficheros pequeños de resultados (tablas, `.db`). Ojo: el schedd rechaza `/eos`
  dentro de los ficheros de submit de Condor — ésa es la razón del par `trials_dir` /
  `submit_dir` del MOBO; el mismo patrón aplica aquí.
- Sigue el estilo del código existente: `Gaudi::Algorithm` con
  `execute(const EventContext&) const`, `DECLARE_COMPONENT`, propiedades con docstring,
  comentarios que explican el *porqué* (mira `MTCSciFiMeasConverter.cpp` como referencia).
- No toques `ACTSProtoTracker` (tiene un bug de signo sin corregir); SND-PFA consume sólo
  la salida de `ACTSSequentialTracker`.
- No arregles los 22 tests viejos de `Optimization/mobo/tests/ship/`; sólo añade los tuyos.

## C. Convenciones fijas del código nuevo

- Directorio `gaudi_source/pfa/`. Módulo CMake `SND_pfa`. Un `.cpp` por algoritmo.
- Headers compartidos (header-only): `SNDCellDecoder.h`, `SNDSamplingWeights.h`,
  `SNDTruthLink.h`.
- `SNDCellDecoder`: de `(nombre de colección, cellID)` devuelve
  `{det, station, layer, plane, view, layerOrdinal}` con `det` 0=SiTarget 1=SiPad 2=MTC,
  `view` 0=X/U 1=Y/V 2=2D, `layerOrdinal` contiguo dentro del detector (MTC:
  `station*MTC_NLayers + layer`, `MTC_NLayers` leído de las constantes DD4hep). Usa
  `dd4hep::DDSegmentation::BitFieldCoder` con los descriptores de `SND_compact.xml` (patrón:
  `EDM4HEP2RNTuple.cpp:119-129`).
- Colecciones (nombres fijos):

  ```
  MCParticles                                     MCParticleCollection        (Fase 0.2)
  SiTargetCaloHits SiPadCaloHits MTCSciFiCaloHits MTCScintCaloHits   CalorimeterHitCollection
  SiTargetCaloHitLinks … (4)                      CaloHitSimCaloHitLinkCollection
  SiTargetLayerClusters SiPadLayerClusters MTCSciFiLayerClusters MTCScintLayerClusters   ClusterCollection
  TruthClusters                                   ClusterCollection
  MuonCandidates                                  ReconstructedParticleCollection
  MuonCorridorClusters                            ClusterCollection
  MuonPID                                         ParticleIDCollection
  HadronicClusters                                ClusterCollection (Fase 4b)
  ReconstructedParticles                          ReconstructedParticleCollection
  RecoMCParticleLinks                             RecoMCParticleLinkCollection
  ```

- Codificación: `CalorimeterHit::type = det*10 + view`; `Cluster::type` = `det*10 + view`
  (capa), `100 + k` (3D), `200` (truth); `Cluster::shapeParameters` de capa =
  `[layerOrdinal, nHits, width_mm, zPlane]`; los reco hits conservan el **mismo orden** que
  los sim hits y un link 1:1.
- Pipelines: `gaudi_jobs/genie_pipeline/job6_pfa.py` (entrada `tracks_seq.edm4hep.root`,
  salida `pfa_seq.edm4hep.root`) y `job7_pfa_rntuple.py` (salida `ShipPFA.root`).

---

# FASE 0 — Prerrequisitos

## Paso 0.1 — Rama y registro

**Haz**
1. Pregunta al usuario si crear `feature/snd-pfa` desde `feature/genie-events`; créala.
2. Abre `docs/snd_pfa_progress.md`, pon la fecha en la primera fila y añade un enlace a
   ambos ficheros al final de `docs/gaudi_pipeline.md`.

**Comprueba** `git branch --show-current` = `feature/snd-pfa`.
**Anota** fecha de inicio.

## Paso 0.2 — `MCParticles` en la salida del overlay

Motivo: `events.edm4hep.root` y `tracks_seq.edm4hep.root` no tienen `MCParticles`; sin
ellas no hay E_ν, E_had ni muón verdaderos.

**Haz**
1. En `gaudi_source/EventOverlay.cpp`, en el bucle sobre fuentes (junto a los `accum`,
   ~l. 140), lee la `edm4hep::MCParticleCollection` `MCParticles` de cada fichero fuente.
2. Escribe en el frame de salida una colección `MCParticles`: con una sola fuente,
   `coll.clone()` (con relaciones); con N fuentes, copia campo a campo (PDG, generatorStatus,
   simulatorStatus, charge, time, mass, vertex, endpoint, momentum, momentumAtEndpoint)
   sin relaciones padre/hija, y deja un comentario diciendo que es una limitación.
3. Añade el parámetro de frame `MCParticleSourceOffsets` (`std::vector<int>`, índice inicial
   de cada fuente en `MCParticles`, en el orden de `SourceIDs`).
4. Añade la propiedad `WriteMCParticles` (bool, default `true`).
5. Reinstala. Desde `gaudi_jobs/genie_pipeline/` ejecuta `k4run job1_overlay.py`.
   **No** ejecutes job4 todavía (se hace una sola vez en 0.3).
6. Crea `Optimization/Analysis/pfa/check_truth.py`: recorre `events.edm4hep.root` y para
   cada evento comprueba (a) `len(MCParticles) > 0`; (b) para cada colección
   `*HitsWindowed`, `max(ContribParticleIDs) + offset < len(MCParticles)`; (c) el PDG de
   `MCParticles[id+offset]` coincide con `ContribPDGs` en todas las contribuciones; (d) hay
   un neutrino (|PDG| ∈ {12,14,16}); (e) en CC hay un leptón cargado primario
   (|PDG| ∈ {11,13,15}, generatorStatus 1). Imprime cuántos eventos tienen neutrino con
   `generatorStatus` 4 (beam) o cualquier otro valor.

**Comprueba** `podio-dump -e 0:3 events.edm4hep.root` lista `MCParticles`;
`check_truth.py` termina `OK` con (c) al 100 %.
**Anota** cómo aparece el neutrino en `MCParticles` (PDG, generatorStatus, si tiene
momento). Si ddsim NO lo conserva, escribe en el registro que E_ν,true se calculará como
`Σ E(primarias status 1) − 0.938 GeV` (aproximación) y úsalo así en las Fases 4–5.

## Paso 0.3 — Persistir q/p del ajuste de trazas

Motivo: el KF de `ACTSSequentialTracker` ajusta q/p en el campo del hierro
(`bv[Acts::eBoundQOverP]`, l. 1921, 2504) pero no lo escribe; el muón necesita su momento.

**Haz**
1. Lee `gaudi_source/ACTSSequentialTracker.cpp` l. 1855-1880 y 3070-3150. Confirma que el
   refit deja libre q/p (covarianza `m_seedVarQOverP`, refit 0.04) y que usa
   `SeqIronSlabBField`. Si no, para y pregunta.
2. En `writeTrack` (l. ~3080), tras recorrer los estados, añade dos `edm4hep::TrackState`
   más: `location = AtLastHit` (estado medido más aguas abajo) y `location = AtFirstHit`
   (más aguas arriba), con `referencePoint` = posición global (vía `toDD`), `phi`,
   `tanLambda`, **`omega = qOverP` [1/GeV]** del vector suavizado, y `covMatrix` con al
   menos los elementos diagonales de (phi, theta, qOverP) del suavizado. Deja los estados
   `AtIP`/`AtOther` exactamente como están.
3. Comprueba en `EDM4HEP2RNTuple.cpp` (l. 547, 567) y `event_display/event_display_eve.py`
   que ambos filtran por `location`; si el display itera todos los estados, filtra ahí
   `AtIP`/`AtOther`.
4. Reinstala. Desde `gaudi_jobs/genie_pipeline/` ejecuta `k4run job4_tracking_seq.py`
   (~21 min, un solo proceso; mide `cpu.pressure` antes y después) y después
   `k4run job5_rntuple_seq.py`.
5. Crea `Optimization/Analysis/pfa/check_qop.py`: para cada `GlobalTrack` con `type & 4`,
   identifica la MC particle dominante de sus hits MTC (misma lógica que
   `job4_purity.py`: medidas → cellID → `*ContribParticleIDs`); si es el muón primario,
   calcula `p_reco = 1/|qOverP|` del estado `AtLastHit` y `p_true` = |p| del muón en
   `MCParticles`. Imprime, por bins de `p_true` (2–5, 5–10, 10–20, >20 GeV): N, mediana y
   σ₆₈ de `(p_reco − p_true)/p_true`, fracción con signo de la carga correcto.

**Comprueba** `ShipHits.root` se sigue produciendo y el display abre un evento.
Criterio: |sesgo| < 15 % y signo correcto > 80 % para p > 5 GeV.
**Anota** la tabla de `check_qop.py`. Si falla el criterio, lee
`gaudi_jobs/muon_pipeline/MTC_OVERCURVATURE_FINDINGS.md` y pregunta antes de tocar el
tracker.

## Paso 0.4 — Esqueleto del módulo `SND_pfa`

**Haz**
1. Crea `gaudi_source/pfa/SNDCellDecoder.h` según §C. Constructor con
   `(const dd4hep::Detector&, const std::map<std::string,std::string>& readoutByCollection)`;
   método `decode(collectionName, cellID)`.
2. Crea `gaudi_source/pfa/SNDCaloHitBuilder.cpp` mínimo: propiedades `InputCollections`
   (4 nombres), `ReadoutNames` (4), lee las colecciones y cuenta hits por `(det, view,
   layerOrdinal)`; imprime el resumen en `finalize()`.
3. Añade a `CMakeLists.txt`:
   `gaudi_add_module(SND_pfa SOURCES gaudi_source/pfa/*.cpp LINK Gaudi::GaudiKernel
   k4FWCore::k4FWCore EDM4HEP::edm4hep podio::podioRootIO DD4hep::DDCore DD4hep::DDRec
   ROOT::ROOTNTuple)` (lista los ficheros explícitamente, como hace `SND_reco`).
4. Crea `gaudi_jobs/genie_pipeline/job6_pfa.py` con el patrón de `job4_tracking_seq.py`
   (`IOSvc.Input = ["tracks_seq.edm4hep.root"]`, `Output = "pfa_seq.edm4hep.root"`,
   `GeoSvc` con el compact, `SNDGeometry` para bitfields).
5. Crea `Optimization/Analysis/pfa/check_decoder.py`: para 1000 hits de cada colección,
   decodifica el cellID con `SNDGeometry.bitfields` en Python y compara con lo que imprime
   el algoritmo (activa un `Debug` que vuelque los 1000 primeros).

**Comprueba** compila e instala; `k4run job6_pfa.py -n 10` corre; `check_decoder.py` = `OK`
(0 discrepancias en layer/plane/station).
**Anota** fecha; nada más.

---

# FASE 1 — Hits calorimétricos calibrados: `SNDCaloHitBuilder`

## Paso 1.1 — Pesos de muestreo desde el material real

**Haz**
1. Crea `gaudi_source/pfa/SNDSamplingWeights.h` con: tabla `(dE/dx)_min` [MeV cm²/g] por
   nombre de material (`TungstenDens1910` 1.145, `Iron` 1.451, `Steel` 1.451, `Silicon`
   1.664, `Copper` 1.403, `Air` 1.815, `Scintillator` 1.956, `Carbon`/`CarbonFiber` 1.745,
   `Kapton`/`PCB`/`G10` 1.70, `Epoxy`/`Glue` 1.90; desconocido → 1.70 y se añade a una
   lista de avisos), y una función que, dada una lista de planos sensibles ordenados en z
   (con su espesor y material activo), recorra con
   `dd4hep::rec::MaterialManager(detector.world().volume()).materialsBetween(p0, p1)` el
   segmento (0,0,z_prev)→(0,0,z_plane) y devuelva por plano `{E_pasivo_MIP, E_activo_MIP,
   w = (E_pasivo + E_activo)/E_activo}`.
2. Obtén los planos sensibles al `initialize()` de `SNDCaloHitBuilder` recorriendo el TGeo
   (volúmenes con `SensitiveDetector` válido) o, si es más simple, desde
   `ISNDGeoSvc::surfacesOf(det)` + el espesor de los slices sensibles de las constantes
   (`SiTarget_SiThickness`, `SiPad_SiThickness`, `MTC_SciFiThickness`,
   `MTC_ScintThickness`; mira los nombres reales en `SND_compact.xml`). **Prohibido**
   escribir espesores numéricos en el código.
3. En `finalize()` imprime la tabla (det, ordinal, plane, z, E_pasivo, E_activo, w) y la
   lista de materiales desconocidos.

**EL PESO ES POR CAPA, no por plano** (corregido el 2026-09-23):
`w = (ΔE_MIP,pasivo,capa + Σ_planos ΔE_MIP,activo) / Σ_planos ΔE_MIP,activo`, el mismo para
todos los planos sensibles de la capa. Atribuir el absorbedor al primer plano que lo sigue
deja la segunda vista con peso ~1 e inútil para energía — medido en 50 eventos GENIE, el
plano V del MTC tenía más hits que U (2307 vs 2251) y aportaba el 0.45 % de la energía
frente al 55 % de U — y hace que una capa pierda casi toda su energía si el plano pesado
cae bajo umbral. Con el reparto por capa las vistas ortogonales coinciden al 0.6 % (U/V) y
1.5 % (X/Y): eso es además una comprobación física gratis, porque dos vistas del mismo paso
deben medir lo mismo.

**Comprueba** ejecuta `k4run job6_pfa.py -n 1` en baseline. Cifras medidas el 2026-09-23,
una por detector: **SiTarget** pasivo 7.656 MeV (3.5 mm W) / activo 0.2326 MeV (2 × 0.3 mm
Si) → **w = 33.91** en X y en Y; **SiPad** 22.50 / 0.2520 → **w = 92.25**; **MTC** pasivo
65.26 MeV (53 mm Fe antes de U + 3 mm antes del scint) / activo 4.582 MeV (0.2725 U +
0.2725 V + 4.037 scint) → **w = 14.96** en U, V y scint. La capa 0 de cada detector sale
distinta (90.27 SiPad, 15.24 MTC) porque tiene menos material delante. **Comprueba también
la periodicidad**: dentro de un detector, todas las capas salvo la 0 deben dar el MISMO w;
dos familias de valores significan que el barrido de material se está truncando (fue el bug
del 2026-09-23, paso nulo de TGeo en una frontera). Luego genera una variante con
`Optimization/Simulation/geometry/config.py` con `SiTarget_WThickness` doble y confirma que
w del SiTarget se duplica.
**Anota** la tabla baseline y la de la variante.

## Paso 1.2 — El algoritmo completo

**Haz**
1. Salidas: las 4 `*CaloHits` + 4 `*CaloHitLinks` (§C). Un reco hit por sim hit que supere
   `ThresholdMIP` (default 0.5), mismo orden.
2. Por hit: `E_MIP = E_dep / E_activo(plane)`; `energy = E_dep × w(capa) × GlobalScale`
   (`GlobalScaleECAL` para det 0/1, `GlobalScaleMTC` para det 2; defaults 1.0); `time` =
   mínimo de las contribuciones; `position` = `dd4hep::rec::CellIDPositionConverter`
   (nunca la posición del sim hit: en el MTC la sobrescribe `SND_SciFiAction`); `type` =
   `det*10 + view`; `energyError` = `E_MIP` (guardamos el valor en MIP ahí, documentado).
3. Propiedades: `ThresholdMIP`, `WeightMode` (`"sampling"`|`"none"`), `GlobalScaleECAL`,
   `GlobalScaleMTC`, `MIPEnergyOverride` (lista `"det:plane:MeV"`, vacía).
4. Añádelo a `job6_pfa.py`.

**Comprueba** `k4run job6_pfa.py -n 50`; `podio-dump` muestra las 8 colecciones; un script
rápido verifica `len(CaloHits) == len(Links)` y que cada link apunta al sim hit del mismo
cellID.
**Anota** hits/evento por colección antes y después del umbral.

## Paso 1.3 — Muestras de calibración (HTCondor)

**Haz**
1. En `simulation/run_script/steer/` crea steers mono-partícula sobre `multiplePG_base.py`
   (lista `particles` de un solo elemento), vértice (0, 0, −1640 mm), 500 eventos cada uno:
   `mu-` 10 GeV; `e-` 1, 2, 5, 10, 20 GeV; `pi+` 2, 5, 10, 20, 40 GeV. Nombra los ficheros
   `calib_<part>_<E>GeV`.
2. Lánzalos con `launch_multiplePG.sh` (añade los casos a su lista) a HTCondor; 11 jobs.
3. Para cada salida ejecuta `job1_overlay.py` y `job6_pfa.py` (sólo hit builder) en una
   carpeta `gaudi_jobs/calib_pipeline/` (copia de `genie_pipeline` con rutas
   parametrizadas por variable de entorno `INPUT_FILE`). Uno a uno.

**Comprueba** 11 ficheros `pfa_<part>_<E>.edm4hep.root` con 500 eventos.
**Anota** rutas de las 11 muestras.

## Paso 1.4 — Calibración

**Haz** crea `Optimization/Analysis/pfa/calib.py` que, sobre las muestras de 1.3:
1. μ 10 GeV: histograma de `E_MIP` (= `energyError`) por hit en SiTarget, SiPad, MTC
   scint, MTC SciFi. Imprime la moda.
2. e⁻: `E_vis = Σ energy` (SiTarget+SiPad) por evento; ajusta `E_vis = a·E`; imprime
   `a`, residuos por punto, σ/E por punto y el ajuste `σ/E = s/√E ⊕ c`.
3. π⁺: `E_vis,all = Σ energy` (4 colecciones); calcula `EoverH = 1/(pendiente/a)`; imprime
   residuos, σ/E, y la fracción media de energía en el MTC por E.
4. Escribe los valores elegidos: `GlobalScaleECAL = 1/a`, `GlobalScaleMTC` (el que
   linealiza los π con el MTC incluido; si no hay sensibilidad, 1.0), `EoverH`.

**Comprueba** moda MIP en 0.8–1.2; residuos EM < 5 % en 2–20 GeV; residuos π < 10 % en
5–40 GeV tras `EoverH`.
**Anota** `GlobalScaleECAL`, `GlobalScaleMTC`, `EoverH`, `s`, `c`, fuga al MTC vs E. Pon
las escalas como defaults en `job6_pfa.py`.

## Paso 1.5 — Cota superior de E_had en GENIE

**Haz** crea `Optimization/Analysis/pfa/ehad_bound.py` sobre `pfa_seq.edm4hep.root`
(regenéralo con las escalas de 1.4): por evento, `E_had,true = Σ E_cin` de las primarias de
estado final que no sean leptones (|PDG| ∉ {11,12,13,14,15,16}) y `E_vis,nomu = Σ energy`
de los reco hits cuya MC particle dominante (link → índice sim → `*ContribParticleIDs` +
offset) no es el muón primario. Imprime correlación, sesgo mediano de
`(EoverH·E_vis,nomu − E_had,true)/E_had,true` y σ₆₈ por bins (0–2, 2–5, 5–10, 10–20,
>20 GeV).

**Comprueba** correlación > 0.9 y |sesgo| < 10 %.
**Anota** la tabla; es el techo de la Fase 4.

---

# FASE 2 — Clustering por capa

## Paso 2.1 — `SNDLayerClusterer`

**Haz**
1. Crea `gaudi_source/pfa/SNDLayerClusterer.cpp`. Entrada: **una** colección `*CaloHits`;
   salida: **una** `ClusterCollection`. En `job6_pfa.py` se instancia 4 veces
   (`LayerClusterer_SiTarget`, `_SiPad`, `_MTCSciFi`, `_MTCScint`) escribiendo
   `SiTargetLayerClusters`, `SiPadLayerClusters`, `MTCSciFiLayerClusters`,
   `MTCScintLayerClusters` (los nombres de §C).
2. Algoritmo: agrupa hits por `(layerOrdinal, view)`; ordena por índice de canal (strip /
   x+y de pad / x+y de tile, decodificado con `SNDCellDecoder`); hits contiguos (hueco ≤
   `MaxGapChannels`: 1 strips, 0 pads/tiles — en 2D "contiguo" = vecino en 8-conectividad)
   forman un proto-cluster; máximos locales = hits con E > ambos vecinos y E >
   `SeedMIP` (2.0); si un proto-cluster tiene ≥2 máximos separados por un mínimo <
   `ValleyFraction` (0.5) × el menor de los dos máximos, divídelo asignando cada hit al
   máximo más cercano (energía de los hits del valle repartida linealmente por distancia).
3. Por cluster: `energy` = Σ, `position` = baricentro ponderado, `type = det*10+view`,
   `shapeParameters = [layerOrdinal, nHits, width_mm, zPlane]`, `addToHits`.
4. Propiedades: `MaxGapChannels`, `SeedMIP`, `ValleyFraction`, `MinClusterMIP` (0.5).

**Comprueba** `k4run job6_pfa.py -n 100`; script `layer_clusters_check.py`: para cada
cluster, Σ E(hits) == energy y todos sus hits comparten `(layerOrdinal, view)`.
**Anota** clusters/evento y `nHits` mediano por colección.

## Paso 2.2 — `SNDTruthClusterer` y `SNDTruthLink.h`

**Haz**
1. `SNDTruthLink.h`: clase que, construida con el frame (los 4 `*CaloHitLinks`, los 4
   `*ContribParticleIDs`/`*ContribPDGs`, `MCParticleSourceOffsets`, las
   `*HitsWindowed_Contributions`), devuelve para un reco hit `(mcIndex dominante,
   fracción de energía)` y para un cluster la MC particle dominante y su pureza (fracción
   de energía). Reproduce la lógica de atribución de `EDM4HEP2RNTuple.cpp:407-456`.
2. `SNDTruthClusterer.cpp`: un `Cluster` por MC particle con reco hits (todas las
   colecciones), `type = 200`, `shapeParameters[0] = mcIndex`, `addToHits`. Salida
   `TruthClusters`.

**Comprueba** Σ energía de `TruthClusters` == Σ energía de todos los `*CaloHits` (±1e-6).
**Anota** nada más que la fecha.

## Paso 2.3 — Validación del clustering de capa (GENIE, 2000 eventos)

**Haz** crea `Optimization/Analysis/pfa/layer_perf.py`:
1. Por colección: clusters/evento, distribución de `nHits` y `width_mm`.
2. SiPad: por capa, nº de clusters vs nº de MC particles distintas con hits (vía
   `SNDTruthLink`); razón mediana.
3. Pureza (fracción de energía de la MC particle dominante) por colección: mediana y
   fracción con pureza = 1.
4. **Aislamiento del muón**: en eventos ν_μ CC, para el muón primario, fracción de capas
   que atraviesa en las que sus hits forman un cluster propio con E ≤ 3 MIP y ≤ 3 canales
   (pureza > 0.8), por detector y separando capas "dentro de la mancha hadrónica" (hay
   otro cluster a < 20 mm) y "fuera".
5. Añade a `event_display/event_display_eve.py` una opción `--draw-clusters` que dibuje
   los `*LayerClusters` de `pfa_seq.edm4hep.root` como segmentos/cajas.

**Comprueba** pureza mediana ≥ 0.9 en SiPad y MTC scint; display en 3 eventos coherente
con los hits.
**Anota** todas las tablas, especialmente el aislamiento del muón (de él depende la Fase 3).

---

# FASE 3 — Muón: `SNDMuonFinder` (prioridad 1)

## Paso 3.1 — Corredor traza–calorímetro

**Haz**
1. Crea `gaudi_source/pfa/SNDMuonFinder.cpp`. Entrada: `GlobalTracks`, `MTCTracks`, las
   colecciones `*LayerClusters`. Salidas: `MuonCandidates`, `MuonCorridorClusters`,
   `MuonPID`.
2. Para cada traza (primero `GlobalTracks` con `type & 4`, luego `MTCTracks` no usadas):
   construye la trayectoria como polilínea por los `referencePoint` de los estados
   `AtOther` (ordenados en z); extrapola en línea recta antes del primero y después del
   último hasta `ExtrapolateLayers` (5) capas.
3. Para cada capa (ordinal) entre los extremos extrapolados, calcula la posición esperada
   y asigna al corredor los clusters de capa con distancia en la coordenada de su vista
   (x para X/U, y para Y/V, 2D para pads/tiles; para U/V del MTC usa la coordenada de la
   vista rotada ±5°, igual que `MTCSciFiMeasConverter`) < `CorridorHalfWidth[det]`
   (SiTarget 2, SiPad 10, MTC SciFi 15, MTC scint 40 mm) **y** E ≤ `CorridorMaxMIP` (3.0)
   MIP. Un cluster sólo puede pertenecer a un corredor (la traza con más medidas gana).
   Cuenta aparte los clusters cercanos que no son MIP-like.
4. Calcula por traza: `nMTCLayersHit`, `nLayersTraversed` (capas MTC entre el primer y el
   último hit del corredor), `fracLayersHit`, `meanCorridorMIP`, `maxCorridorMIP`,
   `ironRange_mm` (Σ del pasivo de hierro de las capas MTC atravesadas, de la tabla de
   pesos de 1.1 convertida a mm de Fe), `chi2ndf` de la traza.

**Comprueba** `k4run job6_pfa.py -n 100`; script: cada cluster de `MuonCorridorClusters`
aparece en un solo candidato.
**Anota** distribución de `nMTCLayersHit` y `meanCorridorMIP` para trazas con muón MC
dominante vs sin él (usa `SNDTruthLink` sobre los clusters del corredor).

## Paso 3.2 — Identificación y momento

**Haz**
1. Candidato = muón si: `type & 4`; `nMTCLayersHit ≥ MinMTCLayers` (8);
   `fracLayersHit ≥ MinFracLayers` (0.6); `meanCorridorMIP ≤ MaxMeanMIP` (2.0). Todo
   propiedades.
2. Momento — **el rango es la medida principal** (decidido el 2026-09-22 a la vista de
   los números del paso 0.3: la curvatura en el MTC da σ₆₈ ≈ 40–50 % para p > 10 GeV,
   dominada por dispersión múltiple en 2.2 m de Fe). Por defecto
   `p = p_range = RangeSlopeGeVPerMm (0.0115) × ironRange_mm` con `p_source = 1`. Se usa
   la curvatura (`p_spec = 1/|qOverP|` del estado `AtLastHit`, `p_source = 0`) sólo
   cuando el muón **sale** del detector — último hit del corredor en la última capa del
   MTC60, donde el rango sólo da una cota inferior. `charge` = signo de `qOverP` siempre
   (el signo sí es fiable: 88–97 % correcto por bins), y se propaga
   `relVar = cov(qOverP)/qOverP²` a `MuonPID`. Dirección = la del último estado (`phi`,
   `tanLambda`).
   - Ojo con `relVar`: las trazas sin campo (SiTarget/SiPad puras) salen con el valor de
     la semilla, `q/p = −1/3` y `relVar = 0.36` exactos; ése no es un momento medido.
     Umbral sugerido `MaxRelQopVar = 0.30` (las medidas de verdad dan ~0.026).
   - El fit **subestima su propio error**: declara relVar 0.026 (16 %) frente a una
     dispersión real del 40–50 %. Si alguna vez se pondera `p_spec`, hay que recalibrar
     ese error con la dispersión medida, no con la covarianza del fit.
3. Rellena `MuonCandidates` (`ReconstructedParticle`: PDG ±13, energy, momentum, charge,
   `addToTracks`, `addToClusters` con los del corredor) y `MuonPID` (`ParticleID` con
   `PDG`, `likelihood = fracLayersHit`, `parameters = [nMTCLayersHit, nLayersTraversed,
   fracLayersHit, meanCorridorMIP, maxCorridorMIP, ironRange_mm, chi2ndf, p_source]`).
4. Propiedad `UseTruthMuon` (bool, default false): en vez de trazas, construye un candidato
   a partir del `TruthCluster` del muón primario (momento = MC) y aplica los mismos cortes
   de ID sobre sus clusters de capa → cota superior.

**Comprueba** `k4run job6_pfa.py` completo (2000 eventos).
**Anota** nº de candidatos/evento.

## Paso 3.3 — Validación del muón (GENIE, 2000 eventos)

**Haz** crea `Optimization/Analysis/pfa/muon_perf.py`:
1. Clasifica eventos por `MCParticles`: ν_μ CC, ν_e CC, ν_τ CC, NC.
2. **Eficiencia**: fracción de ν_μ CC con un `MuonCandidate` cuyo corredor está dominado
   por el muón primario; por clases: todos / muón con ≥20 hits / muón que llega al MTC.
   Compara con "el muón tiene una `GlobalTrack`" (referencia del informe de tracking:
   66.2 % y 85.5 % con ≥20 hits).
3. **Fakes**: fracción de NC y ν_e CC con algún `MuonCandidate`; lista qué MC particle
   domina esos corredores (PDG) y su E.
4. **Momento**: por bins de `p_true`, sesgo mediano y σ₆₈ de `(p − p_true)/p_true`;
   fracción con carga correcta; separa `p_source` 0/1. La referencia de partida es el
   paso 0.3: la curvatura da σ₆₈ 0.39–0.98 frente al momento en el endpoint. El objetivo
   del rango es mejorarlo claramente para los muones que se paran; si no lo consigue, hay
   que decirlo.
5. Cota: repite 2 con `UseTruthMuon=true` (segunda pasada de `job6_pfa.py` a otro fichero).
6. Barrido: `MinMTCLayers` ∈ {6, 8, 10} × `CorridorHalfWidth` MTC SciFi ∈ {10, 15, 25} mm
   (9 pasadas de `job6_pfa.py`, una a una, segundos cada una); tabla eficiencia/fakes.

**Comprueba** eficiencia ≥ (eficiencia de traza del muón − 3 puntos); fakes ≤ 5 % en NC.
**Anota** todas las tablas y los defaults elegidos (fíjalos en `job6_pfa.py`).

---

# FASE 4 — Sistema hadrónico y E_ν (prioridad 2)

## Paso 4.1 — `SNDPFOBuilder`: E_had "todo menos el muón"

**Haz**
1. Crea `gaudi_source/pfa/SNDPFOBuilder.cpp`. Entrada: `*CaloHits`, `*LayerClusters`,
   `MuonCandidates`, `MuonCorridorClusters`, `SiTargetTracks`, `MCParticles` (sólo para los
   links). Salidas: `ReconstructedParticles`, `RecoMCParticleLinks`, parámetros de frame
   `pfa_E_nu, pfa_E_had, pfa_E_had_ecal, pfa_E_had_mtc, pfa_p_mu, pfa_q_mu, pfa_p_mu_source,
   pfa_flavour, pfa_vtx_x, pfa_vtx_y, pfa_vtx_z` (float/int).
2. **Vértice** (propiedad `VertexMode`): `"clusters"` = cluster de capa más aguas arriba
   con E > `VertexSeedMIP` (2) seguido de actividad en ≥2 de las 3 capas siguientes;
   `(x,y)` de los clusters X e Y de ese ordinal, z del plano. `"track"` = primer estado de
   la `SiTargetTrack` más temprana. Implementa ambos.
3. `E_had,ECAL` = Σ energy de los clusters de capa de SiTarget+SiPad no asignados a un
   corredor de muón; `E_had,MTC` = idem MTC. **Corrección de corredor**: para cada cluster
   del corredor del muón con E_MIP > 1, suma al hadrónico `(E_MIP − 1) × E_activo × w`.
   `E_had = EoverH × (E_had,ECAL + E_had,MTC)` (propiedad `EoverH` con el valor de 1.4).
4. Muón: si hay `MuonCandidates`, toma el de más `nMTCLayersHit`; `E_ν = p_μ + E_had`,
   `pfa_flavour = 14`; si no, `E_ν = E_had`, `pfa_flavour = 0`.
5. `ReconstructedParticles`: copia del PFO muón (con `tracks`/`clusters`) y un PFO
   hadrónico (`type = 1`, `PDG = 0`, `energy = E_had`, `momentum` = E_had × dirección
   vértice→baricentro de los clusters hadrónicos, `clusters` = los no usados).
   `RecoMCParticleLinks`: muón → MC muón dominante (peso = pureza del corredor);
   hadrónico → neutrino (peso 1).

**Comprueba** `k4run job6_pfa.py` completo; `podio-dump` muestra los parámetros `pfa_*`.
**Anota** fecha.

## Paso 4.2 — Validación de E_had, E_ν y vértice (GENIE, 2000 eventos)

**Haz** crea `Optimization/Analysis/pfa/enu_perf.py`:
1. `E_had,reco` vs `E_had,true` (definición del paso 1.5: Σ E_cin de los hadrones
   primarios de estado final): sesgo mediano y σ₆₈/E por bins (0–2, 2–5, 5–10, 10–20,
   >20 GeV) **y por bins de z del vértice**, con la cota del paso 1.5 al lado: 0.949 /
   0.888 / 0.790 / 0.728 / 0.638 de z = −2 m al final del MTC. La diferencia con esa cota
   es lo que cuestan corredor y umbrales; el resto es fuga longitudinal, que ningún
   algoritmo recupera.
2. ν_μ CC: `E_ν,reco` vs `E_ν,true` (del neutrino en `MCParticles`, o la aproximación
   anotada en 0.2): sesgo y σ₆₈/E por bins; separado por `pfa_p_mu_source`.
3. Vértice: Δx, Δy, Δz frente al vértice MC (vértice de la primaria) para ambos
   `VertexMode`; mediana y fracción con |Δz| > 3 capas.
4. Fuga: fracción de eventos con `E_had,MTC/E_had > 0.3` por bin de E.

**Comprueba** |sesgo| de E_ν < 10 % en ν_μ CC; diferencia de σ/E respecto a la cota de 1.5
anotada.
**Anota** todas las tablas; elige `VertexMode` y fíjalo. **Para aquí y presenta al
usuario** los resultados junto con la pregunta de si hace falta 4.3.

## Paso 4.3 — `SNDHadronicClusterer` (clusters 3D) — sólo si el usuario lo aprueba tras 4.2

**Haz**
1. Crea `gaudi_source/pfa/SNDHadronicClusterer.cpp`; entrada: clusters de capa no usados
   por el muón; salida `HadronicClusters` (`type = 100+k`, `addToClusters` de los de
   capa, `energy`, `position`, eje en `shapeParameters`).
2. Enlace cono capa a capa (`LongiConeLinking` de CyberPFA): ordena por `layerOrdinal`;
   semilla = cluster de mayor E no asignado; eje = vértice → semilla; añade capa a capa
   (hasta `MaxLayerSkip` = 2 capas vacías) los clusters cuyo ángulo respecto al eje <
   `ConeAngle` (0.3 rad); recalcula el eje con el baricentro; repite hasta agotar
   semillas con E > `SeedMIP` (5). Empareja X e Y (y U/V) del mismo ordinal exigiendo
   `|E_X − E_Y|/(E_X + E_Y) < 0.5`.
3. Cambia `SNDPFOBuilder` para que, si existe `HadronicClusters`, haga un PFO hadrónico
   por cluster 3D y calcule `E_had` sólo con los clusters 3D unidos al vértice (distancia
   del eje al vértice < `VertexLinkMM` = 30).

**Comprueba** nº de `HadronicClusters`/evento; pureza por `SNDTruthLink`; repite 4.2 y
compara.
**Anota** la comparación 4.1 vs 4.3.

---

# FASE 5 — Exportación y métricas MOBO

## Paso 5.1 — `PFASummaryWriter` → `ShipPFA.root`

**Haz**
1. Crea `gaudi_source/pfa/PFASummaryWriter.cpp` con el patrón RNTuple de
   `EDM4HEP2RNTuple.cpp:119-140`. Salida `ShipPFA.root` con:
   - `EventPFA`: `window_id, flavour_true, is_cc, E_nu_true, E_lep_true, E_had_true,
     vtx_true_x/y/z, n_muon_cand, muon_matched, p_mu_reco, q_mu_reco, p_mu_true,
     p_mu_source, E_had_reco_ecal, E_had_reco_mtc, E_had_reco, E_nu_reco, flavour_reco,
     vtx_reco_x/y/z`;
   - `PFOs`: `window_id, type, pdg, E, px, py, pz, n_clusters, n_tracks, mc_index, mc_pdg,
     mc_E`;
   - `LayerClusters` (si `WriteClusters`): `window_id, det, view, layer, E, x, y, z, nhits,
     width, used_by` (0 libre / 1 muón / 2 hadrónico).
   Las columnas `*_true` y `muon_matched` se calculan con `SNDTruthLink` aquí (no en los
   algoritmos de reco).
2. Crea `gaudi_jobs/genie_pipeline/job7_pfa_rntuple.py` (entrada `pfa_seq.edm4hep.root`).
3. Añade a `genie_pipeline_seq.sh` las líneas `k4run job6_pfa.py` y
   `k4run job7_pfa_rntuple.py`.

**Comprueba** `ROOT.RNTupleReader.Open("EventPFA", "ShipPFA.root").GetNEntries()` = 2000.
**Anota** fecha.

## Paso 5.2 — Métricas en el MOBO

**Haz**
1. En `Optimization/mobo/src/mobo/ship/metrics.py` añade `pfa_metrics(path) -> dict`
   leyendo `EventPFA`: `eff_muon_id` (ν_μ CC con `muon_matched` / ν_μ CC), `fake_muon_nc`
   (NC con `n_muon_cand > 0` / NC), `res_pmu` (σ₆₈ de `(p_mu_reco − p_mu_true)/p_mu_true`
   en emparejados), `bias_ehad`, `res_ehad`, `res_enu` (ν_μ CC), `frac_leak_mtc`,
   `nevents_pfa`. Ausencia del fichero → no es error (mismo patrón que `track_metrics`).
2. Copia `job6_pfa.py`/`job7_pfa_rntuple.py` a `Optimization/Analysis/` con las rutas
   parametrizadas como sus vecinos; en `payload.py` (l. ~405-448) añádelos a la cadena tras
   `job5_rntuple_seq` cuando el evaluator tenga `pfa_metrics: true`, y mezcla el dict de
   `pfa_metrics` en las métricas del trial.
3. Añade `pfa_metrics: true` a `conf/evaluator/snd_genie.yaml` **sólo** en un nuevo
   `conf/experiment/snd_genie_pfa.yaml` (copia de `snd_genie.yaml` con objetivos
   `eff_muon_id` max, `res_enu` min, `cost_proxy` min; `ref_from_baseline` como en el
   original). No toques `snd_genie.yaml`.
4. Crea `Optimization/mobo/tests/ship/test_pfa_metrics.py` que escribe un `EventPFA` de
   20 filas con `ROOT.RNTupleWriter` en `tmp_path` y comprueba las 8 claves y dos valores
   calculados a mano.

**Comprueba** `pytest Optimization/mobo/tests/ship/test_pfa_metrics.py` pasa;
`pfa_metrics('gaudi_jobs/genie_pipeline/ShipPFA.root')` reproduce los números de 3.3 y 4.2.
**Anota** el dict de métricas del baseline.

## Paso 5.3 — Un trial local de humo

**Haz** con el venv del MOBO activado después de key4hep, `nice -n 19 timeout 2h mobo-run
experiment=snd_genie_pfa executor=local max_trials=1` (un solo proceso; mide
`cpu.pressure` antes y después). Cualquier campaña real: HTCondor, driver incluido
(`CLAUDE.md`).

**Comprueba** el trial termina y su `metrics.json` contiene las claves `pfa_*`.
**Anota** tiempo del trial y las métricas.

---

# FASE 6 (opcional; sólo si el usuario la pide) — Cascadas EM: ν_e CC y π⁰

**Haz** `gaudi_source/pfa/SNDEMShowerFinder.cpp` sobre los clusters de capa no usados:
semillas = máximos locales densos (E > 5 MIP en ≤ 3 canales) en las 10 capas siguientes al
vértice; enlace recto capa a capa desde el vértice (el `HoughClusteringAlg` de CyberPFA
reducido a rectas desde un punto conocido); ajuste del perfil longitudinal (Γ, `t_max ≈
ln(E/E_c) − 0.5` en X₀ de W) y lateral (R_M(W) ≈ 9 mm); e vs γ por presencia de una
`SiTargetTrack` que empiece en el vértice. `pfa_flavour = 12` si hay cascada EM ≥
`MinEMEnergy` unida al vértice y no hay muón. Validación: eficiencia/pureza ν_e CC vs NC
con π⁰; resolución de E_e con las guns de e⁻ del paso 1.3; muestra GENIE filtrada
`genie_gst_to_hepmc.py --flavour 12`.

---

# Cierre — verificación end-to-end

Desde `gaudi_jobs/genie_pipeline/`:
```
k4run job1_overlay.py            # MCParticles presente
k4run job4_tracking_seq.py       # q/p en AtLastHit (sólo si el tracker cambió)
k4run job6_pfa.py                # hits → clusters → muón → PFOs
k4run job7_pfa_rntuple.py        # ShipPFA.root
python ../../Optimization/Analysis/pfa/check_truth.py
python ../../Optimization/Analysis/pfa/check_qop.py
python ../../Optimization/Analysis/pfa/muon_perf.py
python ../../Optimization/Analysis/pfa/enu_perf.py
python -c "from mobo.ship.metrics import pfa_metrics; print(pfa_metrics('ShipPFA.root'))"
```
Escribe `docs/informe_pfa_<fecha>.md` al estilo de `docs/informe_tracking_mobo_2026-09-11.md`
con: tabla final de métricas y sus cotas truth, defaults elegidos, lo que no funcionó y los
pendientes.

---

# Anexo 1 — Hechos del repo en los que se apoya cada orden

| Hecho | Dónde |
|---|---|
| Los tres subdetectores son `sens.setType("calorimeter")` → `edm4hep::SimCalorimeterHit` | `simulation/detector_plugin/{SiTargetDetector.cpp:18, SiPadDetector.cpp:131, MTCDetector.cpp:63}` |
| SiTarget: 120 capas × (3.5 mm W + Si 0.3 mm plano X + Si 0.3 mm plano Y), strips 75.5 µm; `id: system:8,layer:8,slice:4,plane:1,column:2,row:2,strip:14` | `simulation/geometry/SND_compact.xml:527-544`, `SiTarget.xml` |
| SiPad: 22 capas × (10 mm W + Si 0.65 mm), pads 2D 5.53 mm; `id: system:8,layer:8,slice:4,x:9,y:9,wafer:6` | `SND_compact.xml:174-175, 313-497` |
| MTC: 3 estaciones × 14 capas × (53 mm Fe + SciFi U + SciFi V 1.35 mm + 3 mm Fe + scint 20 mm); SciFi strips 1D estéreo ±5° de 1 mm (`plane` 0/1), scint 50×50 mm (`plane` 2); `id: system:8,station:2,layer:8,slice:4,plane:2,strip:14,x:9,y:9`; −1.9 T sólo en el hierro | `SND_compact.xml:248-278, 499-525`, `MTCDetector.xml:18-41`, `MTCDetector.cpp` |
| `EventOverlay` escribe `*HitsWindowed`, `*_Contributions` y los parámetros de frame `*SourceIDs`, `*ContribPDGs`, `*ContribParticleIDs` (índice de MCParticle dentro de su fuente); divide el MTC por `plane` (bit 22) | `gaudi_source/EventOverlay.cpp:158-177, 256-271, 296-305` |
| `events.edm4hep.root` y `tracks_seq.edm4hep.root` no tienen `MCParticles`; `tracks_seq` sí tiene hits + medidas + trazas | comprobado 2026-09-21 con `podio.root_io.Reader` |
| `ACTSSequentialTracker`: estados `AtIP` y `AtOther` (`omega` = inclinación estéreo); q/p en el ajuste pero no persistido; `Track::type` 1=SiT 2=SiP 4=MTC | `gaudi_source/ACTSSequentialTracker.cpp:450-457, 1863-1874, 1921, 2504, 3070-3150` |
| `EDM4HEP2RNTuple` y el display filtran estados por `location` | `EDM4HEP2RNTuple.cpp:547, 567` |
| Scint del MTC sólo como veto geométrico (30 mm) | `MTCSciFiMeasConverter.cpp:200-213, 369-390` |
| Materiales hard-coded en el tracker (no imitar) | `ACTSGeoSvc.cpp:370-400` |
| `ISNDGeoSvc`: `surfacesOf`, `layerOrdinalOf`, `nLayersOf`, `surfaceByAddress` | `gaudi_source/ISNDGeoSvc.h` |
| `SND_SciFiAction` sobrescribe `hit.position.y` en el MTC | `simulation/ddg4/SND_SciFiAction.cpp:30-113` |
| Atribución hit → MC particle dominante (patrón a reutilizar) | `EDM4HEP2RNTuple.cpp:407-456`, `gaudi_jobs/genie_pipeline/job4_purity.py` |
| Métricas MOBO desde RNTuples; `eff_track`; objetivos | `Optimization/mobo/src/mobo/ship/metrics.py:33-43, 205-260`, `payload.py:405-448`, `conf/experiment/snd_genie.yaml:91-96` |
| Guns: `multiplePG_base.py` + `steer/`, `launch_multiplePG.sh`, `generic_condor_multiplePG.sh` | `simulation/run_script/README_multiplePG.md` |
| GENIE→HepMC: neutrino como *beam*, leptón + hadrones status 1 | `simulation/run_script/genie_gst_to_hepmc.py:24-27` |
| En el stack: `edm4hep/CaloHitSimCaloHitLinkCollection.h`, `RecoMCParticleLinkCollection.h`, `ParticleIDCollection.h`; `TrackState::AtFirstHit=2, AtLastHit=3`; `MutableCluster::addToClusters/addToHits`; `DDRec/MaterialManager.h::materialsBetween` | `/cvmfs/sw.hsf.org/key4hep/releases/2026-02-01/x86_64-almalinux9-gcc14.2.0-opt/` |
| Build: un módulo `SND_reco`; sin ctest; pytest en `Optimization/mobo/tests/` (22 en `ship/` fallan) | `CMakeLists.txt:37-55` |

# Anexo 2 — Mapa CyberPFA → SND-PFA

| CyberPFA (`RecPFACyber/src/…`) | SND-PFA | Paso |
|---|---|---|
| `Tools/CaloHitsCreator` | `SNDCaloHitBuilder` | 1.2 |
| `CaloHalfCluster`, `Calo1DCluster`, `LocalMaxFindingAlg`, `EnergySplittingAlg` | `SNDLayerClusterer` | 2.1 |
| `Truth*Alg` | `SNDTruthClusterer`, `UseTruthMuon` | 2.2, 3.2 |
| `TrackExtrapolatingAlg`, `TrackMatchingAlg`, `TrackClusterConnectingAlg` | `SNDMuonFinder` | 3.1–3.2 |
| `PFOCreatingAlg`, `OutputCreator` | `SNDPFOBuilder` | 4.1 |
| `HcalClusteringAlg::LongiConeLinking`, `ConeClusteringAlg`, `AxisMergingAlg` | `SNDHadronicClusterer` | 4.3 |
| ntuple de depuración | `PFASummaryWriter` | 5.1 |
| `HoughClusteringAlg` (fotones) | `SNDEMShowerFinder` | Fase 6 |

# Anexo 3 — Riesgos

- `materialsBetween` a lo largo de (0,0) cruza el hierro magnetizado central del MTC:
  correcto para este diseño; si la geometría tuviera huecos en el eje, muestrear varios
  (x,y).
- `CellIDPositionConverter` en el MTC estéreo devuelve y=0: los clusters U/V sólo tienen
  su coordenada de vista; no intentar 2D en el MTC salvo vía scint o el emparejamiento U/V
  que ya hace `MTCSciFiMeasConverter` (no duplicar).
- Si el MOBO lleva el W a 2–8 mm, la linealidad de 1.4 debe revisarse en dos variantes
  extremas (dos ddsim en Condor) antes de fiarse de `res_enu` como objetivo.
- Si una etapa reordena o filtra hits, el índice al parámetro de frame se pierde: ir
  siempre por `CaloHitSimCaloHitLink → índice sim`.

# Referencias

- CyberPFA: `Reconstruction/RecPFACyber` en https://code.ihep.ac.cn/cepc/CEPCSW
- Diseño del ECAL de cristal del CEPC: arXiv:2602.21489, arXiv:2602.09836.
- F. Guo, "Particle Flow in CEPC: from CDR to Ref-det TDR" (FCC week 2025).
