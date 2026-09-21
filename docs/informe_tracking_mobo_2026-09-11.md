# Tracking secuencial y optimización MOBO sobre eventos GENIE — informe completo

Estado a 2026-09-11, 16:30 (segunda edición del día: incorpora el bug de signo del marco
ACTS, el retuning del MTC y la limpieza del pool que se hicieron por la tarde). Muestra: 2000
interacciones de neutrino GENIE (todos los sabores, CC+NC), re-simuladas sobre la geometría
baseline corregida de **120 capas de SiTarget**, que es la geometría sobre la que se
generaron los vértices.

Nomenclatura: **local track** = traza ajustada dentro de un solo detector (colecciones
`SiTargetTracks` / `SiPadTracks` / `MTCTracks`); **global track** = objeto final en
`GlobalTracks`, o bien *spliced* (locales de detectores vecinos empalmadas y reajustadas) o
bien *promoted* (una local sin vecino con el que empalmar).

---

## Índice

1. Resumen ejecutivo
2. La geometría: qué estaba mal y qué se corrigió
3. Tracking — el problema del seeder y su rediseño
4. Tracking — resultados
5. Tracking — pureza y fragmentación (verdad MC)
6. Tracking — los barridos de cortes
7. MOBO — infraestructura y espacio de búsqueda
8. MOBO — resultados: Sobol, qLogNEHVI, frente de Pareto
9. El sesgo del muestreo de vértices
10. Tiempos
11. Lo que salió mal en la infraestructura
12. Estado del código y pendientes

---

## 1. Resumen ejecutivo

- **La geometría del repo no era la del MC.** El compact generaba 166 capas de SiTarget;
  los vértices GENIE se muestrearon sobre 120. Se corrigió, se fijó el techo de capas de
  forma permanente y se re-simuló la muestra entera.
- **El seeder del tracker solo podía encontrar trazas paralelas al haz.** Se reescribió
  (Hough de tracklets con dirección).
- **Había un bug de signo en el marco local de ACTS** que dejaba las medidas reflejadas en
  x respecto a las semillas. Era invisible con el seeder viejo (trazas junto al eje) y fatal
  con el nuevo: el muón del evento 1, a x = −138 mm, tenía su semilla a 274 mm de sus
  propios hits. Corregido en el calibrador, en la inversión estéreo U/V del MTC y en la
  escritura de los estados (que también los dibujaba reflejados).
- **Con el marco correcto se re-barrieron los cortes del MTC** con una métrica nueva, la
  **eficiencia de muón** (el muón MC deja hits → sale una global dominada por él): mínimo
  de 3 → 5 medidas y puertas χ² 5/10 → 8/16 cuestan 0.8 puntos de muón y quitan la cola de
  fragmentos.
- **Se limpió el pool del MTC**: aceptancia del plano, confirmación por celda de
  scintillator y purga de densidad con borde. Pool 1086 → 272 medidas/evento, eficiencia
  de muón 65.0 → **66.2 %**, eventos con >5 globales 6.2 → **0.2 %**.
- **Cifras finales** (2000 eventos): 48.6 % de eventos con ≥1 global, **0.67 globales por
  evento**, 36.5 % con exactamente una; **eficiencia de muón 66.2 %, 85.5 % para muones con
  ≥20 hits**. Pureza mediana 1.00 en SiPad y MTC; el MTC pasa de 2.69 a **1.84 locales por
  partícula** y de 56 % a **17 %** de locales mezcladas.
- **MOBO**: 100 trials sobre 8 dimensiones, 96 completados, en **3 h 45 min** de reloj.
  Frente de Pareto de 16 puntos que **domina al baseline** en los tres objetivos. La campaña
  corrió con el tracking de la mañana (con el bug de signo); la conclusión sobre la
  geometría no cambia, la escala de `eff_track` sí.
- **Pero 14 de los 16 puntos del frente están en el mínimo de capas de SiTarget (20).**
  Es el sesgo del muestreo de vértices: un detector corto se acredita interacciones que su
  material no habría producido. El frente es una cota superior, no un diseño.

---

## 2. La geometría: qué estaba mal y qué se corrigió

### 2.1 El diagnóstico

Tres sitios tenían el número de capas del SiTarget y solo dos coincidían:

| fuente | capas SiTarget |
|---|---|
| `SND_compact.xml` (el que usa ddsim y la reconstrucción) | **166** |
| `SND_compact.gdml` (export del 28-jul, el que usó GENIE) | 120 |
| baseline de MOBO | 120 |

El compact calculaba `SiTarget_NLayers = floor(SiTarget_dim_z / SiTarget_spacing)` =
floor(1830/11) = 166. GENIE, en cambio, muestreó los vértices sobre el detector del GDML.

**Verificación en los datos**: el histograma de vértices GENIE en la región del SiTarget es
denso y uniforme (~35 por 100 mm) hasta z ≈ −1690 mm y **se corta en seco ahí**. Con 166
capas el material llegaría a −2196 y habría ~180 vértices más en esa franja: hay **2** de
487. Y −1690 es exactamente 120 × 11 mm contando desde −370.

Todo lo reconstruido hasta ese momento propagaba vértices generados sobre un detector de
120 capas a través de uno de 166, con el SiTarget extendido 506 mm más aguas arriba de donde
el generador creía. La muestra ya simulada tenía índices de capa hasta 165: había que
re-simular, no solo re-reconstruir.

### 2.2 El arreglo

```xml
SiTarget_NLayers_max  = 120                                   <!-- techo duro -->
SiTarget_dim_z        = SiTarget_NLayers_max * SiTarget_spacing   <!-- 1320 mm -->
SiTarget_SiPad_gap    = 2*env_safety                           <!-- input nuevo, 0.002 mm -->
SiTarget_z_position   = -SiTarget_dim_z/2 - SiPad_dim_z - SiTarget_SiPad_gap
SiPad_NLayers_max     = 22
```

Poner 120 no bastaba: con el envelope de 1830 mm las capas quedaban **centradas** en
−1940…−630, desplazadas 250 mm de donde el MC las tiene. El envelope ahora **es** sus capas.

La separación entre detectores no existía como input (la fórmula cableaba cero); ahora es
`SiTarget_SiPad_gap`. Su valor por defecto es 0.002 mm y no 0 porque cada envelope se infla
por `env_safety` y a separación literalmente cero ROOT reporta un solapamiento.

`min()` no está disponible ni en `parse_geometry.py` (eval sin builtins) ni garantizado en
DD4hep, así que el techo es un literal: si `SiTarget_spacing` crece por encima de 15.25 mm
las capas dejan de caber y hay que revisarlo a mano.

### 2.3 Geometría final verificada

| detector | capas | paso | bordes z |
|---|---|---|---|
| SiTarget | 120 | 11.00 mm | **−1690.00 … −370.00** |
| SiPad | 22 | 16.55 mm | −370 … 0 (capas en −367.05 … −2.95) |
| MTC | 42 (3 × 14) | 80.00 mm | +1 … +3361 |

**0 solapamientos**. `SiDetector_total_dim_z = 2200` pasa a ser un presupuesto: 1690 mm
ocupados, **510 mm libres aguas arriba** para grosores, gaps o separación.

Pendiente: el GDML sigue siendo un export del 28 de julio. Coincide en capas pero no se ha
regenerado tras el cambio de posición.

---

## 3. Tracking — el problema del seeder y su rediseño

### 3.1 Lo que había

El `findSeeds` original acumulaba un histograma 2D de las **posiciones transversales**
(x, y) de los puntos, descartando z. Un pico significaba "muchos hits en la misma (x,y) a
distintas z": **solo podía encontrar trazas paralelas al haz**. El seed no llevaba
dirección; el CKF arrancaba con `seedDir = (1, 0, 0)` cableado.

Distribución angular de las partículas cargadas de estado final en GENIE:

| | mediana | <30 mrad | <100 mrad | <200 mrad |
|---|---|---|---|---|
| todas (n=11534) | 921 mrad | 1.1 % | 6.7 % | 14.0 % |
| la más energética por evento | 243 mrad | 4.6 % | 25.0 % | 43.4 % |

La aceptancia del acumulador era de ~30–60 mrad. En el evento 1 (muón de 3.09 GeV a 183
mrad) `findSeeds` devolvía cero seeds, y el `if (seeds.empty()) continue;` de la línea 1625
no lo registraba: el modo de fallo dominante era invisible.

### 3.2 Lo que limita el diseño

El MTC lleva 50 mm de hierro por capa de 80 mm. Residuo RMS de un ajuste de recta a las
posiciones de verdad de 199 muones de >1 GeV:

| capas | brazo | RMS y (med/p90) | RMS x (med/p90) |
|---|---|---|---|
| 3 | 160 mm | 0.51 / 5.74 mm | 1.81 / 6.06 mm |
| 4 | 240 mm | 1.11 / 6.79 mm | 3.27 / 8.01 mm |
| 6 | 400 mm | 2.83 / 9.52 mm | 4.51 / 11.00 mm |
| 10 | 720 mm | 4.28 / 11.19 mm | 6.41 / 14.87 mm |
| 20 | 1520 mm | 6.82 / 12.39 mm | 9.95 / 16.16 mm |

`y` es la proyección sin campo (B ∥ Y): su residuo es dispersión múltiple pura. Campo y
scattering son **comparables**, así que no hay brazo donde la curvatura sea medible y la
asociación fiable a la vez. El seeder ajusta rectas sobre ventanas cortas; el momento es
trabajo del CKF, que propaga con el campo real losa a losa.

### 3.3 El seeder nuevo

- `SeedCandidate` lleva `tx, ty`; el CKF arranca con `seedDir(1, ty, tx)` en la superficie
  donde empieza el tracklet.
- Detectores 2D (SiPad, MTC `AllPairs`): Hough lineal en y–z y luego en x–z solo sobre los
  puntos que apoyan el pico de y. Asociación exacta, sin fantasmas.
- Detectores de tiras (SiTarget): Hough por proyección; desaparecen los cruces X×Y y el
  tope de 500 cruces/estación que descartaba estaciones enteras en silencio.
- Ventanas deslizantes de 4 capas con paso 2 (6 capas era peor: 33 vs 39 eventos con traza
  en 60, porque un muón cambia `dtx ≈ 0.24` en 10 capas y una ventana larga devuelve la
  cuerda, no la tangente).
- Parámetros por detector (`*PerDet`), porque el paso de capa va de 11 mm a 80 mm.
- Interceptos fuera de rango descartados en vez de saturados al bin del borde;
  `HoughHalfSize` 200 → 350 mm (la envolvente del MTC es 300).
- Filtro de duplicados por número de measurements y luego chi2/ndf.
- El `seeds.empty()` mudo registra una línea de debug.

### 3.4 El bug de signo del marco local (encontrado por la tarde)

Todas las superficies planas que construye `ACTSGeoSvc` llevan `rot90Y = AngleAxis(π/2, Y)`.
Sus ejes locales, en el global de ACTS (haz, Y, X de DD4hep), son e0 = (0, 0, −1) y
e1 = (0, 1, 0): **loc0 = −x, loc1 = +y**. Las semillas pasan por
`BoundTrackParameters::create(global → local)` y respetan ese marco; el calibrador recibía
`x` en bruto de los convertidores y lo pasaba como loc0 tal cual. Semilla y medidas vivían en
marcos espejo, a 2|x| de distancia en loc0.

Con el seeder viejo (trazas paralelas junto al eje) 2|x| son unos milímetros y cabían en la
σ de 10 mm de la semilla. Con el nuevo, que encuentra trazas inclinadas en ±300 mm, es
fatal: en el evento 1 la predicción del CKF caía en loc0 = +136.9 mientras el hit del muón
estaba en x = −138.0; el selector tomaba un hit de un π⁺ a 29 mm y ningún seed sobre el muón
producía nunca un segmento. La verificación fue numérica (`R[:,0]` de la rotación) y la
prueba, la traza VERBOSE del CKF.

Corregido en tres sitios: (1) `SeqCalibrator`: 2D → (−x, y), tira X del SiTarget → −x,
1D del MTC → −s; (2) `MTCSciFiMeasConverter`: la inversión estéreo dividía por 2cosα/2sinα
sobre el centro de tira a y = 0 en vez de por 2/2tanα — un sesgo de escala del 0.4 % en x
e y, 1 mm en el borde, mayor que el paso de fibra; (3) la escritura de `TrackStates` (PHASE
4) convertía local → DD4hep con trigonometría a mano que asumía x = +loc0, así que cada
estado se dibujaba reflejado — es lo que se vio en el display como "el muón va hacia el
otro lado en x". Ahora pasa por `surface.localToGlobal()`. El display tenía la misma
fórmula para re-emparejar estados U/V; también corregido.

También se corrigió el mapa de campo de las steerings: usaba ±300 mm en x para las tres
estaciones cuando el núcleo magnetizado de hierro mide 400/500/600 mm (el acero que lo
flanquea no lleva campo en Geant4). Solo, este arreglo no recuperaba el evento 1.

---

## 4. Tracking — resultados

### 4.1 Evolución a lo largo del día (120 capas, mismos 2000 eventos)

| estado | ev ≥1 global | globales/ev | =1 global | >5 globales | eff muón |
|---|---|---|---|---|---|
| mañana: seeder nuevo, bug de signo | 51.7 % | 1.27 | 24.7 % | — | — |
| + signo corregido (MTC min 3, χ² 5/10) | 65.2 % | 2.20 | 28.7 % | 14.4 % | 66.0 % |
| + aceptancia MTC | 63.9 % | 2.30 | 28.6 % | 16.2 % | 65.8 % |
| + MTC min 5, χ² 8/16 | 56.4 % | 1.37 | 32.5 % | 6.2 % | 65.0 % |
| + scint 30 mm + purga MTC 30 mm/8 | 49.7 % | 0.72 | 36.4 % | 0.6 % | 65.4 % |
| **+ borde de la purga (final)** | **48.6 %** | **0.67** | **36.5 %** | **0.2 %** | **66.2 %** |

La lectura correcta de la columna "ev ≥1 global" es que **baja a propósito**: cada paso
desde el 65.2 % quita fragmentos de shower de 3–5 hits, no muones — la eficiencia de muón
es plana dentro del ±1 % en todas las filas. La columna que mejora es la de "exactamente
una global", que es lo que un evento CC típico debería dar.

### 4.2 Configuración final

| | valor |
|---|---|
| eventos que dejan hits | 1965 / 2000 |
| **eventos con ≥1 global** | **48.6 %** (49.5 % de los que dejan hits) |
| globales por evento | **0.67** |
| globales spliced / promoted | 197 / 1135 |
| locales SiT / SiP / MTC (no promovidas) | 239 / 284 / 131 |
| fusiones intra-detector | 389 |
| **eficiencia de muón** (1117 eventos con muón que deja hits) | **66.2 %** |
| eficiencia de muón, ≥20 hits del muón (733 eventos) | **85.5 %** |

Configuración: `MinMeasPerSegment = [5,3,5]`, `SegMaxChi2TightPerDet = [10,3,8]`,
`SegMaxChi2LoosePerDet = [20,6,16]`, `HitPurgeWindow = [1,8,30]`,
`HitPurgeMaxNeighbors = [8,4,8]`, `HitPurgeMaskBorder = True`,
`MTCSciFiMeasConverter.ScintConfirmWindow = 30`, `AcceptanceTolerance = 2`,
`DuplicateSharedCellIDs = 1`, `MergeSegments = True`, `MergeMaxOverlap = 3500`,
`ScaleMinMeasByDepth = True`, `PromoteUnspliced = True`.

### 4.3 Distribución de globales por evento

| nº globales | % eventos |
|---|---|
| 0 | 49.6 % (de ellos 1.75 % sin hits) |
| **1** | **36.5 %** |
| 2 | 8.0 % |
| 3 | 3.0 % |
| 4 | 0.7 % |
| 5–6 | 0.2 % |
| >6 | 0.1 % |

### 4.4 Eficiencia de muón por número de hits que deja el muón

| hits del muón (todos los detectores) | eventos | eff |
|---|---|---|
| <6 | 85 | 1 % |
| 6–10 | 109 | 10 % |
| 10–20 | 190 | 53 % |
| 20–50 | 315 | 85 % |
| ≥50 | 418 | 86 % |

La mitad de los muones perdidos dejan ≤10 hits en total (blandos, que salen de lado o se
paran): eso no es ineficiencia del algoritmo. El plateau del 85 % para muones largos es la
ineficiencia real que queda. Definición: un evento cuenta si alguna local aceptada tiene al
muón como partícula MC dominante (≥50 % de sus medidas). Sin corte de pureza adicional.

---

## 5. Tracking — pureza y fragmentación (verdad MC)

Método: `DumpSegmentCellIDs` vuelca los cellIDs de cada local aceptada; se cruzan contra
el parámetro de frame `*ContribParticleIDs`. **Limitación**: para `AllPairs` el converter
propaga el cellID de la tira U, así que un par fantasma U/V se ve puro. **La pureza del MTC
es una cota superior**; las de SiTarget (tiras sin emparejar) y SiPad (hits 2D) son limpias.

### 5.1 Configuración final, 2000 eventos

| detector | eventos con local | locales | **locales/partícula MC** | **pureza mediana** | ≥2 partículas | nMeas mediana | completeness |
|---|---|---|---|---|---|---|---|
| SiTarget | 245 | 461 | 1.49 | 0.88 | 70 % | 29 | 0.87 |
| SiPad | 438 | 957 | 1.44 | 1.00 | 36 % | 7 | 1.00 |
| **MTC** | 777 | 1551 | **1.84** | **1.00** | **17 %** | **8** | 0.94 |

Comparado con la mañana (bug de signo, MTC min 3): MTC 2.69 → 1.84 locales por partícula,
pureza 0.75 → 1.00, mezcladas 56 % → 17 %, nMeas mediana 4 → 8. SiPad y SiTarget apenas se
mueven (1.50 → 1.44 y 1.56 → 1.49).

### 5.2 Pureza frente a longitud y χ² (tras el signo, antes del retuning del MTC)

Sobre 7415 locales del MTC con mínimo 3:

| nMeas | n | pureza media | % puras |
|---|---|---|---|
| 3 | 2937 | 0.77 | 46 % |
| 4 | 2068 | 0.76 | 40 % |
| 5–6 | 1386 | 0.72 | 30 % |
| 7–8 | 326 | 0.80 | 53 % |
| 9–11 | 237 | 0.94 | 79 % |
| ≥12 | 461 | 0.97 | 79 % |

| χ²/ndf | n | % puras |
|---|---|---|
| [0, 1) | 1391 | 49 % |
| [1, 3) | 2559 | 50 % |
| [3, 5) | 1877 | 43 % |
| [5, 10) | 1588 | 34 % |

La longitud discrimina; el χ² apenas. El 67 % de las locales del MTC eran de 3–4 hits y
menos de la mitad de ellas puras. De ahí el mínimo de 5 — y las puertas laxas, porque con
el marco correcto un muón puro de 11 hits ajusta a χ²/ndf ≈ 6 en 50 mm de Fe por capa, y
la puerta de 5 lo cortaba (el CKF lo encontraba entero y lo rechazaba).

## 6. Tracking — los barridos de cortes

Los barridos de la mañana (§6.1–6.3 de la primera edición) se hicieron con el bug de signo
y quedan obsoletos: sus conclusiones ("subir el mínimo del MTC está dominado", "el veto de
scint pierde eventos") eran artefactos de un tracker que solo veía trazas junto al eje. Se
rehicieron todos por la tarde, sobre 2000 eventos, en HTCondor, y con la eficiencia de
muón como juez.

### 6.1 MTC: mínimo de medidas × puertas χ² (11 puntos)

| MTC min / puertas | ev ≥1 global | globales/ev | =1 | >5 | **eff muón** |
|---|---|---|---|---|---|
| 3 / 5,10 | 63.9 % | 2.30 | 28.6 % | 16.2 % | 65.8 % |
| 3 / 8,16 | 65.2 % | 2.72 | 28.7 % | 19.9 % | — |
| 4 / 5,10 | 59.8 % | 1.67 | 33.3 % | 9.1 % | 63.8 % |
| 5 / 5,10 | 53.8 % | 1.07 | 33.5 % | 3.1 % | 61.2 % |
| **5 / 8,16** | **56.4 %** | **1.37** | **32.5 %** | **6.2 %** | **65.0 %** |
| 6 / 8,16 | 52.0 % | 1.00 | 32.6 % | 2.8 % | 63.4 % |
| 8 / 8,16 | 45.7 % | 0.76 | 31.7 % | 1.4 % | 60.2 % |
| 5 / 5,10, SiPad min 5 | 52.9 % | 1.02 | 34.2 % | 2.8 % | — |

3 → 5 quita 8 puntos de "eventos con global" y 0.8 de muón. Subir SiPad no aporta.

### 6.2 Limpieza del pool del MTC (6 puntos + borde)

Diagnóstico previo, sobre verdad MC: el pool del MTC en eventos con shower es ~96 %
fantasmas U×V (1817 cruces/evento; N fibras por plano dan N² cruces), y los fantasmas no
son densos — el cruce de dos fibras separadas Δx cae a 5.7·Δx en y, repartido por todo el
plano. Una purga de densidad por capa apenas discrimina **por hit** (quita ~9 % de los hits
del propio muón, que sale del vértice dentro de la shower, contra ~17 % del resto), y la
energía del scintillator tampoco (las celdas sin muón tienen mediana 0.12 MeV: shower
blanda en muchas celdas, no pocas celdas calientes).

| limpieza MTC | medidas/ev | ev ≥1 global | globales/ev | =1 | >5 | **eff muón** |
|---|---|---|---|---|---|---|
| ninguna | 1086 | 56.4 % | 1.37 | 32.5 % | 6.2 % | 65.0 % |
| scint 30 mm | 836 | 55.6 % | 1.25 | 33.5 % | 5.3 % | 64.9 % |
| scint 40 mm | 915 | 55.7 % | 1.29 | 33.0 % | 5.3 % | 65.2 % |
| purga 20 mm / 4 | 368 | 50.8 % | 0.76 | 36.0 % | 0.7 % | 64.0 % |
| purga 30 mm / 8 | 364 | 50.6 % | 0.76 | 36.1 % | 0.6 % | 64.5 % |
| scint 30 + purga 30/8 | 309 | 49.7 % | 0.72 | 36.4 % | 0.6 % | 65.4 % |
| **+ borde (DBSCAN)** | **272** | **48.6 %** | **0.67** | **36.5 %** | **0.2 %** | **66.2 %** |

Por hit la purga no discrimina; **por traza sí**: el CKF atraviesa los huecos del muón,
mientras que los fragmentos de 3–5 hits pierden el cúmulo del que se construían. El borde
(quitar también los supervivientes a menos de la ventana de un hit purgado, propuesta del
usuario) es una mejora pequeña en todo y ninguna pérdida.

Detalle de infraestructura: `k4run --Alg.PropBool False` pone **True** (cadena no vacía);
la comparación con/sin borde salió idéntica por eso y se verificó localmente antes de
concluir nada. Los booleanos se cambian en el steering.

### 6.3 Fantasmas U/V: qué se resuelve y qué no

Tres filtros en cascada, todos activos en producción:

1. **Aceptancia del plano** (`AcceptanceTolerance = 2 mm`): el cruce debe caer en la
   superficie ACTS de su capa (±300 mm). Mata los geométricamente imposibles (|y| ≈ 1.7 m
   en el evento 1: 22 de 102 medidas). Coste en eficiencia de muón: 0.2 puntos.
2. **Confirmación por celda de scintillator** (30 mm): el cruce debe caer en una celda
   encendida de su propia capa. Celdas de 50 mm frente a fibras de 1 mm: resuelve los
   fantasmas lejos de cualquier partícula, no los de dentro de la mancha de la shower.
3. **CKF**: los que quedan dentro del plano los rechaza el χ² traza a traza.

Lo que ningún filtro a posteriori resuelve: fantasmas dentro de una celda con varias
partículas. La vía correcta es usar el scint **en el emparejamiento** (para cada celda
encendida, emparejar solo las fibras U y V que la cruzan; la energía de la celda, ~3.5
MeV por MIP, dice cuántos cruces reales esperar). Pendiente.

### 6.4 Fusión intra-detector

Mi primer diseño exigía que un segmento terminara antes de que empezara el otro. **Los
fragmentos se solapan en z** (el purgado por seed hace que el siguiente seed reconstruya la
partícula con los hits que sobraron, intercalados en el mismo rango): 15 012 pares
rechazados por el gap frente a 198 que llegaban al test. Con solapamiento permitido:

| | eficiencia | trazas/ev | fusiones |
|---|---|---|---|
| solo contigüidad | 53.5 % | 1.54 | 1 |
| solapamiento 3500 mm | **53.5 %** | **1.42** | **79** |

Misma eficiencia hasta el evento. Es el único cambio que reduce trazas sin gastar
eficiencia, porque une en vez de descartar. Sobre los 2000 eventos de la geometría correcta, configuración
final: 389 fusiones.

### 6.5 Escalado del mínimo por profundidad

`ScaleMinMeasByDepth` capa el suelo de medidas por las capas disponibles aguas abajo del
punto de entrada. Implementado en el bypass del detector y en la aceptación del segmento.
**Efecto medido: −1 evento, −1 traza.** Es un no-op en GENIE porque los suelos del MTC y
SiPad ya están en el mínimo del CKF. Se deja activado por principio; no hay evidencia de que
aporte en esta muestra.

---

## 7. MOBO — infraestructura y espacio de búsqueda

### 7.1 Cadena GENIE

- `run_sim_genie.py`: ddsim desde el HepMC de GENIE con geometría variable. **El HepMC es
  compartido**: todos los trials replican los mismos primarios.
- `job4_tracking_seq.py` / `job5_rntuple_seq.py`: ports con `$SND_COMPACT` y
  `$SND_OUTPUT_LEVEL`, regenerados desde el steering de la pipeline.
- Métrica `eff_track` = eventos con ≥1 `GlobalTrack` / eventos pedidos, leyendo
  `window_id` del RNTuple. El denominador es el pedido, no el observado; `n_events_seen`
  (de los contenedores de hits) detecta truncamientos.
- El payload borra los ROOT intermedios de cada trial tras escribir `metrics.json`: de
  37 MB a ~1 MB por trial.

### 7.2 Inversión del template

Antes: `dim_z` era input y el gap entre capas se derivaba (`auto`) para rellenar una
profundidad fija. Ahora:

```
dim_z          = NLayers × LayerThickness          (derivado)
LayerThickness = W + offsets + sensores + holgura   (todo inputs)
SiDetector_total_dim_z = 2200                       (TECHO, comprobado en el dry-run)
```

| parámetro | rango | baseline |
|---|---|---|
| `SiTarget_NLayers` | 20 – **120** | 120 |
| `SiPad_NLayers` | 4 – **22** | 22 |
| `SiTarget_WThickness` | 2 – 8 mm | 3.5 |
| `SiPad_WThickness` | 5 – 30 mm | 10 |
| `SiTarget_layer_slack` | 1 – 15 mm | 4.9 |
| `SiPad_layer_gap` | 0.5 – 10 mm | 1.0 |
| `SiTarget_XY_plane_gap_frac` | 0 – 1 | 1 |
| `SiTarget_SiPad_gap` | 0.002 – 200 mm | 0.002 |

Fijos: `SiPad_frame_gap = 0.1`, `SiTarget_module_offset = 1`, `SiPad_env_slack = 5.9`
(los 5.9 mm de envelope sin usar del SiPad, necesarios para reproducir el compact).

El baseline renderizado reproduce **exactamente** el compact: SiTarget −1690…−370, SiPad
−370…0, presupuesto 1690 de 2200.

### 7.3 Objetivos

| objetivo | dirección | referencia |
|---|---|---|
| `eff_track` | max | 0.5× baseline |
| `nhits_sipad` | max | 0.5× baseline |
| `cost_proxy` | min | 2.0× baseline |

El coste es un proxy (30 k€/m² Si, 0.1 k€/kg W). Sin él, los otros dos no están en
tensión y el frente colapsaría a un punto.

---

## 8. MOBO — resultados

**Nota**: la campaña corrió con el tracking de la mañana (seeder nuevo, bug de signo, MTC
min 3). `eff_track` es "eventos con ≥1 global" bajo ese tracker, así que su escala (0.5–0.6)
no es comparable con el 48.6 % de §4.2, y sobreestima a las geometrías que dan fragmentos.
El ranking relativo entre geometrías es el que es; la próxima campaña usa el tracking final
(el port `Optimization/Analysis/job4_tracking_seq.py` ya está alineado).

Campaña `snd_genie_100`: 100 trials, **96 completados**, 3 infeasible (rechazados por el
gate antes de encolar), 1 fallido (job atascado en Condor sin llegar a arrancar).
Diseño inicial: 1 baseline + 18 Sobol (`2×(d+1)` con d=8). Fase bayesiana: 81 trials.

### 8.1 Resumen por fase

| fase | completados | eff_track med (max) | nhits_sipad med (max) | cost med (min) | en el frente |
|---|---|---|---|---|---|
| baseline + Sobol | 16 | 0.545 (0.610) | 1970 (3549) | 1016 (341) | 1 |
| **qLogNEHVI** | 80 | 0.550 (0.595) | **4202** (6164) | **485** (271) | **15** |

La fase bayesiana dobla los hits en SiPad y divide el coste por dos respecto a la Sobol,
con la misma eficiencia mediana. 15 de los 16 puntos del frente son suyos.

### 8.2 Frente de Pareto final

Columnas: capas SiTarget / SiPad, grosor W SiTarget / SiPad (mm), holgura SiTarget (mm),
gap SiPad (mm), fracción XY, separación SiTarget–SiPad (mm).

| trial | fase | eff_track | nhits_sipad | cost_proxy | ST cap | SP cap | ST W | SP W | ST holg | SP gap | frac | ST-SP |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| t00006 | sobol | **0.610** | 2360 | 1381 | **115** | 12 | 3.96 | 6.43 | 2.02 | 9.19 | 0.80 | 49.6 |
| t00045 | qlognehvi | **0.595** | 2780 | 1123 | **92** | 13 | 2.32 | 9.58 | 2.45 | 5.04 | 0.59 | 183.6 |
| t00088 | qlognehvi | 0.585 | 4112 | 418 | 20 | 17 | 2.00 | 7.11 | 4.40 | 2.13 | 0.00 | 200.0 |
| t00022 | qlognehvi | 0.575 | **5682** | 456 | 20 | 22 | 2.00 | 5.00 | 15.00 | 0.50 | 0.00 | 191.6 |
| t00023 | qlognehvi | 0.560 | 1834 | **271** | 20 | 5 | 3.72 | 5.00 | 1.00 | 10.00 | 0.00 | 200.0 |
| t00092 | qlognehvi | 0.555 | 4530 | 447 | 20 | 19 | 2.00 | 7.43 | 3.23 | 2.46 | 0.40 | 200.0 |
| t00020 | qlognehvi | 0.550 | 4568 | 365 | 20 | 14 | 2.24 | 5.00 | 2.90 | 7.62 | 0.00 | 200.0 |
| t00027 | qlognehvi | 0.550 | 3659 | 317 | 20 | 10 | 2.00 | 5.00 | 3.78 | 0.50 | 0.00 | 200.0 |
| t00091 | qlognehvi | 0.550 | 4724 | 410 | 20 | 18 | 2.00 | 5.00 | 3.47 | 0.50 | 0.00 | 186.9 |
| t00031 | qlognehvi | 0.545 | 4059 | 363 | 20 | 14 | 2.00 | 5.00 | 1.24 | 2.20 | 0.93 | 60.6 |
| t00026 | qlognehvi | 0.540 | **6164** | 456 | 20 | 22 | 2.00 | 5.00 | 15.00 | 0.50 | 0.00 | 65.8 |
| t00029 | qlognehvi | 0.540 | 5184 | 433 | 20 | 20 | 2.00 | 5.00 | 15.00 | 0.50 | 1.00 | 0.002 |
| t00090 | qlognehvi | 0.540 | 4873 | 410 | 20 | 18 | 2.00 | 5.00 | 5.78 | 0.50 | 0.00 | 200.0 |
| t00068 | qlognehvi | 0.535 | 4923 | 421 | 20 | 19 | 2.00 | 5.00 | 11.90 | 10.00 | 0.00 | 30.7 |
| t00024 | qlognehvi | 0.510 | 4726 | 375 | 20 | 15 | 2.00 | 5.00 | 15.00 | 10.00 | 0.00 | 200.0 |
| t00036 | qlognehvi | 0.490 | 2295 | 271 | 20 | 6 | 2.00 | 5.00 | 15.00 | 0.50 | 0.57 | 200.0 |

**Baseline (t00000): eff_track 0.535, nhits_sipad 3109, cost_proxy 1585. Está dominado** —
no aparece en el frente.

### 8.3 Lectura del frente

**14 de los 16 puntos tienen `SiTarget_NLayers = 20`, el mínimo del rango.** Y en esos 14,
el grosor del tungsteno está casi siempre en el mínimo (2.0 y 5.0 mm) y la separación
entre detectores en el máximo (200 mm). El optimizador convergió a "el SiTarget más pequeño
y ligero posible, lo más lejos posible del SiPad".

Eso es exactamente la firma del sesgo de la sección 9: un SiTarget de 20 capas se apunta
el 19.4 % de vértices que caen en aire, y esos eventos le llegan sin haber atravesado
material. **No es un resultado físico sobre el tamaño del detector; es el optimizador
encontrando el regalo.**

Los dos puntos que **sí** informan sobre física son los únicos con muchas capas:

- **t00006** (Sobol, 115 capas): la mayor eficiencia del frente, 0.610, a coste 1381.
- **t00045** (qLogNEHVI, 92 capas): 0.595 a coste 1123 — un 28 % más barato que el baseline
  con 6 puntos más de eficiencia. Su receta: SiTarget algo más corto (92 vs 120), tungsteno
  fino (2.3 mm), SiPad de 13 capas gruesas (9.6 mm), y 184 mm de separación.

Lo que sí se puede afirmar con este frente: **el baseline no es óptimo bajo este modelo de
coste** — t00045 lo domina con 92 capas, donde el sesgo vale ~5 %, no 19 %.

Sobre `nhits_sipad`: los puntos a 6164 y 5682 hits (t00026, t00022) tienen 22 capas de SiPad
con 15 mm de holgura en el SiTarget — el SiTarget más largo posible sin material, para que
la ducha llegue intacta al SiPad. Es coherente, pero es la misma trampa vista desde el otro
objetivo.

---

## 9. El sesgo del muestreo de vértices

GENIE generó los vértices **una vez**, muestreándolos según el material del detector de 120
capas. Ese HepMC lo comparten todos los trials, lo que hace justa la comparación en la
propagación. Pero ddsim coloca cada primario en su posición **sea cual sea el material que
haya ahí**: no comprueba nada.

Un SiTarget más corto se pega igualmente al SiPad; su cara aguas arriba se desplaza y los
vértices que quedan detrás caen en aire. Geant4 propaga igual sus estados finales. La
variante **se acredita interacciones que su material no habría producido**, y encima le
llegan limpias, sin ducharse en material aguas arriba.

| capas SiTarget | profundidad | cara aguas arriba | vértices en aire |
|---|---|---|---|
| 120 | 1320 mm | −1690 mm | 0 (0.0 %) |
| 100 | 1100 mm | −1470 mm | 79 (4.0 %) |
| 80 | 880 mm | −1250 mm | 157 (7.8 %) |
| 60 | 660 mm | −1030 mm | 240 (12.0 %) |
| 40 | 440 mm | −810 mm | 321 (16.1 %) |
| **20** | 220 mm | −590 mm | **388 (19.4 %)** |

El sesgo es **simétrico**: añadir capas aguas arriba pondría material donde no hay
vértices y esas capas parecerían inútiles. Y no lo corrige ninguna reducción de trazas ni
de coste: afecta a la **tasa**, que es proporcional a la masa atravesada.

**¿Hay que resimular GENIE por geometría?** Estrictamente sí, por dos motivos: la tasa
total y la distribución espacial de los vértices dependen del material. Pero depende de la
pregunta:

- "Dado un detector de esta profundidad, ¿cómo reparto el muestreo?" → fijar la extensión y
  variar el interior. Sesgo despreciable, sin resimular.
- "¿Cuánto detector necesito?" → la tasa es parte de la respuesta y hay que regenerar. Eso
  necesita la instalación de GENIE, los ficheros de flujo y un GDML por variante.

Un arreglo intermedio, no aplicado en esta campaña: descartar de numerador y denominador
los eventos cuyo vértice caiga fuera del material de la variante. Corrige la extensión a
primer orden (la densidad de vértices es uniforme, ~35/100 mm); no corrige cambios de
grosor de W.

**Decisión tomada para esta campaña**: dejarla correr y explicar el sesgo. El frente es una
cota superior para los puntos de pocas capas.

---

## 10. Tiempos

### 10.1 Reconstrucción de 2000 eventos GENIE (un proceso, lxplus)

| etapa | reloj |
|---|---|
| `job1_overlay` | 18 s |
| `job4_tracking_seq` | **~20 min** |
| `job5_rntuple_seq` | 52 s |
| **total** | **~21 min** |

Cadena completa con ddsim (Condor, un job): ~50 min.

### 10.2 Campaña MOBO de 100 trials

| | valor |
|---|---|
| **reloj de campaña** | **3 h 45 min** (00:43 → 04:29) |
| wall por trial | min 189, **mediana 295**, max 540 s |
| CPU por trial | mediana 275 s |
| cómputo agregado | 8.1 h wall, 7.5 h CPU |
| en vuelo | 4 |
| driver (lxplus, 1 hilo, nice 19) | 3–6 % de CPU, 420 MB |

De los 225 min, unos 20 los pasó el loop esperando al único trial atascado en Condor. Sin
él, ~3 h 25 min. Los trials son más rápidos que en la campaña anterior (295 vs 390 s de
mediana): el SiTarget de 120 capas simula más rápido que el de 166.

---

## 11. Lo que salió mal en la infraestructura

Tres avisos de IT en dos días, todos por tratar lxplus como máquina de cálculo.

**1. Un barrido en paralelo en lxplus** (`xargs -P 5`), 49.95 % de presión en lxplus990.

**2. La cuota de AFS llena** (10 GB al 100 %). Cada trial guardaba 37 MB de ROOT; 200 trials
son 7.4 GB. ddsim simulaba bien y moría al escribir (`SysError in TFile::Flush`), y el loop
lo reportaba como `job1_failed` tres capas por encima de la causa. Arreglo: el payload
borra los intermedios tras `metrics.json`.

**3. Dos drivers de MOBO a la vez a 99 % de CPU**, 44 CPU-horas en 4 horas en lxplus924.
Dos fallos encadenados: (a) el driver no es un sondeador ligero, es BoTorch ejecutando
`optimize_acqf` sobre 8 dimensiones con el kernel fusionado caído a Python; (b) mi
`pkill -f "mobo-resume"` coincidía con mi propio shell, lo mataba a él y el driver
sobrevivía. Arreglo: un solo driver, `OMP_NUM_THREADS=1`, `nice -n 19`, `ulimit -t`, kill
por PID leído a fichero, y verificación en el nodo donde vive el proceso.

**4. `max_trials` no acompañó a la dimensión del espacio.** El diseño inicial es
`2×(d+1)`; al pasar de 6 a 8 dimensiones consumió 19 de los 20 trials de la primera
campaña, dejando **un solo** trial bayesiano.

**5. Un job de Condor atascado** en estado Running con `NumJobStarts=0` durante 2 h 21 min.
El executor no lo declaraba perdido porque `condor_q` seguía devolviendo estado 3 (REMOVED)
tras eliminarlo, y cualquier estado no nulo contaba como "presente". Arreglado en
`htcondor.py`: REMOVED y COMPLETED cuentan como ausentes.

Todo está escrito en `CLAUDE.md` para que no dependa de acordarse.

---

## 12. Estado del código y pendientes

### Modificado (sin commitear)

```
gaudi_source/ACTSSequentialTracker.cpp          seeder nuevo, dedup cellID, fusión, escalado
simulation/geometry/SND_compact.xml             techos de capas, envelope = capas, gap input
Optimization/Simulation/geometry/SND_compact_template.xml   parameterización invertida
Optimization/Simulation/geometry/config.py      presupuesto como techo, auto vacío
Optimization/Simulation/geometry/parameters_template.yaml   baseline nuevo
Optimization/mobo/src/mobo/ship/{metrics,evaluator,payload}.py   eff_track, cadena GENIE, limpieza
Optimization/mobo/src/mobo/exec/htcondor.py     REMOVED/COMPLETED = ausente
Optimization/mobo/conf/experiment/snd_genie.yaml   espacio de 8D, techos, 100 trials
gaudi_jobs/genie_pipeline/job4_tracking_seq.py  configuración final
gaudi_source/MTCSciFiMeasConverter.cpp          inversión estéreo, aceptancia, veto scint activo
gaudi_source/ACTSSequentialTracker.cpp          + signo loc0 en el calibrador, estados por localToGlobal, borde de la purga
event_display/event_display_eve.py              re-emparejado U/V por (x,y) escritos
```

Nuevos: `run_sim_genie.py`, `Optimization/Analysis/job{4,5}_*_seq.py`,
`conf/evaluator/snd_genie.yaml`, `CLAUDE.md`, `simulation/run_script/steer/*.{sh,sub}`.

### Pendientes, por orden de impacto

1. **Emparejamiento U/V guiado por scintillator** en `MTCSciFiMeasConverter` (§6.3): es la
   única vía para los fantasmas dentro de la mancha de la shower, y el 15 % de muones
   largos que aún se pierden vive ahí.
2. **Relanzar MOBO con el tracking final** (con corte de vértice en material si la
   pregunta sigue siendo el tamaño del detector), en Condor también el driver.
3. **`ACTSProtoTracker` (CKF global de las otras pipelines) tiene el mismo bug de signo**
   en su calibrador y en su escritura de estados; no se ha tocado porque no está en la
   cadena GENIE. Corregir antes de volver a usarlo.
4. **22 tests de `Optimization/mobo/tests/ship/` fallan** (parameterización vieja) y
   `geometry.py:limits()` está muerto.
5. **Regenerar el GDML** desde el compact corregido.
6. **Pureza del MTC como cota superior**: una medida honesta necesita la posición 2D del
   cruce en el `SegDump`.
7. **Las otras 5 pipelines siguen sin validar** con el seeder nuevo ni con el signo.

---

## Anexo — todos los trials

### Fase inicial: baseline + Sobol (t00000–t00018)

| trial | eff_track | nhits_sipad | cost_proxy | ST cap | SP cap | ST W | SP W | ST holg | SP gap | frac | ST-SP | Pareto |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| t00000 | 0.535 | 3109 | 1585 | 120 | 22 | 3.50 | 10.00 | 4.90 | 1.00 | 1.00 | 0.002 |  |
| t00001 | 0.560 | 2461 | 989 | 68 | 15 | 4.97 | 12.85 | 14.39 | 1.29 | 0.10 | 69.49 |  |
| t00002 | 0.505 | 1141 | 1087 | 81 | 6 | 7.66 | 21.82 | 5.53 | 8.72 | 0.61 | 163.70 |  |
| t00003 | 0.580 | 3037 | 1252 | 97 | 18 | 2.09 | 10.94 | 3.18 | 6.76 | 0.43 | 122.84 |  |
| t00004 | 0.490 | 1457 | 578 | 34 | 9 | 5.28 | 23.85 | 8.11 | 4.08 | 0.91 | 15.84 |  |
| t00005 | 0.525 | 2186 | 801 | 26 | 20 | 7.17 | 28.68 | 10.88 | 1.80 | 0.29 | 142.19 |  |
| t00006 | 0.610 | 2360 | 1381 | 115 | 12 | 3.96 | 6.43 | 2.02 | 9.19 | 0.80 | 49.60 | ★ |
| t00007 | 0.495 | 1977 | 1349 | 89 | 17 | 5.79 | 20.45 | 6.67 | 6.40 | 0.24 | 90.48 |  |
| t00008 | 0.550 | 1707 | 658 | 50 | 8 | 3.08 | 14.58 | 11.59 | 3.76 | 0.72 | 195.86 |  |
| t00009 | 0.550 | 3549 | 822 | 51 | 19 | 6.37 | 6.70 | 4.70 | 4.90 | 0.36 | 133.78 |  |
| t00010 | 0.555 | 1601 | 1191 | 90 | 10 | 3.25 | 28.18 | 13.56 | 7.58 | 0.87 | 26.78 |  |
| t00012 | 0.530 | 1057 | 341 | 25 | 4 | 3.69 | 17.68 | 3.99 | 0.62 | 0.67 | 180.89 |  |
| t00013 | 0.520 | 1963 | 764 | 43 | 15 | 2.57 | 22.55 | 1.22 | 2.94 | 0.05 | 54.54 |  |
| t00016 | 0.585 | 1916 | 858 | 59 | 12 | 7.48 | 8.42 | 7.50 | 2.47 | 0.98 | 7.41 |  |
| t00017 | 0.540 | 2444 | 1071 | 62 | 21 | 3.40 | 19.82 | 8.80 | 5.48 | 0.58 | 108.63 |  |
| t00018 | 0.585 | 1672 | 1043 | 73 | 12 | 6.22 | 15.51 | 3.23 | 3.43 | 0.06 | 1.91 |  |

**Sobol** — 16 completados: eff_track med 0.545 (max 0.610), nhits_sipad med 1970 (max 3549), cost med 1016 (min 341); en el frente: 1

### Fase bayesiana: qLogNEHVI (t00019–t00099)

| trial | eff_track | nhits_sipad | cost_proxy | ST cap | SP cap | ST W | SP W | ST holg | SP gap | frac | ST-SP | Pareto |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| t00019 | 0.540 | 4442 | 458 | 20 | 22 | 2.38 | 5.00 | 1.40 | 7.47 | 0.62 | 200.00 |  |
| t00020 | 0.550 | 4568 | 365 | 20 | 14 | 2.24 | 5.00 | 2.90 | 7.62 | 0.000 | 200.00 | ★ |
| t00021 | 0.525 | 3663 | 491 | 20 | 22 | 7.81 | 5.00 | 11.80 | 0.50 | 0.000 | 167.90 |  |
| t00022 | 0.575 | 5682 | 456 | 20 | 22 | 2.00 | 5.00 | 15.00 | 0.50 | 0.000 | 191.60 | ★ |
| t00023 | 0.560 | 1834 | 271 | 20 | 5 | 3.72 | 5.00 | 1.00 | 10.00 | 0.000 | 200.00 | ★ |
| t00024 | 0.510 | 4726 | 375 | 20 | 15 | 2.00 | 5.00 | 15.00 | 10.00 | 0.000 | 200.00 | ★ |
| t00025 | 0.510 | 3376 | 1008 | 72 | 20 | 3.54 | 6.37 | 4.26 | 9.14 | 0.000 | 38.84 |  |
| t00026 | 0.540 | 6164 | 456 | 20 | 22 | 2.00 | 5.00 | 15.00 | 0.50 | 0.000 | 65.79 | ★ |
| t00027 | 0.550 | 3659 | 317 | 20 | 10 | 2.00 | 5.00 | 3.78 | 0.50 | 0.000 | 200.00 | ★ |
| t00028 | 0.545 | 5135 | 456 | 20 | 22 | 2.00 | 5.00 | 1.00 | 9.11 | 0.000 | 3.09 |  |
| t00029 | 0.540 | 5184 | 433 | 20 | 20 | 2.00 | 5.00 | 15.00 | 0.50 | 1.00 | 0.002 | ★ |
| t00030 | 0.530 | 3168 | 373 | 20 | 13 | 2.95 | 7.10 | 2.46 | 6.91 | 0.85 | 70.33 |  |
| t00031 | 0.545 | 4059 | 363 | 20 | 14 | 2.00 | 5.00 | 1.24 | 2.20 | 0.93 | 60.62 | ★ |
| t00032 | 0.530 | 4378 | 433 | 20 | 20 | 2.00 | 5.00 | 1.00 | 0.50 | 0.000 | 0.002 |  |
| t00033 | 0.555 | 5407 | 554 | 32 | 20 | 2.00 | 5.00 | 15.00 | 0.50 | 0.000 | 200.00 |  |
| t00034 | 0.560 | 3277 | 1391 | 111 | 19 | 2.90 | 6.69 | 2.23 | 8.41 | 0.93 | 84.97 |  |
| t00035 | 0.545 | 4448 | 386 | 20 | 16 | 2.00 | 5.00 | 15.00 | 2.16 | 0.000 | 200.00 |  |
| t00036 | 0.490 | 2295 | 271 | 20 | 6 | 2.00 | 5.00 | 15.00 | 0.50 | 0.57 | 200.00 | ★ |
| t00037 | 0.560 | 5587 | 456 | 20 | 22 | 2.00 | 5.00 | 15.00 | 0.50 | 1.00 | 200.00 |  |
| t00038 | 0.565 | 5439 | 456 | 20 | 22 | 2.00 | 5.00 | 15.00 | 10.00 | 0.000 | 200.00 |  |
| t00039 | 0.545 | 4425 | 490 | 20 | 22 | 2.00 | 7.78 | 15.00 | 0.50 | 0.000 | 200.00 |  |
| t00040 | *FAILED* | | | 20 | 22 | 2.00 | 5.00 | 15.00 | 10.00 | 0.000 | 0.002 | |
| t00041 | 0.575 | 2694 | 1341 | 113 | 13 | 3.26 | 5.44 | 2.60 | 9.86 | 0.43 | 37.93 |  |
| t00042 | 0.520 | 5266 | 456 | 20 | 22 | 2.00 | 5.00 | 15.00 | 0.50 | 0.000 | 200.00 |  |
| t00043 | 0.530 | 2569 | 811 | 50 | 18 | 6.36 | 8.00 | 10.71 | 0.73 | 0.99 | 0.002 |  |
| t00044 | 0.570 | 3786 | 980 | 82 | 13 | 2.00 | 5.00 | 1.20 | 3.78 | 0.44 | 200.00 |  |
| t00045 | 0.595 | 2780 | 1123 | 92 | 13 | 2.32 | 9.58 | 2.45 | 5.04 | 0.59 | 183.64 | ★ |
| t00046 | 0.535 | 4876 | 911 | 65 | 22 | 2.00 | 5.00 | 15.00 | 10.00 | 1.00 | 200.00 |  |
| t00047 | 0.585 | 3542 | 510 | 20 | 19 | 2.00 | 13.28 | 3.78 | 0.50 | 0.000 | 200.00 |  |
| t00048 | 0.490 | 2113 | 1245 | 99 | 12 | 4.35 | 9.95 | 2.33 | 2.60 | 0.02 | 200.00 |  |
| t00049 | 0.545 | 2739 | 830 | 61 | 14 | 2.00 | 11.49 | 1.26 | 0.50 | 0.000 | 200.00 |  |
| t00050 | 0.525 | 3914 | 1157 | 95 | 17 | 2.00 | 5.00 | 7.92 | 5.80 | 0.60 | 175.12 |  |
| t00051 | 0.535 | 5711 | 456 | 20 | 22 | 2.00 | 5.00 | 1.00 | 5.07 | 0.000 | 200.00 |  |
| t00052 | 0.560 | 3047 | 1217 | 102 | 14 | 2.75 | 5.00 | 2.91 | 7.64 | 0.69 | 107.64 |  |
| t00053 | 0.555 | 5397 | 456 | 20 | 22 | 2.00 | 5.00 | 4.44 | 0.50 | 0.000 | 195.54 |  |
| t00054 | 0.560 | 3980 | 957 | 73 | 19 | 2.00 | 5.00 | 3.39 | 5.03 | 0.42 | 181.82 |  |
| t00055 | 0.540 | 5405 | 456 | 20 | 22 | 2.00 | 5.00 | 13.00 | 10.00 | 1.00 | 159.65 |  |
| t00056 | 0.535 | 4807 | 456 | 20 | 22 | 2.00 | 5.00 | 3.32 | 0.50 | 0.000 | 137.79 |  |
| t00057 | 0.545 | 5095 | 464 | 20 | 22 | 3.30 | 5.00 | 15.00 | 0.50 | 0.000 | 184.92 |  |
| t00058 | 0.555 | 3775 | 532 | 20 | 22 | 2.00 | 11.14 | 4.14 | 5.44 | 1.00 | 200.00 |  |
| t00059 | 0.555 | 3538 | 525 | 20 | 22 | 2.00 | 10.63 | 4.95 | 5.34 | 0.16 | 200.00 |  |
| t00060 | 0.565 | 4634 | 463 | 20 | 22 | 3.27 | 5.00 | 15.00 | 0.50 | 1.00 | 200.00 |  |
| t00061 | 0.550 | 4443 | 470 | 20 | 22 | 4.36 | 5.00 | 14.51 | 0.50 | 0.000 | 50.36 |  |
| t00062 | 0.570 | 5331 | 456 | 20 | 22 | 2.00 | 5.00 | 10.67 | 10.00 | 0.000 | 0.002 |  |
| t00063 | 0.560 | 4346 | 465 | 20 | 22 | 3.60 | 5.00 | 12.89 | 0.50 | 1.00 | 200.00 |  |
| t00064 | 0.530 | 5668 | 456 | 20 | 22 | 2.00 | 5.00 | 8.84 | 0.50 | 0.000 | 0.002 |  |
| t00065 | 0.525 | 3413 | 1173 | 92 | 16 | 4.05 | 5.00 | 1.31 | 10.00 | 0.44 | 29.42 |  |
| t00066 | 0.495 | 5365 | 456 | 20 | 22 | 2.00 | 5.00 | 10.25 | 10.00 | 0.25 | 200.00 |  |
| t00067 | 0.570 | 2547 | 1216 | 102 | 14 | 2.00 | 7.91 | 1.00 | 8.41 | 0.70 | 33.18 |  |
| t00068 | 0.535 | 4923 | 421 | 20 | 19 | 2.00 | 5.00 | 11.90 | 10.00 | 0.000 | 30.74 | ★ |
| t00069 | 0.545 | 5150 | 607 | 35 | 22 | 2.00 | 5.00 | 15.00 | 6.02 | 0.000 | 150.56 |  |
| t00070 | 0.550 | 4886 | 456 | 20 | 22 | 2.00 | 5.00 | 12.09 | 10.00 | 0.24 | 0.002 |  |
| t00071 | 0.555 | 3363 | 1044 | 85 | 16 | 2.00 | 5.00 | 2.28 | 5.49 | 0.60 | 171.25 |  |
| t00072 | 0.575 | 3722 | 483 | 20 | 19 | 2.00 | 10.84 | 15.00 | 0.50 | 0.04 | 200.00 |  |
| t00073 | 0.550 | 4377 | 464 | 20 | 22 | 3.33 | 5.00 | 11.39 | 10.00 | 0.000 | 0.002 |  |
| t00074 | 0.585 | 3232 | 855 | 58 | 19 | 2.00 | 9.63 | 3.38 | 1.63 | 0.58 | 200.00 |  |
| t00075 | 0.540 | 5706 | 456 | 20 | 22 | 2.00 | 5.00 | 15.00 | 0.50 | 0.80 | 200.00 |  |
| t00076 | 0.530 | 2605 | 1129 | 97 | 12 | 2.00 | 6.38 | 7.19 | 6.65 | 0.91 | 63.91 |  |
| t00077 | 0.570 | 2882 | 893 | 65 | 17 | 2.00 | 9.16 | 4.46 | 0.65 | 0.46 | 200.00 |  |
| t00078 | 0.570 | 3925 | 778 | 48 | 21 | 2.00 | 9.31 | 15.00 | 0.50 | 0.48 | 171.59 |  |
| t00079 | 0.540 | 3406 | 522 | 20 | 20 | 2.00 | 12.97 | 2.81 | 0.50 | 0.31 | 127.60 |  |
| t00080 | 0.575 | 2977 | 487 | 20 | 17 | 2.00 | 14.28 | 3.43 | 0.50 | 0.000 | 70.72 |  |
| t00081 | 0.565 | 3711 | 469 | 20 | 19 | 2.84 | 9.00 | 14.45 | 0.50 | 1.00 | 139.29 |  |
| t00082 | 0.525 | 4496 | 466 | 20 | 22 | 3.62 | 5.00 | 15.00 | 8.38 | 1.00 | 200.00 |  |
| t00083 | 0.555 | 3245 | 1339 | 116 | 14 | 2.00 | 5.43 | 5.52 | 9.68 | 0.70 | 53.61 |  |
| t00084 | 0.560 | 3188 | 1136 | 89 | 17 | 2.00 | 9.19 | 3.00 | 6.64 | 0.66 | 200.00 |  |
| t00085 | 0.555 | 4250 | 843 | 57 | 22 | 2.00 | 6.02 | 15.00 | 0.50 | 0.57 | 200.00 |  |
| t00086 | 0.575 | 3463 | 934 | 63 | 22 | 2.00 | 8.51 | 15.00 | 3.86 | 0.57 | 75.02 |  |
| t00087 | 0.565 | 4721 | 719 | 46 | 22 | 2.00 | 5.00 | 15.00 | 0.50 | 0.79 | 200.00 |  |
| t00088 | 0.585 | 4112 | 418 | 20 | 17 | 2.00 | 7.11 | 4.40 | 2.13 | 0.000 | 200.00 | ★ |
| t00089 | 0.530 | 3416 | 442 | 20 | 16 | 2.00 | 11.19 | 4.38 | 2.63 | 0.000 | 200.00 |  |
| t00090 | 0.540 | 4873 | 410 | 20 | 18 | 2.00 | 5.00 | 5.78 | 0.50 | 0.000 | 200.00 | ★ |
| t00091 | 0.550 | 4724 | 410 | 20 | 18 | 2.00 | 5.00 | 3.47 | 0.50 | 0.000 | 186.94 | ★ |
| t00092 | 0.555 | 4530 | 447 | 20 | 19 | 2.00 | 7.43 | 3.23 | 2.46 | 0.40 | 200.00 | ★ |
| t00093 | 0.530 | 4922 | 490 | 20 | 22 | 2.00 | 7.75 | 15.00 | 0.50 | 0.80 | 200.00 |  |
| t00094 | 0.530 | 4282 | 510 | 20 | 22 | 2.00 | 9.37 | 15.00 | 6.80 | 0.44 | 50.37 |  |
| t00095 | 0.565 | 3893 | 525 | 20 | 22 | 2.00 | 10.60 | 13.60 | 0.50 | 0.18 | 106.94 |  |
| t00096 | 0.555 | 2894 | 686 | 45 | 16 | 2.00 | 10.18 | 3.33 | 3.24 | 0.63 | 134.78 |  |
| t00097 | 0.565 | 4154 | 664 | 40 | 22 | 2.51 | 5.00 | 1.00 | 0.97 | 1.00 | 200.00 |  |
| t00098 | 0.555 | 3583 | 1483 | 120 | 17 | 3.97 | 5.00 | 7.94 | 10.00 | 0.76 | 90.17 |  |
| t00099 | 0.545 | 4905 | 739 | 46 | 22 | 3.41 | 5.00 | 15.00 | 0.50 | 0.46 | 200.00 |  |

**qLogNEHVI** — 80 completados: eff_track med 0.550 (max 0.595), nhits_sipad med 4202 (max 6164), cost med 485 (min 271); en el frente: 15

---

## Anexo — eventos de ejemplo para el display

Índices de evento (`window_id`) del `ShipHits.root` de producción del 2026-09-11 16:43
(`gaudi_jobs/genie_pipeline/ShipHits.root`). Hits = SiTarget / SiPad / MTC SciFi.

**Global spliced a través de los tres detectores** (tipo 7; 60 eventos en total, todos con
una sola global):

| evento | hits SiT/SiP/MTC |
|---|---|
| 103 | 406 / 92 / 61 |
| 110 | 45 / 32 / 95 |
| 119 | 397 / 28 / 43 |
| 127 | 169 / 27 / 97 |
| 191 | 535 / 26 / 94 |
| 302 | 188 / 29 / 68 |
| 327 | 312 / 21 / 105 |
| 364 | 153 / 25 / 100 |

**Interacción solo en SiTarget** (global promovida de SiTarget, nada aguas abajo; 85
eventos):

| evento | globales | hits SiT/SiP/MTC |
|---|---|---|
| 9 | 1 | 282 / 0 / 0 |
| 289 | 1 | 84 / 0 / 0 |
| 242 | 1 | 565 / 0 / 0 |
| 274 | 1 | 907 / 0 / 0 |
| 343 | 1 | 1218 / 0 / 0 |
| 92 | 2 | 97 / 0 / 0 |

**Interacción solo en SiPad** (globales promovidas de SiPad; 153 eventos):

| evento | globales | hits SiT/SiP/MTC |
|---|---|---|
| 36 | 2 | 0 / 128 / 0 |
| 159 | 2 | 0 / 79 / 0 |
| 20 | 1 | 0 / 16 / 11 |
| 179 | 1 | 0 / 18 / 72 |
| 6 | 2 | 0 / 26 / 2 |
| 124 | 2 | 0 / 37 / 5 |

Splices de dos detectores (SiPad+MTC, tipo 6, o SiTarget+SiPad, tipo 3): 136 eventos;
ejemplos 12, 23, 35, 40 (SiPad+MTC), 25 (SiTarget+SiPad), 15 (uno de cada).
